import sys
import os
import os.path as op
import argparse
from toolshed import nopen

def faseq(fa, chrom, start, end, cache=[None]):
    """
    this is called by pileup which is ordered by chrom
    so we can speed things up by reading in a chrom at
    a time into memory
    """
    if cache[0] is None or cache[0][0] != chrom:
        seq = "".join(x.strip() for i, x in
            enumerate(nopen("|samtools faidx %s %s" % (fa, chrom))) if i >
            0).upper()
        cache[0] = (chrom, seq)
    chrom, seq = cache[0]
    return seq[start - 1: end]

def get_context(seq5, forward):
    """
    >>> get_context('GACGG', True)
    'CG+'                  
    """
    if forward:
        assert seq5[2] == "C", seq5
        if seq5[3] == "G": return "CG+"
        if seq5[4] == "G": return "CHG+"
        return "CHH+"
    else: # reverse complement
        assert seq5[2] == "G", seq5
        if seq5[1] == "C": return "CG-"
        if seq5[0] == "C": return "CHG-"
        return "CHH-"

def tabulate_main(args):
    __doc__ = """
    tabulate methylation from bwa-meth.py call
    """
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--reference", help="reference fasta", required=True)
    p.add_argument("--region", help="specific region in which to run mpileup",
            default="")
    p.add_argument("--g-only", default=False, action="store_true",
        help="only use information from the reverse strand (G->A conversions)")
    p.add_argument("--map-q", "-q", type=int, default=10, help="only tabulate "
                   "methylation for reads with at least this mapping quality")
    p.add_argument("--base-q", "-Q", type=int, default=10, help="only tabulate "
                   "methylation for bases with at least this quality")
    p.add_argument("--skip-left", type=int, default=3, help="don't use the "
            "first N bases from the start of the read to avoid bias")
    p.add_argument("--skip-right", type=int, default=3, help="don't use the "
            "last N bases from the end of the read to avoid bias")
    p.add_argument("--read-length", type=int, help="length of reads"
            " before trimming for used with skipping", required=True)
    p.add_argument("--rf", action='append', default=[0x2], type=int)
    p.add_argument("--ff", action='append', default=[0x200, 0x400, 0x800],
            type=int)
    p.add_argument("bams", nargs="+")

    a = p.parse_args(args)
    assert os.path.exists(a.reference)
    if a.region:
        a.region = ("-l " if op.exists(a.region) else "-r ") + a.region
    ff = 0
    for flag in set(a.ff): ff |= flag
    rf = 0
    for flag in set(a.rf): rf |= flag
    cmd = ("|samtools mpileup --rf {rf} --ff {ff} -Of {reference} {region}"
           " -d100000 -BQ {base_q} -q {map_q} {bams}")

    cmd = cmd.format(reference=a.reference, map_q=a.map_q, base_q=a.base_q,
                     bams=" ".join(a.bams), region=a.region, rf=rf, ff=ff)
    sys.stderr.write("generating pileup with command: %s\n" % cmd)
    samples = [op.basename(b)[:-4] for b in a.bams]
    tabulate_methylation(cmd, a.reference, samples, a.g_only, a.skip_left,
            a.read_length - a.skip_right)

def tabulate_methylation(fpileup, reference, samples, g_only=False,
        skip_left=1, skip_right=99):

    conversion = {"C": "T", "G": "A"}
    fhs = []
    for sample in samples:
        fhs.append(open('%s.methylation.txt' % sample, 'w'))
        fhs[-1].write("#chrom\tpos0\tpos1\tpct\tn_same\tn_converted\n")

    import re
    # regexp to remove the ^, $ stuff
    end_re = re.compile("\^.|\$")
    deletion_re = re.compile("(-\d+[ACTGNatcgn]+)")

    insertion_re = re.compile("\+(\d+)([ACTGNatcgn]+)")
    deletion_re = re.compile("-(\d+)([ACTGNatcgn]+)")

    def del_sub(m):
        """ with a base string like "-1ACTG" just return CTG
        for -3CCCTG return TG
        """
        n = int(m.groups()[0])
        bases = m.groups()[1]
        return bases[n:]

    for toks in (l.rstrip("\r\n").split("\t") for l in nopen(fpileup)):
        chrom, pos1, ref = toks[:3]
        pos1 = int(pos1)
        pos0 = pos1 - 1
        ref = ref.upper()
        if g_only and ref != "G": continue
        converted = conversion.get(ref)
        if converted is None:
            continue
        s = faseq(reference, chrom, pos1 - 2, pos1 + 2)
        ctx = get_context(s, ref == "C")
        if not ctx.startswith("CG"): continue
        for i, sample in enumerate(samples):
            coverage, bases, quals, posns = toks[3 + (i * 4): 7 + (i * 4)]
            if coverage == '0': continue
            obases = bases[:]
            oquals = quals[:]

            bases = end_re.sub("", bases) # remove $ and ^.
            ebases = bases
            bases = deletion_re.sub(del_sub, bases)
            dbases = bases
            bases = insertion_re.sub(del_sub, bases)

            keep = [skip_left < int(p) < skip_right for p in posns.split(",")]
            if not any(keep): continue
            bases = "".join(b for j, b in enumerate(bases) if keep[j])
            if bases == "": continue

            # . == same on + strand, , == same on - strand
            n_same_plus = bases.count(".")
            n_same_minus = bases.count(",")
            n_same = n_same_plus + n_same_minus

            #n_converted_plus = b.count(converted)
            #n_converted_minus = b.count(converted_lower())
            #n_converted = n_converted_plus + n_converted_minus
            n_converted = bases.upper().count(converted)
            if n_same + n_converted == 0: continue
            # SNP
            n_other = sum(1 for b in bases.lower() if b in "actg" and b !=
                    converted.lower())
            # weak filtering here for a likely snp.
            if n_other > n_converted: continue

            pct = 100. * n_same / float(n_same + n_converted)
            #fhs[i].write("{chrom}\t{pos0}\t{pos1}\t{pct:.1f}\t{n_same}\t{n_converted}\t{ctx}\n".format(**locals()))
            fhs[i].write("{chrom}\t{pos0}\t{pos1}\t{pct:.1f}\t{n_same}\t{n_converted}\n".format(**locals()))

if __name__ == "__main__":
    tabulate_main(sys.argv[1:])
