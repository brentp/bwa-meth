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
    p.add_argument("--reference", help="reference fasta")
    p.add_argument("-t", "--threads", type=int, default=6)
    p.add_argument("--map-q", type=int, default=10, help="only tabulate "
                   "methylation for reads with at least this mapping quality")
    p.add_argument("bams", nargs="+")

    a = p.parse_args(args)
    assert os.path.exists(a.reference)

    cmd = "|samtools mpileup -f {reference} -r chr12 -d100000 -BIQ 20 -q {map_q} {bams}"
    cmd = cmd.format(reference=a.reference, map_q=a.map_q, bams=" ".join(a.bams))
    print >>sys.stderr, "generating pileup with command:", cmd
    samples = [op.basename(b)[:-4] for b in a.bams]
    tabulate_methylation(cmd, a.reference, samples)

def call_single_pileup(chrom, pos1, ref, coverage, bases, quals):
    """
    With directional protocol, we know:
        1. methylable site has 'G' on the opposite strand
        2. methylated site has 'C' on the + strand
        3. un-methylated site has 'T' on the + strand
        4. SNP has base(minus) != complement(reference)

    Return value is:
        {'p_snp': 0.5,
        'p_methylable': 0.5,
        'p_methylated|methylable': 0.5,
        'C': 20,
        'T': 20}

    p_methylable = n(G-) / n(total-)
    p_methylated|methylable = n(C+) / n(C+ + T+)
    p_snp = n(- != comp(ref)) / n(total -)

    """
    from collections import Counter
    print bases

    ref = ref.upper()
    plus  = Counter(b for b in bases if b != "," and b == b.upper())
    minus = Counter(b for b in bases if b != "." and b == b.lower())

    print plus
    print minus
    1/0


# TODO: use samtools snp calling with strand-filter off. 
# samtools mpileup -f /home/brentp/chr11.mm10.fa -d100000 -ugEIQ 20 -q 10
# bwa-meth.bam | bcftools view -Acvgm 0.99 -p 0.8 - | vcfutils.pl  varFilter -1 0

def tabulate_methylation(fpileup, reference, samples):

    conversion = {"C": "T", "G": "A"}
    fhs = []
    for sample in samples:
        fhs.append(open('%s.methylation.txt' % sample, 'w'))
        fhs[-1].write("#chrom\tpos0\tpos1\tpct\tn_same\tn_converted")

    for toks in (l.rstrip("\r\n").split("\t") for l in nopen(fpileup)):
        chrom, pos1, ref = toks[:3]
        pos1 = int(pos1)
        pos0 = pos1 - 1
        ref = ref.upper()
        converted = conversion.get(ref)
        if converted is None:
            continue
        s = faseq(reference, chrom, pos1 - 2, pos1 + 2)
        ctx = get_context(s, ref == "C")
        if not ctx.startswith("CG"): continue
        for i, sample in enumerate(samples):
            coverage, bases, quals = toks[3 + (i * 3): 6 + (i * 3)]
            if coverage == '0': continue

            # . == same on + strand, , == same on - strand
            n_same_plus = sum(1 for b in bases if b in ".")
            n_same_minus = sum(1 for b in bases if b in ",")
            n_same = n_same_plus + n_same_minus

            n_converted_plus = sum(1 for b in bases if b == converted)
            n_converted_minus = sum(1 for b in bases if b == converted.lower())
            n_converted = n_converted_plus + n_converted_minus
            #n_converted = sum(1 for b in bases if b == converted)
            # SNP
            n_other = sum(1 for b in bases.lower() if b in "actg" and b !=
                    converted.lower())
            if n_other > n_converted: continue

            if n_same + n_converted == 0: continue
            pct = 100. * n_same / float(n_same + n_converted)
            #fhs[i].write("{chrom}\t{pos0}\t{pos1}\t{pct:.1f}\t{n_same}\t{n_converted}\t{ctx}\n".format(**locals()))
            fhs[i].write("{chrom}\t{pos0}\t{pos1}\t{pct:.1f}\t{n_same}\t{n_converted}\n".format(**locals()))


if __name__ == "__main__":
    tabulate_main(sys.argv[1:])
