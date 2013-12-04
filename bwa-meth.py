"""
map bisulfite converted reads to an insilico converted genome using bwa mem.
A command to this program like:

    python bwa-meth.py --reference ref.fa A.fq B.fq

Gets converted to:

    bwa mem ref.c2t.fa '<python bwa-meth.py c2t A.fq' '<python bwa-meth.py g2a B.fq'

So that the reference with C converted to T is created and indexed
automatically and no temporary files are written for the fastqs. The output is
a corrected, indexed BAM, and a BED file similar to that output by Bismark with
cs, ts, and percent methylation at each site.
"""

import sys
import os
import os.path as op

from itertools import groupby, izip
from toolshed import nopen, reader
import string

DEBUG = True

def comp(s, _comp=string.maketrans('ATCGNt', 'TAGCNa')):
    return s.translate(_comp)

def wrap(text, width=100): # much faster than textwrap
    for s in xrange(0, len(text), width):
        yield text[s:s+width]

def run(cmd):
    list(nopen("|%s" % cmd.rstrip("|")))

def fasta_iter(fasta_name):
    fh = nopen(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        yield header, "".join(s.strip() for s in faiter.next())

def convert_reads(fq, ga=False, out=sys.stdout, trim=(0, 0)):
    char_a, char_b = 'Ga' if ga else 'Ct'
    if trim == (0, 0): trim = None
    print >>sys.stderr, "converting %s to %s in %s" % (char_a, char_b, fq)
    fq = nopen(fq)
    for (name, seq, _, qual) in izip(*[fq] * 4):
        seq = seq.upper().rstrip('\n')
        if trim is not None:
            seq = seq[trim[0]:-trim[1]]
            qual = qual.rstrip('\n')[trim[0]:-trim[1]] + '\n'

        # keep original sequence as name.
        name = " ".join((name.split(" ")[0],
                        "YS:Z:" + seq +
                        "\tYC:Z:" + char_a + char_b + '\n'))
        out.write("".join((name, seq.replace(char_a, char_b) , "\n+\n", qual)))

def convert_fasta(ref_fasta):
    print >>sys.stderr, "converting c2t in %s" % ref_fasta
    out_fa = op.splitext(ref_fasta)[0] + ".c2t.fa"
    if op.exists(out_fa): return out_fa
    try:
        fh = open(out_fa, "w")
        for header, seq in fasta_iter(ref_fasta):
            seq = seq.upper()
            print >>fh, ">f%s" % header
            for line in wrap(seq.replace("C", "t")):
                print >>fh, line
            print >>fh, ">r%s" % header
            for line in wrap(comp(seq).replace("C", "t")):
                print >>fh, line
        fh.close()
    except:
        fh.close() or os.unlink(out_fa)
        raise
    return out_fa

def bwa_index(fa):
    if op.exists(fa + '.amb'):
        return
    print >>sys.stderr, "indexing: %s" % fa
    try:
        run("bwa index %s" % fa)
    except:
        if op.exists(fa + ".amb"):
            os.unlink(fa + ".bam")
        raise

class Bam(object):
    __slots__ = 'read flag chrom pos mapq cigar chrom_mate pos_mate tlen \
            seq qual other'.split()
    def __init__(self, args):
        for a, v in zip(self.__slots__[:11], args):
            setattr(self, a, v)
        self.other = args[11:]
        self.flag = int(self.flag)
        self.pos = int(self.pos)
        self.tlen = int(float(self.tlen))
        try:
            self.mapq = int(self.mapq)
        except ValueError:
            pass

    def __repr__(self):
        return "Bam({chr}:{start}:{read}".format(chr=self.chrom,
                                                 start=self.pos,
                                                 read=self.read)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__[:11]) \
                         + "\t" + "\t".join(self.other)

    def is_first_read(self):
        return bool(self.flag & 0x40)

    def is_second_read(self):
        return bool(self.flag & 0x80)

    def is_plus_read(self):
        return not (self.flag & 0x10)

    def is_minus_read(self):
        return bool(self.flag & 0x10)

    def is_mapped(self):
        return not (self.flag & 0x4)

    def cigs(self):
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(cig_iter.next()[1])

    def left_shift(self):
        left = 0
        for n, cig in self.cigs():
            if cig == "M": break
            if cig == "H":
                left += n
        return left

    def right_shift(self):
        right = 0
        for n, cig in reversed(list(self.cigs())):
            if cig == "M": break
            if cig == "H":
                right += n
        return -right or None

    @property
    def start(self):
        return self.pos

    @property
    def original_seq(self):
        return next(x for x in self.other if x.startswith("YS:Z:"))[5:]

    @property
    def ga_ct(self):
        return [x for x in self.other if x.startswith("YC:Z:")]

def rname(fq1, fq2):
    def name(f):
        n = op.basename(op.splitext(f)[0])
        if n.endswith('.fastq'): n = n[:-6]
        if n.endswith(('.fq', '.r1', '.r2')): n = n[:-3]
        return n
    return "".join(a for a, b in zip(name(fq1), name(fq2)) if a == b) or 'bm'

def bwa_meth(ref_fasta, fastqs, prefix=".", extra_args="", threads=1, mapq=0,
             rg=None):
    assert len(fastqs) in (1, 2)
    return bwa_mem(ref_fasta, fastqs, extra_args, prefix=prefix,
                   threads=threads, mapq=mapq, rg=rg)


def bwa_mem(fa, fqs, extra_args, prefix='bwa-meth', threads=1, mapq=0, rg=None):
    conv_fa = convert_fasta(fa)
    bwa_index(conv_fa)
    if not rg is None and not rg.startswith('RG'):
        rg = '@RG\tID:{rg}\tSM:{rg}'.format(rg=rg)

    fq_str = " ".join(fqs)

    cmd = ("|bwa mem -U 20 -L 10 -CMR '{rg}' -t {threads} {extra_args} "
           "{conv_fa} {fq_str}").format(**locals())
    print >>sys.stderr, "running: %s" % cmd.lstrip("|")
    tabulate(cmd, fa, prefix, conv_fa=conv_fa, mapq=mapq)


def tabulate(pfile, fa, prefix, mapq=0, debug=DEBUG, conv_fa=None):
    """
    pfile: either a file or a |process to generate sam output
    fa: the reference fasta
    prefix: the output prefix or directory
    mapq: only tabulate methylation for reads with at least this mapping
          quality
    debug: if True, the alignments and internal stuff to calculate c=>t are
           printed
    conv_fa: needed if debug is True
    """
    out = nopen(prefix + ".sam.gz", "w")
    PG = True
    for toks in reader("%s" % (pfile, ), header=False):
        if toks[0].startswith("@"):
            if toks[0].startswith("@SQ"):
                sq, sn, ln = toks
                # we have f and r, only print out f
                sn = sn.split(":")[1]
                if sn.startswith('r'): continue
                toks[1] = toks[1].replace(":f", ":")
            if toks[0].startswith("@PG"): continue
            print >>out, "\t".join(toks)
            continue
        if PG:
            #print >>out, "@PG\tprog:bwa-meth.py"
            PG = False

        aln = Bam(toks)
        orig_seq = aln.original_seq
        # don't need this any more.
        aln.other = [x for x in aln.other if not x.startswith('YS:Z')]
        if aln.chrom == "*": # chrom
            print >>out, str(aln)
            continue

        # first letter of chrom is 'f' or 'r'
        direction = aln.chrom[0]
        aln.chrom = aln.chrom.lstrip('fr')
        aln.other.append('YD:Z:' + direction)

        if not aln.is_mapped():
            aln.seq = orig_seq
            print >>out, str(aln)
            continue
        assert direction in 'fr', (direction, toks[2], aln)

        if aln.chrom_mate[0] not in "*=":
            aln.chrom_mate = aln.chrom_mate[1:]

        # adjust the original seq to the cigar
        l, r = aln.left_shift(), aln.right_shift()
        if aln.is_plus_read():
            aln.seq = orig_seq[l:r]
        else:
            aln.seq = comp(orig_seq)[::-1][l:r]
        #if direction == 'r':
        #    aln.flag ^= 0x10
        print >>out, str(aln)

    out.close()
    try:
        run("zless {sam} | samtools view -hbS - \
                | samtools sort -m 3G - {bam} \
                && samtools index {bam}.bam".format(sam=out.name, bam=prefix))
    except:
        if op.exists(prefix + ".bam"):
            os.unlink(prefix + ".bam")
            if op.exists(prefix + ".bam.bai"):
                os.unlink(prefix + ".bam.bai")
        raise

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

def summarize_pileup(fpileup, reference):

    conversion = {"C": "T", "G": "A"}
    print "#chrom\tpos1\tn_same\tn_converted\tcontext"

    for toks in (l.rstrip("\r\n").split("\t") for l in nopen(fpileup)):
        chrom, pos1, ref, coverage, bases, quals = toks
        #if int(coverage) < 4: continue
        pos1 = int(pos1)
        if coverage == '0': continue
        ref = ref.upper()
        converted = conversion.get(ref)
        if converted is None:
            continue

        s = faseq(reference, chrom, pos1 - 2, pos1 + 2)
        ctx = get_context(s, ref == "C")
        if not ctx.startswith("CG"): continue

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

        if n_same < 10 or n_converted < 10: continue
        pct = n_same / float(n_same + n_converted)
        print bases

        print "{chrom}\t{pos1}\t{pct}\t{n_same_plus}\t{n_same_minus}\t{n_converted_plus}\t{n_converted_minus}\t{ctx}\t{s}\t{ref}".format(**locals())


def main(args):

    #summarize_pileup('|samtools mpileup -f {reference} -BIQ 20 -q {map_q} {bams}'.format(
    #    bams="bwa-meth.bam", map_q=20, reference="~/chr11.mm10.fa"), "~/chr11.mm10.fa")
    #1/0

    if len(args) > 0 and args[0] in ("c2t", "g2a"):
        # catch these args to convert reads on the fly and stream to bwa
        trim = map(int, args[2].split(",")) if len(args) > 2 else (0, 0)
        sys.exit(convert_reads(args[1], ga=args[0] == "g2a", trim=trim))

    import argparse
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--reference", help="reference fasta")
    p.add_argument("-t", "--threads", type=int, default=6)
    p.add_argument("--trim", default="0,0", help="bases to trim from"
    "start and end of each read to avoid bias. '2,2' is recommended")
    p.add_argument("-p", "--prefix", default="bwa-meth")
    p.add_argument("--read-group", help="read-group to add to bam in same"
            " format as to bwa: '@RG\tID:foo\tSM:bar'")
    p.add_argument("--map-q", type=int, default=10, help="only tabulate "
                   "methylation for reads with at least this mapping quality")
    p.add_argument("fastqs", nargs="+", help="bs-seq fastqs to align")

    args = p.parse_args(args)
    # for the 2nd file. use G => A and bwa's support for streaming.
    conv_fqs = ["'<python %s %s %s %s'" % (__file__, conv, fq, args.trim) for conv, fq in
                                zip(('c2t', 'g2a'), args.fastqs)]
    bwa_meth(args.reference, conv_fqs, prefix=args.prefix,
            threads=args.threads, mapq=args.map_q, rg=args.read_group or
            rname(*args.fastqs))

if __name__ == "__main__":
    main(sys.argv[1:])
