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
import argparse

from itertools import groupby, izip
from toolshed import nopen, reader, is_newer_b
import string

class BWAMethException(Exception): pass

def comp(s, _comp=string.maketrans('ATCG', 'TAGC')):
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
        yield header, "".join(s.strip() for s in faiter.next()).upper()

def convert_reads(fq1, fq2, out=sys.stdout):
    print >>sys.stderr, "converting reads in %s,%s" % (fq1, fq2)
    fq1, fq2 = nopen(fq1), nopen(fq2)
    q1_iter = izip(*[fq1] * 4)
    q2_iter = izip(*[fq2] * 4)
    for pair in izip(q1_iter, q2_iter):
        for read_i, (name, seq, _, qual) in enumerate(pair):
            seq = seq.upper().rstrip('\n')
            char_a, char_b = ['CT', 'GA'][read_i]
            # keep original sequence as name.
            name = " ".join((name.split(" ")[0],
                            "YS:Z:" + seq +
                            "\tYC:Z:" + char_a + char_b + '\n'))
            out.write("".join((name, seq.replace(char_a, char_b) , "\n+\n", qual)))

def convert_fasta(ref_fasta, just_name=False):
    out_fa = op.splitext(ref_fasta)[0] + ".c2t.fa"
    if just_name:
        return out_fa
    msg = "c2t in %s to %s" % (ref_fasta, out_fa)
    if is_newer_b(ref_fasta, out_fa):
        print >>sys.stderr, "already converted", msg
        return out_fa
    print >>sys.stderr, "converting", msg
    try:
        fh = open(out_fa, "w")
        for header, seq in fasta_iter(ref_fasta):
            print >>fh, ">r%s" % header
            for line in wrap(seq.replace("G", "A")):
                print >>fh, line

            print >>fh, ">f%s" % header
            for line in wrap(seq.replace("C", "T")):
                print >>fh, line
        fh.close()
    except:
        fh.close(); os.unlink(out_fa)
        raise
    return out_fa


def bwa_index(fa):
    if is_newer_b(fa, (fa + '.amb', fa + '.sa')):
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

    def cig_len(self):
        return sum(c[0] for c in self.cigs() if c[1] in
                   ("M", "D", "N", "EQ", "X", "P"))

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


def bwa_mem(fa, mfq, extra_args, prefix='bwa-meth', threads=1, rg=None,
            calmd=False):
    conv_fa = convert_fasta(fa, just_name=True)
    if not is_newer_b(conv_fa, (conv_fa + '.amb', conv_fa + '.sa')):
        raise BWAMethException("first run bwa-meth.py index %s" % fa)

    if not rg is None and not rg.startswith('RG'):
        rg = '@RG\tID:{rg}\tSM:{rg}'.format(rg=rg)

    cmd = ("|bwa mem -L 25 -pCMR '{rg}' -t {threads} {extra_args} "
           "{conv_fa} {mfq}").format(**locals())
    print >>sys.stderr, "running: %s" % cmd.lstrip("|")
    as_bam(cmd, fa, prefix, calmd)


def as_bam(pfile, fa, prefix, calmd=False):
    """
    pfile: either a file or a |process to generate sam output
    fa: the reference fasta
    prefix: the output prefix or directory
    """
    calmd = (" samtools sort -m3G -@3 -o {bam}.fix.bam {bam}.tmp"
             " | samtools calmd -AbEr - {fa} > {bam}.bam "
             ).format(fa=fa, bam=prefix) \
             if calmd \
             else "samtools sort -m3G -@3 {bam}.fix.bam {bam} ".format(bam=prefix)

    cmd1 = ("set -eo pipefail; samtools view -bS - "
           "| samtools sort -nm 3G -@3 - {bam} ").format(bam=prefix)

    cmd2 = ("samtools fixmate {bam}.bam {bam}.fix.bam; "
           "{calmd} "
           "&& samtools index {bam}.bam "
           "&& rm -f {bam}.fix.bam").format(bam=prefix, calmd=calmd)
    print >>sys.stderr, "writing to:\n", cmd1.replace("&&", "\\\n\t&&")\
                                            .replace("|", "\\\n\t|")
    p = nopen("|" + cmd1, 'w')
    out = p.stdin
    PG = True
    lengths = {}
    for toks in reader("%s" % (pfile, ), header=False):
        if toks[0].startswith("@"):
            if toks[0].startswith("@SQ"):
                sq, sn, ln = toks  # @SQ    SN:fchr11    LN:122082543
                # we have f and r, only print out f
                chrom = sn.split(":")[1]
                if chrom.startswith('r'): continue
                chrom = chrom[1:]
                lengths[chrom] = int(ln.split(":")[1])
                toks = ["%s\tSN:%s\t%s" % (sq, chrom, ln)]
            if toks[0].startswith("@PG"): continue
            out.write("\t".join(toks) + "\n")
            continue
        if PG:
            #print >>out, "@PG\tprog:bwa-meth.py"
            PG = False

        aln = Bam(toks)

        orig_seq = aln.original_seq
        # don't need this any more.
        aln.other = [x for x in aln.other if not x.startswith('YS:Z')]
        if aln.chrom == "*":  # chrom
            print >>out, str(aln)
            continue

        # first letter of chrom is 'f' or 'r'
        direction = aln.chrom[0]
        aln.chrom = aln.chrom.lstrip('fr')

        if not aln.is_mapped():
            aln.seq = orig_seq
            print >>out, str(aln)
            continue
        assert direction in 'fr', (direction, toks[2], aln)
        aln.other.append('YD:Z:' + direction)

        if aln.chrom_mate[0] not in "*=":
            mate_direction = aln.chrom_mate[0]
            aln.chrom_mate = aln.chrom_mate[1:]
            if mate_direction == 'r':
                aln.flag ^= 0x20

        # adjust the original seq to the cigar
        l, r = aln.left_shift(), aln.right_shift()
        if aln.is_plus_read():
            aln.seq = orig_seq[l:r]
        else:
            aln.seq = comp(orig_seq[::-1][l:r])
        if direction == 'r':
            aln.flag ^= 0x10
            aln.flag ^= 0x20

            aln.seq = comp(aln.seq[::-1])
            aln.qual = aln.qual[::-1]
            aln.pos = lengths[aln.chrom] - aln.pos - aln.cig_len() + 2
            aln.cigar = "".join(["%s%s" % c for c in aln.cigs()][::-1])

            if aln.chrom_mate != "*":
                aln.pos_mate = int(aln.pos_mate)
                # this only works when the other read has same cigar as this one.
                # so we just send to samtools fixmate.
                aln.pos_mate = str(lengths[aln.chrom] - aln.pos_mate + 2 - aln.cig_len())

        print >>out, str(aln)
    out.close()
    print >>sys.stderr, "running:\n", cmd2.replace("&&", "\\\n\t&&")\
                                            .replace("|", "\\\n\t|")
    run(cmd2)

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

def tabulate_main(args):
    __doc__ = """
    tabulate methylation from bwa-meth.py call
    """
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--reference", help="reference fasta")
    p.add_argument("-t", "--threads", type=int, default=6)
    p.add_argument("--dbsnp", help="optional dbsnp for GATK calibration")
    p.add_argument("--prefix", help="output prefix", default='bmeth-tab')
    p.add_argument("--trim", help="left, right trim to avoid bias",
                        default="2,2")
    p.add_argument("--map-q", type=int, default=10, help="only tabulate "
                   "methylation for reads with at least this mapping quality")
    p.add_argument("--bissnp", help="path to bissnp jar")
    p.add_argument("bams", nargs="+")

    a = p.parse_args(args)
    assert os.path.exists(a.reference)
    trim = map(int, a.trim.split(","))

    cmd = """\
    java -Xmx15g -jar {bissnp}
        -R {reference}
        -I {bams}
        -T BisulfiteGenotyper
         -stand_emit_conf 0
        --trim_5_end_bp {trim5}
        --trim_3_end_bp {trim3}
        -vfn1 {prefix}.cpg.vcf -vfn2 {prefix}.snp.vcf
        -stand_call_conf 10
        -mmq {mapq} {dbsnp}
        -nt {threads}""".format(
            threads=a.threads,
            dbsnp=("--dbsnp " + a.dbsnp) if a.dbsnp else "",
            bissnp=a.bissnp,
            trim5=trim[0],
            trim3=trim[1],
            prefix=a.prefix,
            reference=a.reference,
            mapq=a.map_q,
            bams=" -I ".join(a.bams)).replace("\n", " \\\n")


    print cmd


def main(args):

    if len(args) > 0 and args[0] == "index":
        sys.exit(bwa_index(convert_fasta(args[1])))

    if len(args) > 0 and args[0] == "c2t":
        sys.exit(convert_reads(args[1], args[2]))

    if len(args) > 0 and args[0] == "tabulate":
        sys.exit(tabulate_main(args[1:]))

    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--reference", help="reference fasta")
    p.add_argument("-t", "--threads", type=int, default=6)
    p.add_argument("-p", "--prefix", default="bwa-meth")
    p.add_argument("--calmd", default=False, action="store_true")
    p.add_argument("--read-group", help="read-group to add to bam in same"
            " format as to bwa: '@RG\\tID:foo\\tSM:bar'")
    p.add_argument("fastqs", nargs="+", help="bs-seq fastqs to align")

    args = p.parse_args(args)
    # for the 2nd file. use G => A and bwa's support for streaming.
    conv_fqs = "'<python %s c2t %s %s'" % (__file__, args.fastqs[0],
                                                      args.fastqs[1])

    bwa_mem(args.reference, conv_fqs, "", prefix=args.prefix,
             threads=args.threads, rg=args.read_group or
             rname(*args.fastqs), calmd=args.calmd)

if __name__ == "__main__":
    main(sys.argv[1:])
