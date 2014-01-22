#!/usr/bin/env python
"""
map bisulfite converted reads to an insilico converted genome using bwa mem.
A command to this program like:

    python bwameth.py --reference ref.fa A.fq B.fq

Gets converted to:

    bwa mem -pCMR ref.c2t.fa '<python bwameth.py c2t A.fq B.fq'

So that A.fq has C's converted to T's and B.fq has G's converted to A's
and both are streamed directly to the aligner without a temporary file.
The output is a corrected, sorted, indexed BAM.
"""

import sys
import os
import os.path as op
import argparse
from subprocess import check_call

try:
    from itertools import groupby, izip
except ImportError: # python3
    izip = zip
from toolshed import nopen, reader, is_newer_b
import string

__version__ =  "0.06"

def checkX(cmd):
    for p in os.environ['PATH'].split(":"):
        if os.access(os.path.join(p, cmd), os.X_OK):
            break
    else:
        raise Exception("executable for '%s' not found" % cmd)

checkX('samtools')
checkX('bwa')

class BWAMethException(Exception): pass

def comp(s, _comp=string.maketrans('ATCG', 'TAGC')):
    return s.translate(_comp)

def wrap(text, width=100): # much faster than textwrap
    for s in xrange(0, len(text), width):
        yield text[s:s+width]

def run(cmd):
    list(nopen("|%s" % cmd.lstrip("|")))

def fasta_iter(fasta_name):
    fh = nopen(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        yield header, "".join(s.strip() for s in faiter.next()).upper()

def convert_reads(fq1, fq2, out=sys.stdout):
    sys.stderr.write("converting reads in %s,%s\n" % (fq1, fq2))
    fq1, fq2 = nopen(fq1), nopen(fq2)
    q1_iter = izip(*[fq1] * 4)
    q2_iter = izip(*[fq2] * 4)

    jj = 0
    for pair in izip(q1_iter, q2_iter):

        for read_i, (name, seq, _, qual) in enumerate(pair):
            name = name.rstrip("\r\n").split(" ")[0]
            if name.endswith(("_R1", "_R2")):
                name = name[:-3]
            elif name.endswith(("/1", "/2")):
                name = name[:-2]

            seq = seq.upper().rstrip('\n')
            char_a, char_b = ['CT', 'GA'][read_i]
            # keep original sequence as name.
            name = " ".join((name,
                            "YS:Z:" + seq +
                            "\tYC:Z:" + char_a + char_b + '\n'))
            seq = seq.replace(char_a, char_b)
            out.write("".join((name, seq, "\n+\n", qual)))

        jj += 1
        #if jj > 190000: break
    out.flush()
    out.close()

def convert_fasta(ref_fasta, just_name=False):
    out_fa = op.splitext(ref_fasta)[0] + ".c2t.fa"
    if just_name:
        return out_fa
    msg = "c2t in %s to %s" % (ref_fasta, out_fa)
    if is_newer_b(ref_fasta, out_fa):
        sys.stderr.write("already converted: %s\n" % msg)
        return out_fa
    sys.stderr.write("converting %s\n" % msg)
    try:
        fh = open(out_fa, "w")
        for header, seq in fasta_iter(ref_fasta):
            ########### Reverse ######################
            fh.write(">r%s\n" % header)

            #if non_cpg_only:
            #    for ctx in "TAG": # use "ATC" for fwd
            #        seq = seq.replace('G' + ctx, "A" + ctx)
            #    for line in wrap(seq):
            #        print >>fh, line
            #else:
            for line in wrap(seq.replace("G", "A")):
                fh.write(line + '\n')

            ########### Forward ######################
            fh.write(">f%s\n" % header)
            for line in wrap(seq.replace("C", "T")):
                fh.write(line + '\n')
        fh.close()
    except:
        fh.close(); os.unlink(out_fa)
        raise
    return out_fa


def bwa_index(fa):
    if is_newer_b(fa, (fa + '.amb', fa + '.sa')):
        return
    sys.stderr.write("indexing: %s\n" % fa)
    try:
        run("bwa index %s" % fa)
    except:
        if op.exists(fa + ".amb"):
            os.unlink(fa + ".amb")
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
        raise BWAMethException("first run bwameth.py index %s" % fa)

    if not rg is None and not rg.startswith('RG'):
        rg = '@RG\tID:{rg}\tSM:{rg}'.format(rg=rg)

    # penalize clipping and unpaired. lower penalty on mismatches (-B)
    cmd = ("|bwa mem -T 40 -B 3 -L 25 -U 100 -pCMR '{rg}' -t {threads} {extra_args} "
           "{conv_fa} {mfq}").format(**locals())
    sys.stderr.write("running: %s\n" % cmd.lstrip("|"))
    as_bam(cmd, fa, prefix, calmd)


def as_bam(pfile, fa, prefix, calmd=False):
    """
    pfile: either a file or a |process to generate sam output
    fa: the reference fasta
    prefix: the output prefix or directory
    """
    view = "samtools view -bS - | samtools sort -@3 - "
    if calmd:
        cmds = [
            view + "{bam}.tmp",
            "samtools calmd -AbEr {bam}.tmp.bam {fa} > {bam}.bam 2>/dev/null",
            "rm {bam}.tmp.bam"]
    else:
        cmds = [view + "{bam}"]

    cmds.append("samtools index {bam}.bam")
    cmds = [c.format(bam=prefix, fa=fa) for c in cmds]

    sys.stderr.write("writing to:\n%s\n" % cmds[0])

    p = nopen("|" + cmds[0], 'w')
    out = p.stdin

    for toks in reader("%s" % (pfile, ), header=False):
        if toks[0].startswith("@"):
            handle_header(toks, out)
            continue

        aln = handle_read(Bam(toks))
        out.write(str(aln) + '\n')

    p.stdin.flush()
    p.stdout.flush()
    p.communicate()
    out.close()
    for cmd in cmds[1:]:
        sys.stderr.write("running: %s\n" % cmd.strip())
        assert check_call(cmd.strip(), shell=True) == 0

def handle_header(toks, out):
    if toks[0].startswith("@SQ"):
        sq, sn, ln = toks  # @SQ    SN:fchr11    LN:122082543
        # we have f and r, only print out f
        chrom = sn.split(":")[1]
        if chrom.startswith('r'): return
        chrom = chrom[1:]
        toks = ["%s\tSN:%s\t%s" % (sq, chrom, ln)]
    if toks[0].startswith("@PG"):
        toks = ["@PG\tID:bwa-meth\tPN:bwa-meth\tVN:%s\tCL:%s" % (
                         __version__,
                         " ".join(x.replace("\t", "\\t") for x in sys.argv))]
    out.write("\t".join(toks) + "\n")


def handle_read(aln):

    orig_seq = aln.original_seq
    # don't need this any more.
    aln.other = [x for x in aln.other if not x.startswith('YS:Z')]
    if aln.chrom == "*":  # chrom
        return aln

    # first letter of chrom is 'f' or 'r'
    direction = aln.chrom[0]
    aln.chrom = aln.chrom.lstrip('fr')

    if not aln.is_mapped():
        aln.seq = orig_seq
        return aln

    assert direction in 'fr', (direction, aln)
    aln.other.append('YD:Z:' + direction)

    mate_direction = aln.chrom_mate[0]
    if mate_direction not in "*=":
        aln.chrom_mate = aln.chrom_mate[1:]

    # adjust the original seq to the cigar
    l, r = aln.left_shift(), aln.right_shift()
    if aln.is_plus_read():
        aln.seq = orig_seq[l:r]
    else:
        aln.seq = comp(orig_seq[::-1][l:r])

    return aln


def cnvs_main(args):
    __doc__ = """
    calculate CNVs from BS-Seq bams or vcfs
    """
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--regions", help="optional target regions", default='NA')
    p.add_argument("bams", nargs="+")

    a = p.parse_args(args)
    r_script = """
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(cn.mops))
suppressPackageStartupMessages(library(snow))
args = commandArgs(TRUE)
regions = args[1]
bams = args[2:length(args)]
n = length(bams)
if(is.na(regions)){
    bam_counts = getReadCountsFromBAM(bams, parallel=min(n, 4), mode="paired")
    res = cn.mops(bam_counts, parallel=min(n, 4), priorImpact=20)
} else {
    segments = read.delim(regions, header=FALSE)
    gr = GRanges(segments[,1], IRanges(segments[,2], segments[,3]))
    bam_counts = getSegmentReadCountsFromBAM(bams, GR=gr, mode="paired", parallel=min(n, 4))
    res = exomecn.mops(bam_counts, parallel=min(n, 4), priorImpact=20)
}
res = calcIntegerCopyNumbers(res)

df = as.data.frame(cnvs(res))
write.table(df, row.names=FALSE, quote=FALSE, sep="\t")
"""
    import tempfile
    with tempfile.NamedTemporaryFile(delete=True) as rfh:
        rfh.write(r_script + '\n')
        rfh.flush()
        for d in reader('|Rscript {rs_name} {regions} {bams}'.format(
            rs_name=rfh.name, regions=a.regions, bams=" ".join(a.bams)),
            header=False):
            print("\t".join(d))


def tabulate_main(args):
    __doc__ = """
    tabulate methylation from bwameth.py call
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
    assert os.path.exists(a.reference + ".fai"), ("run samtools faidx %s" \
                                                  % a.reference)
    trim = map(int, a.trim.split(","))

    cmd = """\
    java -Xmx15g -jar {bissnp}
        -R {reference}
        -I {bams}
        -T BisulfiteGenotyper
        --trim_5_end_bp {trim5}
        --trim_3_end_bp {trim3}
        -vfn1 {prefix}.cpg.vcf -vfn2 {prefix}.snp.vcf
        -mbq 20
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
    sys.stderr.write(cmd + '\n')

    run(cmd)
    fmt = "{CHROM}\t{START}\t{POS}\t{PCT}\t{CS}\t{TS}\t{CTX}\n"
    sys.stderr.write(a.prefix + ".cpg.vcf\n")
    for i, d in enumerate(reader(a.prefix + ".cpg.vcf",
                       skip_while=lambda toks: toks[0] != "#CHROM",
                       header="ordered")):
        if i == 0:
            samples = d.keys()[9:]
            fhs = {}
            for sample in samples:
                fhs[sample] = open("{prefix}{sample}.cpg.bed"\
                        .format(prefix=a.prefix, sample=sample), "w")

                fhs[sample].write("#" + fmt.replace("}", "").replace("{", ""))
        if not d['FILTER'] in (".", "PASS"): continue
        d['START'] = str(int(d['POS']) - 1)
        for sample in samples:
            info = dict(zip(d['FORMAT'].split(":"), d[sample].split(":")))
            try:
                d['CS'] = int(info['CM'])  # (M)ethylated
                d['TS'] = int(info['CU'])  # (U)n
            except ValueError:
                continue
            d['CTX'] = info['CP']
            if d['CS'] + d['TS'] == 0:
                continue
            else:
                d['PCT'] = "%i" % (100 * d['CS'] / float(d['CS'] + d['TS']))
            fhs[sample].write(fmt.format(**d))


def main(args=sys.argv[1:]):

    if len(args) > 0 and args[0] == "index":
        assert len(args) == 2, ("must specify fasta as 2nd argument")
        sys.exit(bwa_index(convert_fasta(args[1])))

    if len(args) > 0 and args[0] == "c2t":
        sys.exit(convert_reads(args[1], args[2]))

    if len(args) > 0 and args[0] == "tabulate":
        sys.exit(tabulate_main(args[1:]))

    if len(args) > 0 and args[0] == "cnvs":
        sys.exit(cnvs_main(args[1:]))

    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--reference", help="reference fasta", required=True)
    p.add_argument("-t", "--threads", type=int, default=6)
    p.add_argument("-p", "--prefix", default="bwa-meth")
    p.add_argument("--calmd", default=False, action="store_true")
    p.add_argument("--read-group", help="read-group to add to bam in same"
            " format as to bwa: '@RG\\tID:foo\\tSM:bar'")
    p.add_argument("fastqs", nargs="+", help="bs-seq fastqs to align")

    args = p.parse_args(args)
    # for the 2nd file. use G => A and bwa's support for streaming.
    script = __file__
    conv_fqs = "'<python %s c2t %s %s'" % (script, args.fastqs[0],
                                                   args.fastqs[1])

    bwa_mem(args.reference, conv_fqs, "", prefix=args.prefix,
             threads=args.threads, rg=args.read_group or
             rname(*args.fastqs), calmd=args.calmd)

if __name__ == "__main__":
    main(sys.argv[1:])
