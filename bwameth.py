#!/usr/bin/env python
"""
map bisulfite converted reads to an insilico converted genome using bwa mem.
A command to this program like:

    python bwameth.py --reference ref.fa A.fq B.fq

Gets converted to:

    bwa mem -pCMR ref.fa.bwameth.c2t '<python bwameth.py c2t A.fq B.fq'

So that A.fq has C's converted to T's and B.fq has G's converted to A's
and both are streamed directly to the aligner without a temporary file.
The output is a corrected, sorted, indexed BAM.
"""
from __future__ import print_function
import tempfile
import sys
import os
import os.path as op
import argparse
from subprocess import check_call
from operator import itemgetter
from itertools import groupby, repeat, chain
import re

try:
    from itertools import izip
    import string
    maketrans = string.maketrans
except ImportError: # python3
    izip = zip
    maketrans = str.maketrans
from toolshed import nopen, reader, is_newer_b

__version__ = "0.10"

def checkX(cmd):
    for p in os.environ['PATH'].split(":"):
        if os.access(os.path.join(p, cmd), os.X_OK):
            break
    else:
        raise Exception("executable for '%s' not found" % cmd)

checkX('samtools')
checkX('bwa')

class BWAMethException(Exception): pass

def comp(s, _comp=maketrans('ATCG', 'TAGC')):
    return s.translate(_comp)

def wrap(text, width=100): # much faster than textwrap
    try: xrange
    except NameError: xrange = range
    for s in xrange(0, len(text), width):
        yield text[s:s+width]

def run(cmd):
    list(nopen("|%s" % cmd.lstrip("|")))

def fasta_iter(fasta_name):
    fh = nopen(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        yield header, "".join(s.strip() for s in next(faiter)).upper()

def convert_reads(fq1s, fq2s, out=sys.stdout):

    for fq1, fq2 in zip(fq1s.split(","), fq2s.split(",")):
        sys.stderr.write("converting reads in %s,%s\n" % (fq1, fq2))
        fq1 = nopen(fq1)
        if fq2 != "NA":
            fq2 = nopen(fq2)
            q2_iter = izip(*[fq2] * 4)
        else:
            sys.stderr.write("WARNING: running bwameth in single-end mode\n")
            q2_iter = repeat((None, None, None, None))
        q1_iter = izip(*[fq1] * 4)

        lt80 = 0
        for pair in izip(q1_iter, q2_iter):
            for read_i, (name, seq, _, qual) in enumerate(pair):
                if name is None: continue
                name = name.rstrip("\r\n").split(" ")[0]
                if name[0] != "@":
                    sys.stderr.write("""ERROR!!!!
ERROR!!! FASTQ conversion failed
ERROR!!! expecting FASTQ 4-tuples, but found a record %s that doesn't start with "@"
""" % name)
                    sys.exit(1)
                if name.endswith(("_R1", "_R2")):
                    name = name[:-3]
                elif name.endswith(("/1", "/2")):
                    name = name[:-2]

                seq = seq.upper().rstrip('\n')
                if len(seq) < 80:
                    lt80 += 1

                char_a, char_b = ['CT', 'GA'][read_i]
                # keep original sequence as name.
                name = " ".join((name,
                                "YS:Z:" + seq +
                                "\tYC:Z:" + char_a + char_b + '\n'))
                seq = seq.replace(char_a, char_b)
                out.write("".join((name, seq, "\n+\n", qual)))

    out.flush()
    if lt80 > 50:
        sys.stderr.write("WARNING: %i reads with length < 80\n" % lt80)
        sys.stderr.write("       : this program is designed for long reads\n")
    return 0

def convert_fasta(ref_fasta, just_name=False):
    out_fa = ref_fasta + ".bwameth.c2t"
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
        try:
            fh.close()
        except UnboundLocalError:
            pass
        os.unlink(out_fa)
        raise
    return out_fa


def bwa_index(fa):
    if is_newer_b(fa, (fa + '.amb', fa + '.sa')):
        return
    sys.stderr.write("indexing: %s\n" % fa)
    try:
        run("bwa index -a bwtsw %s" % fa)
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
            yield int("".join(n)), "".join(next(cig_iter)[1])

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
        try:
            return next(x for x in self.other if x.startswith("YS:Z:"))[5:]
        except:
            sys.stderr.write(repr(self.other) + "\n")
            sys.stderr.write(self.read + "\n")
            raise

    @property
    def ga_ct(self):
        return [x for x in self.other if x.startswith("YC:Z:")]

    def longest_match(self, patt=re.compile("\d+M")):
        return max(int(x[:-1]) for x in patt.findall(self.cigar))


def rname(fq1, fq2=""):
    fq1, fq2 = fq1.split(",")[0], fq2.split(",")[0]
    def name(f):
        n = op.basename(op.splitext(f)[0])
        if n.endswith('.fastq'): n = n[:-6]
        if n.endswith(('.fq', '.r1', '.r2')): n = n[:-3]
        return n
    return "".join(a for a, b in zip(name(fq1), name(fq2)) if a == b) or 'bm'


def bwa_mem(fa, mfq, extra_args, prefix='bwa-meth', threads=1, rg=None,
            calmd=False, paired=True, set_as_failed=None, mismatch_ratio=1):
    conv_fa = convert_fasta(fa, just_name=True)
    if not is_newer_b(conv_fa, (conv_fa + '.amb', conv_fa + '.sa')):
        raise BWAMethException("first run bwameth.py index %s" % fa)

    if not rg is None and not rg.startswith('@RG'):
        rg = '@RG\tID:{rg}\tSM:{rg}'.format(rg=rg)

    # penalize clipping and unpaired. lower penalty on mismatches (-B)
    cmd = "|bwa mem -T 40 -B 2 -L 10 -CM "

    if paired:
        cmd += ("-U 100 -p ")
    cmd += "-R '{rg}' -t {threads} {extra_args} {conv_fa} {mfq}"
    cmd = cmd.format(**locals())
    sys.stderr.write("running: %s\n" % cmd.lstrip("|"))
    as_bam(cmd, fa, prefix, calmd, set_as_failed, mismatch_ratio)


def as_bam(pfile, fa, prefix, calmd=False, set_as_failed=None, mismatch_ratio=1):
    """
    pfile: either a file or a |process to generate sam output
    fa: the reference fasta
    prefix: the output prefix or directory
    set_as_failed: None, 'f', or 'r'. If 'f'. Reads mapping to that strand
                      are given the sam flag of a failed QC alignment (0x200).
    """
    view = "samtools view -bS - | samtools sort -m 2415919104 - -T {bam} -o "
    if calmd:
        cmds = [
            view + "{bam}.tmp",
            "samtools calmd -AbEr {bam}.tmp.bam {fa} > {bam}.bam 2>/dev/null",
            "rm {bam}.tmp.bam"]
    else:
        cmds = [view + "{bam}.bam"]

    cmds.append("samtools index {bam}.bam")
    cmds = [c.format(bam=prefix, fa=fa) for c in cmds]

    sys.stderr.write("writing to:\n%s\n" % cmds[0])

    p = nopen("|" + cmds[0], 'w')
    out = p.stdin
    #out = sys.stdout # useful for debugging
    bam_iter = reader("%s" % (pfile,), header=False, quotechar=None)
    for toks in bam_iter:
        if not toks[0].startswith("@"): break
        handle_header(toks, out)
    else:
        sys.stderr.flush()
        raise Exception("bad or empty fastqs")
    bam_iter2 = chain([toks], bam_iter)
    for read_name, pair_list in groupby(bam_iter2, itemgetter(0)):
        pair_list = [Bam(toks) for toks in pair_list]

        for aln in handle_reads(pair_list, set_as_failed, mismatch_ratio):
            out.write(str(aln) + '\n')

    p.stdin.flush()
    p.stdout.flush()
    p.stdin.close()
    assert p.wait() == 0
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
        #out.write("\t".join(toks) + "\n")
        toks = ["@PG\tID:bwa-meth\tPN:bwa-meth\tVN:%s\tCL:\"%s\"" % (
                         __version__,
                         " ".join(x.replace("\t", "\\t") for x in sys.argv))]
    out.write("\t".join(toks) + "\n")


def handle_reads(alns, set_as_failed, mismatch_ratio=1):

    for aln in alns:
        orig_seq = aln.original_seq
        assert len(aln.seq) == len(aln.qual), aln.read
        # don't need this any more.
        aln.other = [x for x in aln.other if not x.startswith('YS:Z')]
        
        # first letter of chrom is 'f' or 'r'
        direction = aln.chrom[0]
        aln.chrom = aln.chrom.lstrip('fr')

        if not aln.is_mapped():
            aln.seq = orig_seq
            continue

        assert direction in 'fr', (direction, aln)
        aln.other.append('YD:Z:' + direction)

        if mismatch_ratio < 1:
            mismatches = int(next((s for s in aln.other if "NM:i:" in s), None).split(":")[2])
            if mismatches > (aln.longest_match() * mismatch_ratio):
                aln.flag |= 0x200
                aln.other.append('MC:Z:' + aln.chrom)
                aln.other.append('MP:Z:' + str(aln.pos))
                aln.chrom = "*"
                aln.pos = 0


        if set_as_failed == direction:
            aln.flag |= 0x200

        # here we have a heuristic that if the longest match is not 44% of the
        # sequence length, we mark it as failed QC and un-pair it. At the end
        # of the loop we set all members of this pair to be unmapped
        if aln.longest_match() < (len(orig_seq) * 0.44):
            aln.flag |= 0x200  # fail qc
            aln.flag &= (~0x2) # un-pair
            aln.mapq = min(int(aln.mapq), 1)

        mate_direction = aln.chrom_mate[0]
        if mate_direction not in "*=":
            aln.chrom_mate = aln.chrom_mate[1:]

        # adjust the original seq to the cigar
        l, r = aln.left_shift(), aln.right_shift()
        if aln.is_plus_read():
            aln.seq = orig_seq[l:r]
        else:
            aln.seq = comp(orig_seq[::-1][l:r])

    if any(aln.flag & 0x200 for aln in alns):
        for aln in alns:
            aln.flag |= 0x200
            aln.flag &= (~0x2)
    return alns

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
    p.add_argument("-t", "--threads", type=int, default=3)
    p.add_argument("--dbsnp", help="optional dbsnp for GATK calibration")
    p.add_argument("--prefix", help="output prefix", default='bmeth-tab')
    p.add_argument("--trim", help="left, right trim to avoid bias",
                        default="2,2")
    p.add_argument("--map-q", type=int, default=25, help="only tabulate "
                   "methylation for reads with at least this mapping quality")
    p.add_argument("--bissnp", help="path to bissnp jar", required=True)
    p.add_argument("--region", help="region to call variants, e.g. 'chr1'"
            "can be used to aid in parallelization")
    p.add_argument("--format", help="format for output summary."
            " default matches bismark begraph output with: %(default)r"
            " Where cs and ts are the counts of C's and T's and pct is the"
            " percent methylation. Available variables in addition to those"
            " are 'end' and 'ctx', where ctx is the  cpg context (CG/CHH/CHG).",
            default="{chrom}\t{start}\t{start}\t{pct}\t{cs}\t{ts}")
    group = p.add_mutually_exclusive_group()
    group.add_argument("--context", choices=("all", "CG", "CG-strict"),
            default="CG-strict", help="which methylation context to output to the"
            "summary (BED) file. Default 'CG-strict' follows bissnp and only"
            " pulls those sites that are CG in all samples (no CC samples,"
            " for example) 'CG' will pull anything that is CG in any sample:"
            " (MG,SG, YG, CR, CK). 'all' will pull all contexts. Automatically"
            " set to 'all' when '--nome' is specified.")
    group.add_argument("--nome", action='store_true', help="Whether the assay performed"
            " is NOMe and therefore both GpC and CpG sites should be output."
            " Implies '--context all'")
    p.add_argument("bams", nargs="+")

    a = p.parse_args(args)
    assert os.path.exists(a.reference)
    if not os.path.exists(a.reference + ".fai"):
        sys.stderr.write("ERROR: run 'samtools faidx %s' before tabulation\n"
                         % a.reference)
        sys.exit(1)
    trim = list(map(int, a.trim.split(",")))
    if not a.prefix.endswith(('/', '.')):
        a.prefix += "."
    if a.region:
        if op.exists(a.region):
            name = op.basename(a.region)
            if name.endswith('.gz'): name = name[:-3]
            if name.endswith('.bed'): name = name[:-4]
        else:
            name = a.region.replace(":", "-")
        a.prefix += name + "."
    if a.nome:
        a.context = 'all'

    cmd = """\
    java -Xmx24g -jar {bissnp}
        -R {reference}
        -I {bams}
        -T BisulfiteGenotyper
        --trim_5_end_bp {trim5}
        --trim_3_end_bp {trim3}
        -vfn1 {prefix}meth.vcf -vfn2 {prefix}snp.vcf
        -mbq 12
        -minConv 1
        -toCoverage 1000
        -mmq {mapq} {dbsnp} {region} {nome}
        -nt {threads}""".format(
            threads=a.threads,
            dbsnp=("--dbsnp " + a.dbsnp) if a.dbsnp else "",
            region=("-L " + a.region) if a.region else "",
            bissnp=a.bissnp,
            trim5=trim[0],
            trim3=trim[1],
            prefix=a.prefix,
            reference=a.reference,
            mapq=a.map_q,
            bams=" -I ".join(a.bams),
            nome=("-sm GM -out_modes EMIT_VARIANT_AND_CYTOSINES -C WCG,2 -C GCH,2") if a.nome else ""
            ).replace("\n", " \\\n")
    sys.stderr.write(cmd + '\n')
    if not "--format" in args and a.context == "all":
        a.format = a.format.rstrip('\n') + "\t{ctx}"

    contexts = {'all': None,
                'CG-strict': ('CG',),
                'CG': "CG YG SG MG CR CS CK".split()}[a.context]

    run(cmd)

    fmt = a.format.rstrip('\n') + '\n'
    sys.stderr.write(a.prefix + "meth.vcf\n")
    for i, d in enumerate(reader(a.prefix + "meth.vcf",
                       skip_while=lambda toks: toks[0] != "#CHROM",
                       header="ordered")):
        if i == 0:
            samples = list(d.keys())[9:]
            fhs = {}
            for sample in samples:
                fhs[sample] = open("{prefix}{sample}.meth.bed"\
                        .format(prefix=a.prefix, sample=sample), "w")

                fhs[sample].write("#" + fmt.replace("}", "").replace("{", ""))
        if not d['FILTER'] in (".", "PASS"): continue

        d['start'], d['end'] = str(int(d['POS']) - 1), d['POS']
        d['chrom'] = d['CHROM']
        for sample in samples:
            if d[sample] == "./.": continue
            sinfo = dict(zip(d['FORMAT'].split(":"), d[sample].split(":")))
            try:
                d['cs'] = int(sinfo['CM'])  # (M)ethylated
                d['ts'] = int(sinfo['CU'])  # (U)n
            except ValueError:
                continue
            except:
                sys.stderr.write("\nline:%i\t%s\t%s\n" % (i, d, sinfo))
                sys.stderr.write("%s\n" % d[sample])
                raise
            d['ctx'] = sinfo['CP']
            if contexts is not None:
                if not d['ctx'] in contexts: continue

            if d['cs'] + d['ts'] == 0:
                continue
            else:
                d['pct'] = "%.1f" % (100 * d['cs'] / float(d['cs'] + d['ts']))
            fhs[sample].write(fmt.format(**d))

    with open(a.prefix + "command.sh", 'w') as fh:
        fh.write(cmd)

def convert_fqs(fqs):
    script = __file__
    return "'<%s %s c2t %s %s'" % (sys.executable, script, fqs[0],
               fqs[1] if len(fqs) > 1
                      else ','.join(['NA'] * len(fqs[0].split(","))))

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
    p.add_argument("-m", "--mismatch-ratio", type=float, default=1, help="Maximum ratio of mismatches to "
                                                                            "alignment length. Hits with more mismatches"
                                                                            "will be reported as qc-failed (0x200)")
    p.add_argument("-p", "--prefix", default="bwa-meth")
    p.add_argument("--calmd", default=False, action="store_true")
    p.add_argument("--read-group", help="read-group to add to bam in same"
            " format as to bwa: '@RG\\tID:foo\\tSM:bar'")
    p.add_argument('--set-as-failed', help="flag alignments to this strand"
            " as not passing QC (0x200). Targetted BS-Seq libraries are often"
            " to a single strand, so we can flag them as QC failures. Note"
            " f == OT, r == OB. Likely, this will be 'f' as we will expect"
            " reads to align to the original-bottom (OB) strand and will flag"
            " as failed those aligning to the forward, or original top (OT).",
        default=None, choices=('f', 'r'))

    p.add_argument("fastqs", nargs="+", help="bs-seq fastqs to align. Run"
            "multiple sets separated by commas, e.g. ... a_R1.fastq,b_R1.fastq"
            " a_R2.fastq,b_R2.fastq note that the order must be maintained.")

    args, pass_through_args = p.parse_known_args(args)

    # for the 2nd file. use G => A and bwa's support for streaming.
    conv_fqs_cmd = convert_fqs(args.fastqs)

    bwa_mem(args.reference, conv_fqs_cmd, ' '.join(map(str, pass_through_args)), prefix=args.prefix,
             threads=args.threads, rg=args.read_group or
             rname(*args.fastqs), calmd=args.calmd,
             paired=len(args.fastqs) == 2,
             set_as_failed=args.set_as_failed,
             mismatch_ratio=args.mismatch_ratio)


def test():
    ref = "example/ref.fa"
    bwa_index(ref)
    with tempfile.NamedTemporaryFile(delete=True) as rfh:
        prefix=rfh.name
        bwa_mem(ref, convert_fqs(("example/t_R1.fastq.gz",
            "example/t_R2.fastq.gz")), "",
            prefix=rfh.name, threads=2, rg="ex", paired=True)
        assert op.exists("%s.bam" % prefix)
        assert op.exists("%s.bam.bai" % prefix)

if __name__ == "__main__":
    main(sys.argv[1:])
