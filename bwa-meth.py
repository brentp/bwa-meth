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

from itertools import groupby, izip, count
from toolshed import nopen, reader
from collections import defaultdict
import string

DEBUG = False

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

def convert_reads(fq, ga=False, out=sys.stdout):
    char_a, char_b = 'Ga' if ga else 'Ct'
    print >>sys.stderr, "converting %s to %s in %s" % (char_a, char_b, fq)
    fq = nopen(fq)
    for (name, seq, _, qual) in izip(*[fq] * 4):
        seq = seq.upper()
        # keep original sequence as name.
        name = " ".join((name.split(" ")[0],
                        "YS:Z:" + seq.rstrip('\r\n') +
                        "\tYC:Z:" + char_a + char_b + '\n'))
        out.write("".join((name, seq.replace(char_a, char_b) , "+\n", qual)))

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

    @classmethod
    def cig_len(self, cigs):
        return sum(c[0] for c in cigs if c[1] in ("M", "D", "N", "EQ", "X", "P"))

    def right_end(self):
        # http://www.biostars.org/p/1680/#1682
        length = sum(c[0] for c in self.cigs() if c[1] in ("M", "D", "N", "EQ", "X", "P"))
        return self.start + length - 1

    # This slow...
    def match_string(self, oseq, rseq):
        o, r = [], []
        opos, rpos = 0, 0
        for n, cig in self.cigs():
            if cig == "I":
                r.append('.' * n)
                o.append(oseq[opos:opos + n])
                opos += n
            elif cig in "SH":
                o.append(oseq[opos:opos + n])
                r.append("." * n)
                opos += n

            elif cig == "D":
                o.append("." * n)
                r.append(rseq[rpos:rpos + n])
                rpos += n
            elif cig in ("M",):
                o.append(oseq[opos:opos + n])
                r.append(rseq[rpos:rpos + n])
                opos += n
                rpos += n
            else:
                1/0
        return "".join(r), "".join(o)

    def left_shift(self):
        left = 0
        for n, cig in self.cigs():
            if cig == "M": break
            if cig != "H":
                left += n
        return left

    @property
    def start(self):
        return self.pos

    @property
    def len(self):
        return abs(self.tlen)

    @property
    def original_seq(self):
        return next(x for x in self.other if x.startswith("YS:Z:"))[5:]

    @property
    def ga_ct(self):
        return [x for x in self.other if x.startswith("YC:Z:")]

def rname(fq1, fq2):
    name = lambda f:op.basename(op.splitext(f)[0])
    return "".join(a for a, b in zip(name(fq1), name(fq2)) if a == b)

def bwa_meth(ref_fasta, fastqs, prefix=".", extra_args="", threads=1, mapq=0,
        rg=None):
    assert len(fastqs) in (1, 2)
    return bwa_mem(ref_fasta, fastqs, extra_args, prefix=prefix,
                   threads=threads, mapq=mapq, rg=rg)


def bwa_mem(fa, fqs, extra_args, prefix='bwa-meth', threads=1, mapq=0, rg=None):
    conv_fa = convert_fasta(fa)
    bwa_index(conv_fa)
    if rg is None:
        rg = rname(fqs[0], fqs[1])
        rg = '@RG\tID:{rg}\tSM:{rg}'.format(rg=rg)

    fq_str = " ".join(fqs)

    cmd = ("|bwa mem -C -M -R '{rg}' -t {threads} {extra_args} {conv_fa} {fq_str}"
                    ).format(**locals())
    print >>sys.stderr, "running: %s" % cmd
    cs, ts = tabulate(cmd, fa, prefix, conv_fa=conv_fa, mapq=mapq)


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
    cs = defaultdict(lambda: defaultdict(int))
    ts = defaultdict(lambda: defaultdict(int))

    out = nopen(prefix + ".sam.gz", "w")
    lengths = {}
    for toks in reader("%s" % (pfile, ), header=False):
        if toks[0].startswith("@"):
            if toks[0].startswith("@SQ"):
                sq, sn, ln = toks
                # we have f and r, only print out f
                sn = sn.split(":")[1]
                if sn.startswith('r'): continue
                lengths[sn[1:]] = int(ln.split(":")[1])
                toks[1] = toks[1].replace(":f", ":")
                print >>out, "\t".join(toks)
            continue
        if toks[2] == "*": # chrom
            print >>out, "\t".join(toks)

        aln = Bam(toks)
        orig_seq = aln.original_seq

        # first letter of chrom is 'f' or 'r'
        direction = aln.chrom[0]
        aln.chrom = aln.chrom.lstrip('fr')

        if not aln.is_mapped():
            print >>out, str(aln)
            continue
        assert direction in 'fr', (direction, toks[2], aln)

        if direction == 'r':
            aln.seq = comp(aln.seq)
            orig_seq = comp(orig_seq)

        if aln.chrom_mate[0] != "*":
            aln.chrom_mate = aln.chrom_mate[1:]

        if aln.mapq < mapq:
            print >>out, str(aln)
            continue

        #"""
        left, right = aln.pos, aln.right_end() + 1
        gseq = faseq(fa, aln.chrom, left, right).upper()

        oseq = comp(orig_seq)[::-1] if aln.is_minus_read() else orig_seq

        if debug:
            cseq = faseq(conv_fa, direction + aln.chrom, left, right).upper()

            print >>sys.stderr, "==", direction, aln.cigar, aln.flag, aln.is_minus_read(), aln.mapq, len(aln.seq)
            print >>sys.stderr, "original seq     >>", oseq
            print >>sys.stderr, "mapped   seq     >>", aln.seq
            print >>sys.stderr, "converted genome >>", (aln.left_shift() * " ") + (comp(cseq) if direction == "r" else cseq)
            print >>sys.stderr, "original  genome >>", (aln.left_shift() * " ") + gseq
            print >>sys.stderr, ""
            print >>sys.stderr, "\n".join(aln.match_string(oseq, gseq))
            print >>sys.stderr, ""

        # This section is quite slow. can optimize quite a bit.
        ref, match = aln.match_string(oseq, gseq)
        ref_letter = "G" if direction == "R" else "C"
        for lpos, r, m in (lrm for lrm in izip(ref, match, count(left))
                           if lrm[1] == ref_letter):
            # TODO: check mismatches here since we know more
            if r == m:
                cs[aln.chrom][lpos] += 1
            else: # should check C=>T, G=>A
                ts[aln.chrom][lpos] += 1
        #"""

        print >>out, str(aln)
    #"""
    out.close()
    fmt = "{chrom}\t{pos}\t{pos}\t{pct_meth}\t{c}\t{t}"

    chroms = sorted(set(ts.keys() + cs.keys()))
    try:
        with nopen(prefix + ".summary.bed.gz", "w") as fh:
            for chrom in chroms:
                tsc, csc = ts[chrom], cs[chrom]
                posns = set(tsc.keys())
                posns.update(csc.keys())
                ccs, tts = 0, 0
                for pos in sorted(posns):
                    c, t = tsc.get(pos, 0), csc.get(pos, 0)
                    pct_meth = c / float(c + t)
                    print >>fh, fmt.format(**locals())
                    ccs += c
                    tts += t
                print >>sys.stderr, chrom, ccs / float(ccs + tts), ccs, tts
        print >>sys.stderr, "wrote", fh.name
    except:
        os.unlink(prefix + ".summary.bed.gz")
        raise


    #"""
    try:
        run("zless %s | samtools view -hbS - \
                | samtools sort -m 2G - %s \
                && samtools index %s.bam" \
                % (out.name, out.name[:-4], out.name[:-4]))
    except:
        os.unlink(out.name[:-4] + ".bam")
        if op.exists(out.name[:-4] + ".bam.bai"):
            os.unlink(out.name[:-4] + ".bam.bai")
        raise
    return cs, ts

def faseq(fa, chrom, start, end):
    end = int(end) - 1
    return "".join(x.strip() for i, x in \
            enumerate(nopen("|samtools faidx %s %s:%s-%s" \
            % (fa, chrom, int(start), end))) if i > 0)

def main(args):

    if len(args) > 0 and args[0] in ("c2t", "g2a"):
        # catch these args to convert reads on the fly and stream to bwa
        sys.exit(convert_reads(args[1], ga=args[0] == "g2a"))

    import argparse
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--reference", help="reference fasta")
    p.add_argument("-t", "--threads", type=int, default=6)
    p.add_argument("-p", "--prefix", default="bwa-meth")
    p.add_argument("--read-group", help="read-group to add to bam in same"
            " format as to bwa: '@RG\tID:foo\tSM:bar'")
    p.add_argument("--map-q", type=int, default=10, help="only tabulate "
                   "methylation for reads with at least this mapping quality")
    p.add_argument("fastqs", nargs="+", help="bs-seq fastqs to align")

    args = p.parse_args(args)
    # for the 2nd file. use G => A and bwa's support for streaming.
    conv_fqs = ["'<python %s %s %s'" % (__file__, conv, fq) for conv, fq in
                                zip(('c2t', 'g2a'), args.fastqs)]
    bwa_meth(args.reference, conv_fqs, prefix=args.prefix,
            threads=args.threads, mapq=args.map_q, rg=args.read_group)

if __name__ == "__main__":
    main(sys.argv[1:])
