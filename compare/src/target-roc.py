import os
import os.path as op
import sys
import tempfile
from toolshed import reader

import numpy as np
from itertools import cycle
import pylab as pl

def name(bam):
    return op.basename(bam).rsplit(".", 1)[0].split("-")[0]

def counter(fname):
    fname = fname[0] if not isinstance(fname, basestring) else fname
    print >>sys.stderr, fname
    qual_count = [0] * 256
    for sam_line in (x.rstrip().split("\t") for x in nopen(fname) if not
            x.startswith('@')):
        qual = int(sam_line[4])
        qual_count[qual] += 1
    # each qual should be the sum of all quals with a lower qual than it
    return np.cumsum(qual_count[::-1])[::-1]

from toolshed import nopen, pmap

def count_both(bam, regions, flags):
    if not "last" in bam: flags += " -f0x2"
    off, on = OFF.format(**locals()), ON.format(**locals())
    return bam, map(counter, [off, on])

def count_bam(bam, flags):
    return int(next(nopen("|samtools view -c {flags} {bam}".format(**locals()))))


ON  = "| samtools view {bam} -L {regions} {flags}"
OFF = "| bedtools intersect -ubam -abam {bam} -b {regions} -wa -v \
       | samtools view - {flags}"

def main(regions, bams, reads=None, flags="-F%i" % (0x100 | 0x4 | 0x200),
        pad=100):
    r2 = open(tempfile.mktemp(), 'w')
    for toks in reader(regions, header=False):
        if toks[0][0] == "@" or not (toks[1] + toks[2]).isdigit(): continue
        toks[1] = str(int(toks[1]) - pad)
        toks[2] = str(int(toks[2]) + pad)
        print >>r2, "\t".join(toks)
    r2.flush()
    regions = r2.name

    if reads != "bam":
        reads = int(nopen("|bioawk -c fastx 'END { print NR }' %s" % reads).next()) * 2.0

    counts = {}
    colors = cycle('rgbkmy')
    bam_reads = {}

    counts = dict(pmap(count_both, ((bam, regions, flags)
                            for bam in bams)))

    for bam in bams:
        nreads = count_bam(bam, flags) if reads == "bam" else reads
        bam_reads[bam] = nreads
        symbol = 'o' if len(set(counts[bam][0])) < 3 else '.'
        pl.plot(counts[bam][0] / float(nreads), counts[bam][1] / float(nreads),
                '%s%s' % (colors.next(), symbol), label=name(bam))

    pl.xlabel('off target')
    pl.ylabel('on target')
    pl.legend(loc='lower right')
    pl.xlim(xmin=0)
    pl.ylim(ymin=0)
    pl.show()
    os.unlink(r2.name)

    out = sys.stdout
    print >>out, "qual\tmethod\toff\ton"

    for qual in range(0, 256):
        for b in bams:
            print >>out, "{qual}\t{bam}\t{off}\t{on}".format(
            qual=qual, bam=name(b),
            off=counts[b][0][qual] / bam_reads[bam],
            on=counts[b][1][qual] / bam_reads[bam])
    print >>sys.stderr, "wrote", out.name

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--reads", help="reads file")
    p.add_argument("regions")
    p.add_argument("bams", nargs="+")

    a = p.parse_args()

    main(a.regions, a.bams, reads=a.reads)
