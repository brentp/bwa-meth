import os
import os.path as op
import sys
import tempfile
from toolshed import reader

import numpy as np
from itertools import cycle
import pylab as pl

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

ON  = "| samtools view {bam} -L {regions} {flags}"
OFF = "| bedtools intersect -ubam -abam {bam} \
            -b <(grep -v '@' {regions} | sort -k1,1 -k2,2n) -wa -v \
       | samtools view - {flags}"

def main(regions, bams, reads=None, flags="-F4", pad=50):
    r2 = open(tempfile.mktemp(), 'w')
    for toks in reader(regions, header=False):
        if toks[0][0] == "@" or not (toks[1] + toks[2]).isdigit(): continue
        print >>r2, "\t".join(toks)
    r2.flush()
    regions = r2.name

    reads = nopen("|bioawk -c fastx 'END { print NR }' %s" % reads)

    counts = {}
    colors = cycle('rgbk')
    for bam in bams:
        off, on = OFF.format(**locals()), ON.format(**locals())
        counts[bam] = list(pmap(counter, [off, on]))
        if not isinstance(reads, int):
            # should divide by 4, but we also have pairs so...
            reads = int(reads.next()) * 2


        symbol = 'o' if len(set(counts[bam][0])) < 3 else '.'
        pl.plot(counts[bam][0] / float(reads), counts[bam][1] / float(reads),
                '%s%s' % (colors.next(), symbol), label=op.basename(bam))

    pl.xlabel('off target')
    pl.ylabel('on target')
    pl.legend(loc='lower right')
    pl.xlim(xmin=0)
    pl.ylim(ymin=0)
    print reads
    pl.show()
    os.unlink(r2.name)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--reads", help="reads file")
    p.add_argument("regions")
    p.add_argument("bams", nargs="+")

    a = p.parse_args()

    main(a.regions, a.bams, reads=a.reads)
