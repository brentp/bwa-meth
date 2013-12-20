import os
import os.path as op
import sys
from toolshed import reader, nopen
from collections import defaultdict

import numpy as np
from itertools import cycle
import pylab as pl

def name(bam):
    return op.basename(bam).rsplit(".", 1)[0].split("-")[0]

def count_on_off(bam, flags, pad):

    on_count = [0] * 256
    off_count = [0] * 256
    rcounts = defaultdict(int)

    # chr1b:3001315:+__chr1b:3001467:-        99      chr1    3001316 60 100M 
    print bam
    for toks in reader("|samtools view {flags} {bam}".format(**locals()),
            header=False):
        if toks[0][0] == "@": continue
        # chr1b:3001315:+__chr1b:3001467:-
        rname, flag = toks[0], int(toks[1])
        rcounts[rname] += 1
        if rcounts[rname] > 2:
            raise Exception("%s:%s", (bam, rname))
        try:
            lname, rname = rname.split("__")
        except:
            print >>sys.stderr, toks
            print >>sys.stderr, rname
            raise
        name = lname if flag & 0x40 else rname if flag & 0x80 else None
        assert name is not None
        chrom, pos = name.split(":")[:2]
        chrom, pos = chrom[:-1], int(pos)
        on = chrom == toks[2] and abs(pos - int(toks[3])) <= pad
        qual = int(toks[4])
        if on:
            on_count[qual] += 1
        else:
            off_count[qual] += 1

    on_count = np.cumsum(on_count[::-1])[::-1]
    off_count = np.cumsum(off_count[::-1])[::-1]
    return off_count, on_count



def main(bams, reads=None, flags="-f2 -F0x100", pad=101):
    reads = 2 * float(nopen("|bioawk -c fastx 'END { print NR }' %s" % reads).next())
    counts = {}
    colors = cycle('rgbk')
    for bam in bams:
        counts[bam] = count_on_off(bam, flags, pad)

        symbol = 'o' if len(set(counts[bam][0])) < 3 else '.'
        pl.plot(counts[bam][0] / float(reads), counts[bam][1] / float(reads),
                '%s%s' % (colors.next(), symbol), label=name(bam))

    pl.xlabel('off target')
    pl.ylabel('on target')
    pl.legend(loc='lower right')
    pl.xlim(xmin=0)
    pl.ylim(ymin=0)
    print reads
    pl.show()

    out = open('sim-qual-summary.txt', 'w')
    print >>out, "qual\tmethod\toff\ton"

    for qual in range(0, 256):
        for b in bams:
            print >>out, "{qual}\t{bam}\t{off}\t{on}".format(
            qual=qual, bam=name(b),
            off=counts[b][0][qual] / reads,
            on=counts[b][1][qual] / reads)
    print >>sys.stderr, "wrote", out.name

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--reads", help="reads file")
    p.add_argument("bams", nargs="+")

    a = p.parse_args()

    main(a.bams, reads=a.reads)
