import os.path as op
import re
import sys
from toolshed import reader, nopen
from collections import defaultdict

import numpy as np
from itertools import cycle
import pylab as pl
import seaborn

colors = cycle(seaborn.color_palette('Set1', 8))

BASES = False

def name(bam):
    return op.basename(bam).rsplit(".", 1)[0].split("-")[0]
    return op.basename(bam).rsplit(".", 1)[0].split("-")[0] + ("-trim" if "trim" in bam else "")
    return op.dirname(bam)

def count_bases(cigar, patt=re.compile("\d+[IM]")):
    return sum(int(p[:-1]) for p in patt.findall(cigar))

def count_on_off(bam, flags, pad, bases=BASES):

    #if not "last" in bam: flags = flags + " -f 0x2"

    on_count = [0] * 256
    off_count = [0] * 256
    rcounts = defaultdict(int)

    posn = re.compile(".*?_(.+):(\d+)-(\d+)")

    # chr1b:3001315:+__chr1b:3001467:-        99      chr1    3001316 60 100M 
    print >>sys.stderr, bam
    for toks in reader("|samtools view {flags} {bam}".format(**locals()),
            header=False):
        if toks[0][0] == "@": continue
        if toks[0].endswith(("_R1", "_R2")):
            toks[0] = toks[0][:-3]

        cigar = toks[5]

        qual = int(toks[4])
        if toks[0].endswith("_BAD"):
            off_count[qual] += (count_bases(cigar) if bases else 1)
            continue

        rname, flag = toks[0], int(toks[1])
        rcounts[rname] += 1
        # @2_chr4:72040172-72040686_R1
        if rcounts[rname] > 2:
            raise Exception("%s:%s", (bam, rname))

        try:
            chrom, start, end = re.match(posn, rname).groups(0)
        except:
            print >>sys.stderr, rname
            raise

        pos = int(start if flag & 0x40 else end)
        on = chrom == toks[2] and abs(pos - int(toks[3])) <= pad
        if on:
            on_count[qual] += (count_bases(cigar) if bases else 1)
        else:
            off_count[qual] += (count_bases(cigar) if bases else 1)

    on_count = np.cumsum(on_count[::-1])[::-1]
    off_count = np.cumsum(off_count[::-1])[::-1]
    return off_count, on_count

FLAGS="-F%i" % (0x4 | 0x100 | 0x200)
def main(bams, reads=None, flags=FLAGS, pad=12002):
    if not reads.isdigit():
        reads = 2 * float(nopen("|bioawk -c fastx 'END { print NR }' %s" % reads).next())
    else:
        reads = 2.0 * int(reads)

    if any('trim' in b for b in bams):
        assert all('trim' in b for b in bams), [b for b in bams if not 'trim'
                in b]
    counts = {}
    names, pts = [], []

    denom = float(reads)
    if BASES: denom *= 100.
    for bam in bams:
        counts[bam] = count_on_off(bam, flags, pad)


    out = sys.stdout
    print >>out, "qual\tmethod\toff\ton"

    for qual in range(0, 256):
        for b in bams:
            print >>out, "{qual}\t{bam}\t{off}\t{on}".format(
            qual=qual, bam=name(b),
            off=counts[b][0][qual] / denom,
            on=counts[b][1][qual] / denom)
    print >>sys.stderr, "wrote", out.name

    for bam in bams:

        # multiply by 100 bases per read.
        color = next(colors)

        symbol = 'o' if len(set(counts[bam][0])) < 3 else '.'
        nmax = max(i for i, v in enumerate(counts[bam][1]) if v > 0)
        p, = pl.plot(counts[bam][0][1:nmax + 1] / denom,
                counts[bam][1][1:nmax + 1] / denom,
                color=color,
                mec=color,
                mfc=color,
                marker=symbol, label=name(bam))
        pts.append(p)
        names.append(name(bam))
    pl.xlabel('off target')
    pl.ylabel('on target')
    pl.legend(pts, names, loc='lower right')
    pl.xlim(xmin=0)
    pl.ylim(ymin=0)
    print >>sys.stderr, reads

    pl.show()

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--reads", help="reads file", required=True)
    p.add_argument("bams", nargs="+")

    a = p.parse_args()

    main(a.bams, reads=a.reads)
