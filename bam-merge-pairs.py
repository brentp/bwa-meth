"""
In order to avoid double counting of bases in overlapping paired end reads,
this script accept mapped reads in SAM format, use the position to determine
overlap, and set the lower base-quality of pair overlap to 0.
It is assumed that the reads are name-sorted as when they come from the
aligner.

Usage:

    bwa mem reference.fa R1.fq R2.fq \
            | bam-merge-pairs.py \
            | samtools view -bS - > merged.bam
"""
from __future__ import print_function
from __future__ import division
import itertools as it
from operator import itemgetter
from bwameth import Bam
import sys

line_iter = (x.rstrip("\r\n").split("\t") for x in sys.stdin)


for g, pair in it.groupby(line_iter, itemgetter(0)):
    if g.startswith("@"):
        print("\n".join("\t".join(line) for line in pair))
        continue
    pair = list(pair)
    if len(pair) == 1:
        print "\t".join(pair[0])
        continue

    assert len(pair) == 2, pair
    left, right = [Bam(b) for b in pair]

    # TODO: use left_shift(), right_shift()
    if left.pos + left.tlen < right.pos:
        print(str(left))
        print(str(right))
        continue

    # there is overlap 
    # L ------------------>
    # R              <------------------
    ovl_bases = right.pos, left.pos + left.tlen




