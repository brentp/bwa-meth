import sys

seq_lens = {}

for line in sys.stdin:
    if line.startswith('@'):
        print line,
        if line.startswith('@SQ'):
            #@SQ    SN:chr1 LN:195471971
            _, sn, ln = line.rstrip().split()
            sn, ln = sn.split(":")[1], int(ln.split(":")[1])
            seq_lens[sn] = ln
        continue
    toks = line.rstrip().split("\t")
    chrom = toks[2]
    pos = int(toks[3])
    pos = seq_lens[chrom] - pos + 1
    toks[3] = str(pos)
    print "\t".join(toks)
