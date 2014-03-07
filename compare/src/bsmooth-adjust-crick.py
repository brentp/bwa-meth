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
    flag = int(toks[1])
    if chrom != "*" and not flag & 0x4:
        pos = int(toks[3])
        pos = max(1, seq_lens[chrom] - pos - len(toks[9]))
        toks[3] = str(pos)
    print "\t".join(toks)
