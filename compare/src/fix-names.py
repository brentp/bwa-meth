import itertools as it
import sys
import gzip


fh1 = open(sys.argv[1])
fh2 = open(sys.argv[2])

out1 = gzip.open('../dnemsim_R1.fastq.gz', 'w')
out2 = gzip.open('../dnemsim_R2.fastq.gz', 'w')

it1 = it.izip(*[fh1] * 4)
it2 = it.izip(*[fh2] * 4)

for i, (r1, r2) in enumerate(it.izip(it1, it2), 1):

    rn1 = r1[0].split(":")
    rn2 = r2[0].split(":")
    #n = "@%s-PAIR-%s\n" % (":".join(rn1[3:-1]).strip(), ":".join(rn2[3:-1]).strip())

    #@:0/1:100:chr12a:10030901:+ original
    chrom, pos1 = rn1[3].rstrip('ab'), int(rn1[4])
    chrom, pos2 = rn2[3].rstrip('ab'), int(rn2[4])

    pos1, pos2 = sorted((pos1, pos2))


    n = "@%i_%s:%i-%i\n" % (i, chrom, pos1, pos2)

    r1, r2 = list(r1), list(r2)
    r1[0] = r2[0] = n

    out1.write("".join(r1))
    out2.write("".join(r2))


