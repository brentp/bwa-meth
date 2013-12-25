import itertools as it
import sys


fh1 = open(sys.argv[1])
fh2 = open(sys.argv[2])

out1 = open('sim_R1.fastq', 'w')
out2 = open('sim_R2.fastq', 'w')

it1 = it.izip(*[fh1] * 4)
it2 = it.izip(*[fh2] * 4)

for r1, r2 in it.izip(it1, it2):

    rn1 = r1[0].split(":")
    rn2 = r2[0].split(":")
    n = "@%s__%s" % (":".join(rn1[3:-1]).strip(), ":".join(rn2[3:-1]).strip())

    r1, r2 = list(r1), list(r2)
    r1[0] = r2[0] = n
    #r1[0] += "/1\n"
    #r2[0] += "/2\n"
    r1[0] += "_R1\n"
    r2[0] += "_R2\n"

    out1.write("".join(r1))
    out2.write("".join(r2))


