import sys
from toolshed import nopen

read_id = sys.argv[1]
r1 = sys.argv[2]
r2 = sys.argv[3]


o1 = open('o_R1.fq', 'w')
o2 = open('o_R2.fq', 'w')

def print_it(f, of, read_id, a, b):
    found = False
    for i, line in enumerate(nopen(f)):
        if line.startswith("@" + read_id):
            of.write(line)
            found = i
        elif found and i == found + 1:
            of.write(line.replace(a, b))
            #of.write(line)
        elif found:
            of.write(line)
        if found and i > found + 2: break

print_it(r1, o1, read_id, "C", "T")
print_it(r2, o2, read_id, "G", "A")
