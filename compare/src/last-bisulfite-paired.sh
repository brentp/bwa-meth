#! /bin/bash

# Align paired bisulfite-converted DNA reads to a genome.

# This assumes that reads1.fastq are all from the converted strand
# (i.e. they have C->T conversions) and reads2.fastq are all from the
# reverse-complement (i.e. they have G->A conversions).

# "GNU parallel" needs to be installed.

[ $# -eq 5 ] || {
    cat <<EOF
Typical usage:

  lastdb -w 2 -u bisulfite_f.seed my_f mygenome.fa
  lastdb -w 2 -u bisulfite_r.seed my_r mygenome.fa

  $(basename $0) my_f my_r reads1.fastq reads2.fastq readgroup_id > results.maf

EOF
    exit 2
}

# Try to get the LAST programs into the PATH, if they aren't already:
PATH=$PATH:$(dirname $0)/../src:$(dirname $0)/../scripts

tmp=/scratch/brentp/$$
trap 'rm -f $tmp.*' EXIT

cat > $tmp.fmat << 'EOF'
    A   C   G   T
A   6 -18 -18 -18
C -18   6 -18   3
G -18 -18   6 -18
T -18 -18 -18   3
EOF

cat > $tmp.rmat << 'EOF'
    A   C   G   T
A   3 -18 -18 -18
C -18   6 -18 -18
G   3 -18   6 -18
T -18 -18 -18   6
EOF

cat > $tmp.script << 'EOF'
t=$1.$$

lastal -m10 -p $1.fmat -s1 -Q1 -e120 -i1 "$2" "$4" > $t.t1f
lastal -m10 -p $1.rmat -s0 -Q1 -e120 -i1 "$3" "$4" > $t.t1r
last-merge-batches.py $t.t1f $t.t1r > $t.t1
rm $t.t1f $t.t1r

lastal -m10 -p $1.fmat -s0 -Q1 -e120 -i1 "$2" "$5" > $t.t2f
lastal -m10 -p $1.rmat -s1 -Q1 -e120 -i1 "$3" "$5" > $t.t2r
last-merge-batches.py $t.t2f $t.t2r > $t.t2
rm $t.t2f $t.t2r

last-pair-probs.py -m0.1 $t.t1 $t.t2 |
perl -F'(\s+)' -ane '$F[12] =~ y/ta/CG/ if /^s/ and $s++ % 2; print @F'
rm $t.t1 $t.t2
EOF

# use less as it does zless too
# Convert C to t, and all other letters to uppercase:
perl -pe 'y/Cca-z/ttA-Z/ if $. % 4 == 2' <(less $3) | split -l1800000 -a5 - $tmp.1

# Convert G to a, and all other letters to uppercase:
perl -pe 'y/Gga-z/aaA-Z/ if $. % 4 == 2' <(less $4) | split -l1800000 -a5 - $tmp.2

parallel --noswap -j 4 --gnu --xapply sh $tmp.script $tmp "$1" "$2" ::: $tmp.1* ::: $tmp.2* > $tmp.merge.maf

rg=$5
maf-convert.py sam -d $tmp.merge.maf -r "ID:$rg SM:$rg"
