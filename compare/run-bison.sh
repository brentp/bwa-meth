. ./common.sh

PATH=$PATH:/opt/mpich2/gnu/bin/

GEN=$(dirname $REF)/
cmd="bison --directional -N 1 --very-sensitive-local -p 12 -g $GEN -o"
prog=bison

rm -rf results/$prog/
mkdir -p results/$prog/
rm -f logs/$prog-$name.err logs/$prog-$name.out

mv="mv results/$prog/${name}_R1.bam results/$prog-$name.bam"
echo "mpiexec $cmd results/$prog/ -1 $FQ1 -2 $FQ2; $mv" \
        | bsub -J $prog-$name \
               -n 3 \
               -R "span[ptile=5]" \
               -a mvapich \
               -q shared \
               -e logs/$prog-$name.err \
               -o logs/$prog-$name.out 

rm -rf results/trim/$prog/
mkdir -p results/trim/$prog/
rm -f logs/trim-$prog-$name.err logs/trim-$prog-$name.out

mv="mv results/trim/$prog/${name}_R1.bam results/trim/$prog-$name.bam"
echo "mpiexec $cmd results/trim/$prog/ -1 $TRIM_FQ1 -2 $TRIM_FQ2; $mv" \
        | bsub -J trim-$prog-$name \
               -n 3 \
               -R "span[ptile=5]" \
               -a mvapich \
               -q shared \
               -e logs/trim-$prog-$name.err \
               -o logs/trim-$prog-$name.out
