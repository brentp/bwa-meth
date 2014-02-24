. ./common.sh

PATH=$PATH:/opt/mpich2/gnu/bin/

GEN=$(dirname $REF)/
cmd="bison --directional -N 1 --very-sensitive -p 12 -g $GEN -o"
prog=bison

rm -rf results/$prog/
mkdir -p results/$prog/
rm -f logs/$prog-$name.err logs/$prog-$name.out

echo "mpiexec $cmd results/$prog/ -1 $FQ1 -2 $FQ2" \
        | bsub -J $prog-$name \
               -R "span[ptile=5]" \
               -a openmpi \
               -q shared \
               -e logs/$prog-$name.err \
               -o logs/$prog-$name.out -n 5
exit

rm -rf results/trim/$prog/
mkdir -p results/trim/$prog/
rm -f logs/trim-$prog-$name.err logs/trim-$prog-$name.out

echo "mpiexec $cmd results/trim/$prog/ -1 $TRIM_FQ1 -2 $TRIM_FQ2" \
        | bsub -J trim-$prog-$name \
               -R "span[ptile=12]" \
               -a mvapich \
               -q shared \
               -e logs/trim-$prog-$name.err \
               -o logs/trim-$prog-$name.out -n 60
