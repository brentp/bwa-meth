. ./common.sh

PATH=/opt/mpich2/gnu/bin/:$PATH

GEN=$(dirname $REF)/

cmd="bison --quiet --directional --very-sensitive-local -N 1 -p 6 -g $GEN -o"
prog=bison

rm -rf results/$prog-$name/
mkdir -p results/$prog-$name/
rm -f logs/$prog-$name.err logs/$prog-$name.out

echo "mpiexec -np 3 $cmd results/$prog-$name/ -1 $FQ1 -2 $FQ2" \
        | bsub -J $prog$name \
               -R "span[ptile=1]" \
               -a mvapich \
               -q shared \
               -e logs/$prog-$name.err \
               -o logs/$prog-$name.out \
               -n 3

rm -rf results/trim/$prog-$name/
mkdir -p results/trim/$prog-$name/
rm -f logs/trim-$prog-$name.err logs/trim-$prog-$name.out

echo "mpiexec -np 3 $cmd results/trim/$prog-$name/ -1 $TRIM_FQ1 -2 $TRIM_FQ2" \
        | bsub -J trim-$prog-$name \
               -R "span[ptile=1]" \
               -a mvapich \
               -q shared \
               -e logs/trim-$prog-$name.err \
               -o logs/trim-$prog-$name.out \
               -n 3

