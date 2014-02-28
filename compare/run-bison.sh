. ./common.sh


PATH=/opt/mpich2/gnu/bin/:$PATH

GEN=$(dirname $REF)/


cmd="bison --quiet --directional --local --very-sensitive-local --score-min L,-0.6,-0.6 -N 1 -p 4 -g $GEN -o"
prog=bison

rm -rf results/$prog/
mkdir -p results/$prog/
rm -f logs/$prog-$name.err logs/$prog-$name.out

mv="mv results/$prog/${name}_R1.bam results/$prog-$name.bam"
echo "mpiexec -np 3 $cmd results/$prog/ -1 $FQ1 -2 $FQ2" \
        | bsub -J $prog-$name \
               -R "span[ptile=1]" \
               -a mvapich \
               -q shared \
               -e logs/$prog-$name.err \
               -o logs/$prog-$name.out \
               -n 3

echo $mv | bsub -w "done('$prog-$name')" -e /dev/null -o /dev/null -J mv


rm -rf results/trim/$prog/
mkdir -p results/trim/$prog/
rm -f logs/trim-$prog-$name.err logs/trim-$prog-$name.out

mv="mv results/trim/$prog/${name}_R1.bam results/trim/$prog-$name.bam"
echo "mpiexec -np 3 $cmd results/trim/$prog/ -1 $TRIM_FQ1 -2 $TRIM_FQ2" \
        | bsub -J trim-$prog-$name \
               -R "span[ptile=1]" \
               -a mvapich \
               -q shared \
               -e logs/trim-$prog-$name.err \
               -o logs/trim-$prog-$name.out \
               -n 3

echo $mv | bsub -w "done('trim-$prog-$name')" -e /dev/null -o /dev/null -J mv-trim
