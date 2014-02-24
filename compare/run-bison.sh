. ./common.sh


cmd="bison -N 1 --very-sensitive -p 12 -o"
name=bison

rm -f logs/$name.err logs/$name.out
echo "$cmd results/$name -1 $FQ1 -2 $FQ2" \
        | bsub -J $name \
               -e logs/$name.err \
               -o logs/$name.out -n 12


rm -f logs/trim-$name.err logs/trim-$name.out

echo "$cmd results/trim/$name -1 $TRIM_FQ1 -2 $TRIM_FQ2" \
        | bsub -J trim-$name \
               -e logs/trim-$name.err \
               -o logs/trim-$name.out -n 12
