. ./common.sh

cmd="python ../bwameth.py --reference $REF -t 22 --set-as-failed f -p "
prog=bwar

rm -f logs/${prog}-$name.err logs/${prog}-$name.out
echo "$cmd results/${prog}-$name $FQ1 $FQ2" \
        | bsub -J ${prog}-$name \
               -e logs/${prog}-$name.err \
               -o logs/${prog}-$name.out -n 12


rm -f logs/trim-${prog}-$name.err logs/trim-${prog}-$name.out
echo "$cmd results/trim/${prog}-$name $TRIM_FQ1 $TRIM_FQ2" \
        | bsub -J trim-${prog}-$name \
               -e logs/trim-${prog}-$name.err \
               -o logs/trim-${prog}-$name.out -n 12
