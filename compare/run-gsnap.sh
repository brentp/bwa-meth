. ./common.sh


cmd="set -e; python src/gsnap-meth.py --stranded -r $REF -t 17 "


rm -f logs/gsnap-$name.err logs/gsnap-$name.out
echo "$cmd --prefix results/gsnap-$name $FQ1 $FQ2" \
        | bsub -J gsnap-$name \
               -e logs/gsnap-$name.err \
               -o logs/gsnap-$name.out -n 12


rm -f logs/trim-gsnap-$name.err logs/trim-gsnap-$name.out
echo "$cmd --prefix results/trim/gsnap-$name $TRIM_FQ1 $TRIM_FQ2" \
        | bsub -J trim-gsnap-$name \
               -e logs/trim-gsnap-$name.err \
               -o logs/trim-gsnap-$name.out -n 12
