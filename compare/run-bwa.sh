. ./common.sh

echo "python ../bwa-meth.py --reference $REF \
        -t 22 -p results/bwa-$name $FQ1 $FQ2" \
        | bsub -J bwa-$name \
               -e logs/bwa-$name.err \
               -o logs/bwa-$name.out -n 12


echo "python ../bwa-meth.py --reference $REF \
        -t 22 -p results/trim/bwa-$name $TRIM_FQ1 $TRIM_FQ2" \
        | bsub -J trim-bwa-$name \
               -e logs/trim-bwa-$name.err \
               -o logs/trim-bwa-$name.out -n 12
