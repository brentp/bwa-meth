#lastdb -w 2 -u ~/src/last-hg/examples/bisulfite_f.seed $REF.last_f $REF 
#lastdb -w 2 -u ~/src/last-hg/examples/bisulfite_r.seed $REF.last_r $REF

. ./common.sh


rm -f logs/last-$name.err logs/last-$name.out logs/trim-last-$name.err logs/trim-last-$name.out

echo "
./src/last-bisulfite-paired.sh $REF.last_f $REF.last_r $FQ1 $FQ2 $name \
    | samtools view -bS - \
    | samtools sort - $OUT/last-$name
samtools fixmate $OUT/last-$name.bam $OUT/last-$name.fix.bam
samtools sort $OUT/last-$name.fix.bam $OUT/last-$name
rm $OUT/last-$name.fix.bam
samtools index $OUT/last-$name.bam
" | bsub -J last-$name \
         -R "span[hosts=1]" \
         -e logs/last-$name.err \
         -o logs/last-$name.out -n 8

echo "
./src/last-bisulfite-paired.sh $REF.last_f $REF.last_r $TRIM_FQ1 $TRIM_FQ2 $name \
    | samtools view -bS - \
    | samtools sort - $OUT/trim/last-$name
samtools fixmate $OUT/trim/last-$name.bam $OUT/trim/last-$name.fix.bam
samtools sort $OUT/trim/last-$name.fix.bam $OUT/trim/last-$name
rm $OUT/trim/last-$name.fix.bam
samtools index $OUT/trim/last-$name.bam
" | bsub -J trim-last-$name \
         -R "span[hosts=1]" \
         -e logs/trim-last-$name.err \
         -o logs/trim-last-$name.out -n 8
