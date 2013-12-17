#lastdb -w 2 -u ~/src/last-hg/examples/bisulfite_f.seed $REF.last_f $REF 
#lastdb -w 2 -u ~/src/last-hg/examples/bisulfite_r.seed $REF.last_r $REF

. ./common.sh
OUTDIR=$OUT

echo "
./src/last-bisulfite-paired.sh $REF.last_f $REF.last_r $FQ1 $FQ2 $name \
    | samtools view -bS - \
    | samtools sort - $OUTDIR/last-$name
samtools fixmate $OUTDIR/last-$name.bam $OUTDIR/last-$name.fix.bam
samtools sort $OUTDIR/last-$name.fix.bam $OUTDIR/last-$name
rm $OUTDIR/last-$name.fix.bam
samtools index $OUTDIR/last-$name.bam
" | bsub -J last-$name \
         -e logs/last-$name.err \
         -o logs/last-$name.out -n 8

echo "
./src/last-bisulfite-paired.sh $REF.last_f $REF.last_r $TRIM_FQ1 $TRIM_FQ2 $name \
    | samtools view -bS - \
    | samtools sort - $OUTDIR/trim/last-$name
samtools fixmate $OUTDIR/trim/last-$name.bam $OUTDIR/trim/last-$name.fix.bam
samtools sort $OUTDIR/trim/last-$name.fix.bam $OUTDIR/trim/last-$name
rm $OUTDIR/trim/last-$name.fix.bam
samtools index $OUTDIR/trim/last-$name.bam
" | bsub -J trim-last-$name \
         -e logs/trim-last-$name.err \
         -o logs/trim-last-$name.out -n 8
