#lastdb -w 2 -u ~/src/last-hg/examples/bisulfite_f.seed $REF.last_f $REF 
#lastdb -w 2 -u ~/src/last-hg/examples/bisulfite_r.seed $REF.last_r $REF

source ./common.sh

echo "
./src/last-bisulfite-paired.sh $REF.last_f $REF.last_r $FQ1 $FQ2 $name \
    | samtools view -bS - \
    | samtools sort -n@3 - $OUTDIR/$name
samtools fixmate $OUTDIR/$name.bam $OUTDIR/$name.fix.bam
samtools sort -@3 $OUTDIR/$name.fix.bam $OUTDIR/$name
rm $OUTDIR/$name.fix.bam
samtools index $OUTDIR/$name.bam
" | bsub -J last-$name \
         -e logs/last-$name.err \
         -o logs/last-$name.out -n 8

echo "
./src/last-bisulfite-paired.sh $REF.last_f $REF.last_r $TRIM_FQ1 $TRIM_FQ2 $name \
    | samtools view -bS - \
    | samtools sort -n@3 - $OUTDIR/trim/$name
samtools fixmate $OUTDIR/trim/$name.bam $OUTDIR/trim/$name.fix.bam
samtools sort -@3 $OUTDIR/trim/$name.fix.bam $OUTDIR/trim/$name
rm $OUTDIR/trim/$name.fix.bam
samtools index $OUTDIR/trim/$name.bam
" | bsub -J trim-last-$name \
         -e logs/trim-last-$name.err \
         -o logs/trim-last-$name.out -n 8
