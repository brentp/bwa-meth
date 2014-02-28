. ./common.sh

set -x

perl $BSMOOTH/bin/bswc_bowtie2_index.pl --name=$(dirname $REF)/$(basename $REF .fa) $REF
exit;

mkdir -p results/
prog=bsmooth

BAM=results/$prog-$name.bam

rm logs/$prog-$name.err logs/$prog-$name.out
echo "
perl $BSMOOTH/bin/bswc_bowtie2_align.pl \
    --metrics results/bsmooth.metrics.txt \
    --out results/$prog-$name/ \
    --stop-after-alignment \
    --bam $BAM \
    -- /home/brentp/src/bsmooth-align/mm10 -- $REF -- --very-sensitive -p 6 -- $FQ1 -- $FQ2

samtools view -h $BAM.crick.bam | python src/bsmooth-adjust-crick.py | samtools view -bS - > $BAM.crick.fix.bam
java -jar $PICARD/MergeSamFiles.jar I=$BAM.watson.bam I=$BAM.crick.fix.bam O=$BAM AS=false " | bsub -J $prog-$name -e logs/$prog-$name.err -o logs/$prog-$name.out -n 12


