. ./common.sh

REF=/data/Schwartz/brentp/mm10/ # bowtie2

method=bis2

mkdir -p $OUT/$method/trim/

rm -f $OUT/$method/$name*.bam $OUT/trim/$method/$name*.bam
FLAGS="--score_min L,0,-0.3 --bowtie2 --maxins 1000 --bam"

cmd="bismark $FLAGS --temp_dir $TEMP --output_dir $OUT/$method/ --prefix $name $REF -1 $FQ1 -2 $FQ2;
samtools sort $OUT/$method/$name*.bam $OUT/$method-$name;
samtools index $OUT/$method-$name.bam"

rm -f logs/$method-$name.err logs/$method-$name.out
echo $cmd  | bsub -J $method-$name -e logs/$method-$name.err -o logs/$method-$name.out -n 2

mkdir -p $TEMP/trim/

cmd="bismark $FLAGS --temp_dir $TEMP/trim/ --output_dir $OUT/$method/trim/ --prefix $name $REF -1 $TRIM_FQ1 -2 $TRIM_FQ2;
samtools sort $OUT/$method/trim/$name*.bam $OUT/trim/$method-$name;
samtools index $OUT/trim/$method-$name.bam"

rm -f logs/trim-$method-$name.err logs/trim-$method-$name.out
echo $cmd | bsub -J trim-$method-$name -e logs/trim-$method-$name.err -o logs/trim-$method-$name.out -n 2
