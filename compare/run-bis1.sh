. ./common.sh

REF=/data/Schwartz/brentp/mm10/ # bowtie
REF=/data/Schwartz/brentp/mm10/ref/ # bowtie2

mkdir -p $OUT/bis1/trim/

cmd="bismark --gzip --bam --temp_dir $TEMP --output_dir $OUT/bis1/ --prefix $name $REF -1 $FQ1 -2 $FQ2;
samtools sort $OUT/bis1/*.bam $OUT/bis-$name; samtools index $OUT/bis-$name.bam"

echo $cmd  | bsub -J bis1-$name -e logs/bis1-$name.err -o logs/bis1-$name.out -n 2



cmd="bismark --gzip --bam --temp_dir $TEMP --output_dir $OUT/bis1/trim/ --prefix $name $REF -1 $TRIM_FQ1 -2 $TRIM_FQ2;
samtools sort $OUT/bis1/trim/*.bam $OUT/trim/bis-$name; samtools index $OUT/trim/bis-$name.bam"

echo $cmd | bsub -J trim-bis1-$name -e logs/trim-bis1-$name.err -o logs/trim-bis1-$name.out -n 2
