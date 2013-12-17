. ./common.sh

REF=/data/Schwartz/brentp/mm10/ref/ # bowtie2
REF=/data/Schwartz/brentp/mm10/ # bowtie

mkdir -p $OUT/bis1/trim/

cmd="bismark --gzip --maxins 1000 -n 2 -l 24 --bam --temp_dir $TEMP --output_dir $OUT/bis1/ --prefix $name $REF -1 $FQ1 -2 $FQ2;
samtools sort $OUT/bis1/*.bam $OUT/bis1-$name;
samtools index $OUT/bis1-$name.bam"

echo $cmd  | bsub -J bis1-$name -e logs/bis1-$name.err -o logs/bis1-$name.out -n 2



cmd="bismark --gzip -maxins 1000 -n 2 -l 24 --bam --temp_dir $TEMP --output_dir $OUT/bis1/trim/ --prefix $name $REF -1 $TRIM_FQ1 -2 $TRIM_FQ2;
samtools sort $OUT/bis1/trim/*.bam $OUT/trim/bis1-$name;
samtools index $OUT/trim/bis1-$name.bam"

echo $cmd | bsub -J trim-bis1-$name -e logs/trim-bis1-$name.err -o logs/trim-bis1-$name.out -n 2
