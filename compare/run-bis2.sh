. ./common.sh

method=bis2

mkdir -p $OUT/$method/trim/
TEMP=/tmp/

rm -f $OUT/$method/$name*.bam $OUT/trim/$method/$name*.bam
FLAGS="--gzip -p 8 --score_min L,-0.4,-0.5 --bowtie2 --bam -N 1"

cmd="bismark $FLAGS --temp_dir $TEMP --output_dir $OUT/$method/ --prefix $name $BIS_REF -1 $FQ1 -2 $FQ2;
rm -f $OUT/$method-$name.bam;
samtools sort $OUT/$method/${name}.${name}*.bam $OUT/$method-$name;
samtools index $OUT/$method-$name.bam"

rm -f logs/$method-$name.err logs/$method-$name.out
echo $cmd  | bsub -J $method-$name -e logs/$method-$name.err -o logs/$method-$name.out -n 12 -R "span[hosts=1]"

mkdir -p $TEMP/trim/

cmd="bismark $FLAGS --temp_dir $TEMP/trim/ --output_dir $OUT/$method/trim/ --prefix $name $BIS_REF -1 $TRIM_FQ1 -2 $TRIM_FQ2;
rm -f $OUT/trim/$method-$name.bam;
samtools sort $OUT/$method/trim/${name}.${name}*.bam $OUT/trim/$method-$name;
samtools index $OUT/trim/$method-$name.bam"

rm -f logs/trim-$method-$name.err logs/trim-$method-$name.out
echo $cmd | bsub -J trim-$method-$name -e logs/trim-$method-$name.err -o logs/trim-$method-$name.out -n 12 -R "span[hosts=1]"
