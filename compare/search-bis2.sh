. ./common.sh

OUT=results/bismark2-search/
mkdir -p $OUT/logs

FLAGS="--gzip -p 8 --bowtie2 --bam "
name=default
cmd="bismark $FLAGS --temp_dir $TEMP --output_dir $OUT/$name/ --prefix $name $BIS_REF -1 $FQ1 -2 $FQ2;"
mkdir -p $OUT/$name/trim/
rm -f $OUT/$name/*.bam $OUT/trim/$name/*.bam
rm -f $OUT/logs/$name.err $OUT/logs/$name.out
echo $cmd  | bsub -J $name -e $OUT/logs/$name.err \
                                   -o $OUT/logs/$name.out -n 12 -R "span[hosts=1]"

cmd="bismark $FLAGS --temp_dir $TEMP/trim/ --output_dir $OUT/$name/trim/ --prefix $name $BIS_REF -1 $TRIM_FQ1 -2 $TRIM_FQ2;"
rm -f $OUT/logs/trim-$name.err $OUT/logs/trim-$name.out
echo $cmd | bsub -J trim-$name -e $OUT/logs/trim-$name.err -o $OUT/logs/trim-$name.out -n 12 -R "span[hosts=1]"

for L in 18 20 22; do
    for s0 in -0.5 -0.25  0; do
        for s1 in -0.5  -0.25  0; do

            TEMP=/tmp/
            name="L-$L-score_${s0}_$s1"
            mkdir -p $OUT/$name/trim/
            rm -f $OUT/$name/*.bam $OUT/trim/$name/*.bam

            FLAGS="--gzip -p 8 --score_min L,$s0,$s1 --bowtie2 --bam -N 1 -L $L"

            cmd="bismark $FLAGS --temp_dir $TEMP --output_dir $OUT/$name/ --prefix $name $BIS_REF -1 $FQ1 -2 $FQ2;"
            rm -f $OUT/logs/$name.err $OUT/logs/$name.out
            echo $cmd  | bsub -J $name -e $OUT/logs/$name.err \
                                               -o $OUT/logs/$name.out -n 12 -R "span[hosts=1]"

            mkdir -p $TEMP/trim/

            cmd="bismark $FLAGS --temp_dir $TEMP/trim/ --output_dir $OUT/$name/trim/ --prefix $name $BIS_REF -1 $TRIM_FQ1 -2 $TRIM_FQ2;"
            rm -f $OUT/logs/trim-$name.err $OUT/logs/trim-$name.out
            echo $cmd | bsub -J trim-$name -e $OUT/logs/trim-$name.err -o $OUT/logs/trim-$name.out -n 12 -R "span[hosts=1]"
        done
    done
done
