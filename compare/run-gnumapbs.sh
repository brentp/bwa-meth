. ./common.sh

module load gnumap

# args taken from gnumap-bs paper
cmd="gnumaps.pl --genome $REF --lib_type wt1 --acc 0.9 --nt_conv bs --num_threads 1 --outdir"

zless $FQ1 > $$.r1.tmp
zless $FQ2 > $$.r2.tmp

rm -f logs/gnumap-$name.err logs/gnumap-$name.out
echo "$cmd results/gnumapbs-$name/ --pair_1 $$.r1.tmp --pair_2 $$.r2.tmp" \
    | bsub -J gnumap-$name \
                 -R "rusage[mem=42000]" \
                 -e logs/gnumap-$name.err \
                 -o logs/gnumap-$name.out -n 8

exit;

rm -f logs/trim-gnumap-$name.err logs/trim-gnumap-$name.out
echo "$cmd results/trim/gnumapbs-$name/ --pair_1 $$.r1.tmp --pair_2 $$.r2.tmp" \
    | bsub -J trim-gnumap-$name \
                 -e logs/trim-gnumap-$name.err \
                 -o logs/trim-gnumap-$name.out -n 8

