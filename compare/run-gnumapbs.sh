. ./common.sh

module load gnumap/3.0.2

export GNUMAPBS=/home/brentp/src/gnumap/gnumaps
PATH=$PATH:$GNUMAPBS/bin:$GNUMAPBS/scripts/
export PATH
# taken from gnumap-bs paper
cmd="gnumaps.pl --genome $REF --lib_type wt1 --acc 0.9 --nt_conv bs --num_threads 8 --outdir"

zless $FQ1 > $$.r1.tmp
zless $FQ2 > $$.r2.tmp

rm -f logs/gnumap-$name.err logs/gnumap-$name.out
echo "$cmd results/gnumapbs-$name/ --pair_1 $$.r1.tmp --pair_2 $$.r2.tmp" \
    | bsub -J gnumap-$name \
                 -e logs/gnumap-$name.err \
                 -o logs/gnumap-$name.out -n 8

exit;

rm -f logs/trim-gnumap-$name.err logs/trim-gnumap-$name.out
echo "$cmd results/trim/gnumapbs-$name/ --pair_1 <(zless $TRIM_FQ1) --pair_2 <(zless $TRIM_FQ2)" \
    | bsub -J trim-gnumap-$name \
                 -e logs/trim-gnumap-$name.err \
                 -o logs/trim-gnumap-$name.out -n 8

