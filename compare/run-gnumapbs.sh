. ./common.sh

# taken from gnumap-bs paper
cmd="gnumap -b -m 17 -s 1 -T 20 -a 0.90 -b -G $REF --illumina ???? -o"

rm -f logs/gnumap-$name.err logs/gnumap-$name.out
echo "$cmd results/gnumapbs-$name.sam" | bsub -J gnumap-$name \
                 -e logs/gnumap-$name.err \
                 -o logs/gnumap-$name.out -n 8


rm -f logs/trim-gnumap-$name.err logs/trim-gnumap-$name.out
echo "$cmd results/gnumapbs-$name.sam" | bsub -J trim-gnumap-$name \
                 -e logs/trim-gnumap-$name.err \
                 -o logs/trim-gnumap-$name.out -n 8
