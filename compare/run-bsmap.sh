. ./common.sh

cmd="bsmap -u -a $FQ1 -b $FQ2 -d $REF -o results/bsmap-$name.bam -v 3 \
             -p 8 -m 0 -x 1000 -S 42 -n 0 -s 12 -I 1"
# run bsmap with default parameters
cmd="bsmap -a $FQ1 -b $FQ2 -d $REF -o results/defaultbsmap-$name.bam -p 8 "

rm -f logs/bsmap-$name.err logs/bsmap-$name.out
echo $cmd | bsub -J bsmap-$name \
                 -R "span[hosts=1]" \
                 -e logs/bsmap-$name.err \
                 -o logs/bsmap-$name.out -n 8


cmd="bsmap -u -a $TRIM_FQ1 -b $TRIM_FQ2 -d $REF -o results/trim/bsmap-$name.bam -v 3 \
             -p 8 -m 0 -x 1000 -S 42 -n 0 -s 12 -I 1"

# run bsmap with default parameters
cmd="bsmap -a $TRIM_FQ1 -b $TRIM_FQ2 -d $REF -o results/trim/defaultbsmap-$name.bam -p 8 "

rm -f logs/trim-bsmap-$name.err logs/trim-bsmap-$name.out
echo $cmd | bsub -J trim-bsmap-$name \
                 -R "span[hosts=1]" \
                 -e logs/trim-bsmap-$name.err \
                 -o logs/trim-bsmap-$name.out -n 8
