. ./common.sh

#    -n 1000000 \
cd data/
Sherman \
    --genome_folder $BIS_REF \
    --length 100 \
    -n 100000 \
    -pe \
    -X 500 \
    --CG_conversion 25 \
    --CH_conversion 95 \
    --error_rate 1

sed 's/_R1$//' simulated_1.fastq | gzip -c > sim_R1.fastq.gz &
sed 's/_R2$//' simulated_2.fastq | gzip -c > sim_R2.fastq.gz
wait
rm -f simulated_1.fastq simulated_2.fastq
