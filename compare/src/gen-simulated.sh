. ./common.sh

cd data/

N=1000000

Sherman \
    --genome_folder $BIS_REF \
    --length 100 \
    -n $N \
    -pe \
    -X 500 \
    --CG_conversion 25 \
    --CH_conversion 95 \
    --error_rate 1.0

sed 's/_R1$//' simulated_1.fastq | gzip -c > sim_R1.fastq.gz &
sed 's/_R2$//' simulated_2.fastq | gzip -c > sim_R2.fastq.gz
wait

# wget -O data/e_coli_k12_mg1655.fa "http://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3?report=fasta&log$=seqview&format=text"
# and remove html

rm -f simulated_1.fastq simulated_2.fastq

Sherman \
    --genome_folder . \
    --length 100 \
    -n $((N / 20)) \
    -pe \
    -X 500 \
    --CG_conversion 98 \
    --CH_conversion 98 \
    --error_rate 1.0

sed 's/_R1$/_BAD/' simulated_1.fastq | gzip -c >> sim_R1.fastq.gz &
sed 's/_R2$/_BAD/' simulated_2.fastq | gzip -c >> sim_R2.fastq.gz
wait
rm -f simulated_1.fastq simulated_2.fastq


Sherman \
    --genome_folder $BIS_REF \
    --length 100 \
    -n $N \
    -pe \
    -X 500 \
    --CG_conversion 25 \
    --CH_conversion 95 \
    --error_rate 0.0

sed 's/_R1$//' simulated_1.fastq | gzip -c > noerrorsim_R1.fastq.gz &
sed 's/_R2$//' simulated_2.fastq | gzip -c > noerrorsim_R2.fastq.gz
wait

rm -f simulated_1.fastq simulated_2.fastq
