. ./common.sh
cd data/
Sherman \
    --genome_folder /data/Schwartz/brentp/mm10/ \
    --length 100 \
    -n 1000000 \
    -pe \
    -X 800 \
    --CG_conversion 25 \
    --CH_conversion 95 \
    --error_rate 0.01 

mv simulated_1.fastq sim_R1.fastq && gzip -f sim_R1.fastq
mv simulated_2.fastq sim_R2.fastq && gzip -f sim_R2.fastq
