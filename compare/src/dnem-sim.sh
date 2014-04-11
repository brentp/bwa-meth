# wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/snp137Common.txt.gz
DNM=/home/brentp/src/dnemulator/dnemulator-16
PATH=$PATH:$DNM
set -e
H=$(pwd)

REF=/data/Schwartz/brentp/mm10/ref/mm10.fa
R1=$H/data/real_R1.fastq.gz
R2=$H/data/real_R2.fastq.gz

mkdir -p data/simmed/
cd data/simmed

rm -f ./*
fasta-methyl-sim $REF > mm10.meth.fa
fasta-polymorph $DNM/snp137Common.txt.gz mm10.meth.fa mm10.meth.fa > mm10.poly.fa

fasta-paired-chunks -n 1000000 -l 100 mm10.poly.fa sim_R1.fa sim_R2.fa
fasta-bisulf-sim sim_R1.fa > sim.R1.bs.fa
fasta-bisulf-sim sim_R2.fa > sim.R2.bs.fa
fastq-sim sim.R1.bs.fa $R1 | awk 'NR % 4 == 1 { gsub(/ /, ":"); print  }(NR % 4 != 1)' > sim_R1_bs.fastq
fastq-sim sim.R2.bs.fa $R2 | awk 'NR % 4 == 1 { gsub(/ /, ":"); print  }(NR % 4 != 1)' > sim_R2_bs.fastq

python ../../src/fix-names.py sim_R1_bs.fastq sim_R2_bs.fastq

gzip -f sim_R1.fastq
gzip -f sim_R2.fastq
mv sim_R1.fastq.gz ../dnemsim_R1.fastq.gz
mv sim_R2.fastq.gz ../dnemsim_R2.fastq.gz

