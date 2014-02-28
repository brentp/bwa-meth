set -e

R1=$1
R2=$2


O1=$(dirname $1)/$(basename $1 .fastq.gz).trim.fastq.gz;
O2=$(dirname $2)/$(basename $2 .fastq.gz).trim.fastq.gz;

trim_galore -t --paired --quality 15 $R1 $R2

mv *_R1_val_1.fq.gz $O1
mv *_R2_val_2.fq.gz $O2
