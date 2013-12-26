set -e

R1=$1
R2=$2


O1=$(dirname $1)/$(basename $1 .fastq.gz).trim.fastq;
O2=$(dirname $2)/$(basename $2 .fastq.gz).trim.fastq;

sickle pe -f $R1 -r $R2 \
          -o $O1 -p $O2 \
          -s singles.fastq \
          -t sanger
gzip -f $O1 &
gzip -f $O2 &
rm singles.fastq
wait
