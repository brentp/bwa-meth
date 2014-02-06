REF=/data/Schwartz/brentp/mm10/ref/mm10.fa
TEMP=/scratch/brentp/

PATH=$PATH:~/src/sherman/

set -eo pipefail

name=real
name=sim


FQ1=data/${name}_R1.fastq.gz
FQ2=${FQ1/_R1/_R2}


TRIM_FQ1=${FQ1/.fastq/.trim.fastq}
TRIM_FQ2=${FQ2/.fastq/.trim.fastq}

OUT=results/


mkdir -p $OUT/trim/ logs/ data/

#if [ ! -e $TRIM_FQ1 ]; then
#    exit 1
#fi
