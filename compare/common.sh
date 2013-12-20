REF=/data/Schwartz/brentp/mm10/ref/mm10.fa

set -eo pipefail

name=real
name=sim


FQ1=data/${DATA}_R1.fastq.gz
FQ2=${FQ1/_R1/_R2}


TRIM_FQ1=${FQ1/.fastq/.trim.fastq}
TRIM_FQ2=${FQ2/.fastq/.trim.fastq}


OUT=results/

TEMP=/scratch/brentp/

mkdir -p $OUT/trim/ logs/


