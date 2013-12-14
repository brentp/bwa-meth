set -eo pipefail
REF=/data/Schwartz/brentp/mm10/ref/mm10.fa

FQ1=/proj/Schwartz/brentp/2013/ken-rrbs/pilot/38379_CAGATC_L003_R1_001.fastq.gz
FQ2=${FQ1/_R1_/_R2_}

name=s38379

TRIM_FQ1=/proj/Schwartz/brentp/2013/ken-rrbs/pilot/trimmed/38379_R1.fastq.gz
TRIM_FQ2=${TRIM_FQ1/_R1/_R2}

OUTDIR=results/

TEMP=/scratch/brentp/

mkdir -p $OUTDIR/trim/ logs/

