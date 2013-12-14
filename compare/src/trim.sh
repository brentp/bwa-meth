files=(/proj/Schwartz/brentp/2013/ken-rrbs/pilot/38372_ACAGTG_L003_R1_001.fastq.gz
/proj/Schwartz/brentp/2013/ken-rrbs/pilot/38378_TGACCA_L003_R1_001.fastq.gz
/proj/Schwartz/brentp/2013/ken-rrbs/pilot/38379_CAGATC_L003_R1_001.fastq.gz
/proj/Schwartz/brentp/2013/ken-rrbs/pilot/38381_GCCAAT_L003_R1_001.fastq.gz
/proj/Schwartz/brentp/2013/ken-rrbs/pilot/38387_CGATGT_L003_R1_001.fastq.gz
/proj/Schwartz/brentp/2013/ken-rrbs/pilot/38484_CTTGTA_L003_R1_001.fastq.gz)


TRIM_DIR=/proj/Schwartz/brentp/2013/ken-rrbs/pilot/trimmed/
mkdir -p $TRIM_DIR
#files=(/proj/Schwartz/brentp/2013/ken-rrbs/pilot/38379_CAGATC_L003_R1_001.fastq.gz)
for R1 in "${files[@]}"; do
    R2=${R1/_R1_/_R2_}
    name=$(basename $R1 | perl -pe 's/(\d+).*/$1/')
    sickle pe -f $R1 -r $R2 -o $TRIM_DIR/${name}_R1.fastq \
                            -p $TRIM_DIR/${name}_R2.fastq \
                            -s $TRIM_DIR/${name}_singles.fastq \
                            -t sanger
    for suff in _R1 _R2 _singles; do
        gzip $TRIM_DIR/${name}${suff}.fastq
    done
done
