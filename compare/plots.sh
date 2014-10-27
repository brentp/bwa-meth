set -ex
<<DONE
python src/target-roc.py data/mm10.capture-regions.bed.gz \
    --reads data/real_R1.fastq.gz \
    results/*-real.bam \
    results/bison-real/real_R1.bam \
    results/bsmooth/bsmooth-real.bam \
    > real-quals.txt

python src/target-roc.py data/mm10.capture-regions.bed.gz \
    --reads data/real_R1.fastq.gz \
    results/trim/*-real.bam \
    results/trim/bison-real/real_R1.trim.bam \
    results/trim/bsmooth/bsmooth-real.bam \
    > real-trim-quals.txt
DONE

for sim in sim noerrorsim; do

    python src/sim-roc.py \
        --reads data/${sim}_R1.fastq.gz \
        results/*-${sim}.bam \
        results/bison-$sim/${sim}_R1.bam \
        results/bsmooth/bsmooth-$sim.bam \
        > ${sim}-quals.txt

    python src/sim-roc.py \
        --reads data/${sim}_R1.fastq.gz \
        results/trim/*-${sim}.bam \
        results/trim/bison-$sim/${sim}_R1.trim.bam \
        results/trim/bsmooth/bsmooth-$sim.bam \
        > ${sim-}trim-quals.txt
done
# make table of time and mem from logs
python src/parse-time-mem.py logs/*-{noerrorsim,sim,real}.out
