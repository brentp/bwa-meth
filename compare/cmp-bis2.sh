<<DONE
# get alignments from bis2
samtools view -f2 results/bis2-noerrorsim.bam | cut -f1 | uniq > bis2-noerror.txt
# get alignments from bwa-meth that are not in bis2
samtools view -f2 results/bwa-noerrorsim.bam | awk '$5 > 20' | cut -f 1 \
    | uniq | grep -vFf bis2-noerror.txt | sort -u > bwa-uniq.txt


exit;
DONE

python src/extract-read.py "10001_chr6:125453407-125453496" \
    data/noerrorsim_R1.fastq.gz \
    data/noerrorsim_R2.fastq.gz


bowtie2 -q -N 1 --score-min L,-0.4,-0.5 --reorder --ignore-quals --no-mixed \
    --maxins 600 -x \
    /data/Schwartz/brentp/mm10/ref/Bisulfite_Genome/CT_conversion/BS_CT \
    -1 o_R1.fq -2 o_R2.fq -q | grep -v "^@"

bowtie2 -q -N 1 --score-min L,-0.4,-0.5 --reorder --ignore-quals --no-mixed \
    --maxins 600 -x \
    /data/Schwartz/brentp/mm10/ref/Bisulfite_Genome/GA_conversion/BS_GA \
    -1 o_R1.fq -2 o_R2.fq -q | grep -v "^@"

