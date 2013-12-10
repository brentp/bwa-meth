bwa-meth
========

align BS-Seq reads and tabulate methylation without intermediate temp files.
This only works for reads from the directional protocol (most common).

usage
=====

Index
-----

One time only, you need to index a reference sequence.

    python bwa-meth.py index $REFERENCE

If you reference is some.fasta, this will create some.c2t.fasta
and all of the bwa indexes associated with it.

Align
-----

    python bwa-meth.py --threads 16 \
         --prefix $PREFIX \
         --reference $REFERENCE \
         $READ1 $READ2
         
This will create $PREFIX.bam and $PREFIX.bam.bai the output will pass
Picard-tools ValidateSam.

Tabulate
--------

In progress. Currently, the best way is to run the align step above for each
bam and then use Bis-SNP to call variants or methylation for all BAMs together.

E.g.:

    java -Xmx15g -jar BisSNP-0.82.2.jar \
        -R $REFERENCE \
        -I $PREFIX.bam \  # though best to use multiple BAM's.
        -T BisulfiteGenotyper  \
        --trim_3_end_bp 3 \
        --trim_5_end_bp 3 \
        -vfn1 $PREFIX.cpg.vcf -vfn2 $PREFIX.snp.vcf \
        -stand_call_conf 10 -mmq 20 \
        --dbsnp $DBSNP \
        -nt 16

That will give both CpG's and SNPs. BisSNP recommends performing
the usual GATK preprocessing steps before this (BSQR, indel realignment).

This software may use the samtools SNP calling pipeline at some point
to do something similar.

