bwa-meth
========

align BS-Seq reads and tabulate methylation without intermediate temp files.
This only works for reads from the directional protocol (most common).

Uses the method employed by methylcoder and Bismark of using *in silico*
conversion of all C's to T's in both reference and reads.

Recovers the original read (needed to tabulate methylation) by attaching it
as a comment which bwa appends as a tag to the read.

usage
=====

Index
-----

One time only, you need to index a reference sequence.

    python bwa-meth.py index $REFERENCE

If your reference is `some.fasta`, this will create `some.c2t.fasta`
and all of the bwa indexes associated with it.

Align
-----

    python bwa-meth.py --threads 16 \
         --prefix $PREFIX \
         --reference $REFERENCE \
         $FQ1 $FQ2
         
This will create $PREFIX.bam and $PREFIX.bam.bai. The output will pass
Picard-tools ValidateSam and will have the
reads in the correct location (flipped from G => A reference).

Handles clipped alignments and indels correctly. Fastqs can be gzipped
or not.

The command above will be sent to BWA to do the work as something like:

    bwa mem -L 25 -pCM -t 15  $REFERENCE.c2t.fa \
            '<python bwa-meth.py c2t $FQ1 $FQ2'

So the converted reads are streamed directly to bwa and **never written
to disk**. The output from that is modified by `bwa-meth.py` and streamed
straight to a bam file.

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
You can also align with `bwa-meth.py` with the `--calmd` flag to have samtools
do it's BAQ.


This software may use the samtools SNP calling pipeline at some point
to do something similar to BisSNP.

