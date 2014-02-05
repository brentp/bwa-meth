bwa-meth
========

align BS-Seq reads and tabulate methylation without intermediate temp files.
This only works for **paired-end reads from the directional protocol** 
(most common).

Uses the method employed by methylcoder and Bismark of using *in silico*
conversion of all C's to T's in both reference and reads.

Recovers the original read (needed to tabulate methylation) by attaching it
as a comment which **bwa** appends as a tag to the read.

Performs favorably to existing aligners gauged by number of on and off-target reads for a capture method that targets CpG-rich region. Some off-target regions may be enriched, but all aligners are be subject to the same assumptions.
See manuscript: http://arxiv.org/abs/1401.1129 for details.
Optimal alignment is the upper-left corner. Curves are drawn by varying the
mapping quality cutoff for alingers that use it.


Vertical dotted line is mapping quality of 60 for bwa.

Run.sh scripts for each method are here: https://github.com/brentp/bwa-meth/tree/master/compare
I have done my best to have each method perform optimally, but no doubt there
could be improvements.

QuickStart
==========

Without installation, you can use as `python bwameth.py` with install, the
command is `bwameth.py`.

The commands:

    bwameth.py index $REF
    bwameth.py --reference $REF some_R1.fastq.gz some_R2.fastq.gz --prefix some.output

will create `some.output.bam` and `some.output.bam.bai`

See the **full example** at: https://github.com/brentp/bwa-meth/tree/master/example/

Installation
============

`bwa-meth` depends on 

 + python 2.7+ (including python3)
   - `toolshed` library. can be installed with: 
      * `easy_install toolshed` or
      * `pip install toolshed`

   - or run: `python setup.py install` from this directory

 + samtools command on the `$PATH` (https://github.com/samtools/samtools)

 + bwa mem from: https://github.com/lh3/bwa


usage
=====

Index
-----

One time only, you need to index a reference sequence.

    bwameth.py index $REFERENCE

If your reference is `some.fasta`, this will create `some.c2t.fasta`
and all of the bwa indexes associated with it.

Align
-----

    bwameth.py --threads 16 \
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
            '<python bwameth.py c2t $FQ1 $FQ2'

So the converted reads are streamed directly to bwa and **never written
to disk**. The output from that is modified by `bwa-meth` and streamed
straight to a bam file.

Tabulate
--------

Currently, `bwa-meth` calls Bis-SNP to call methylation for CpG's and genotypes 
for SNPs.

E.g.:

    bwameth.py tabulate \
                --trim 3,3 \
                --map-q 60 \
                --bissnp BisSNP-0.82.2.jar \
                --prefix out \
                -t 12 \
                --reference $REF \
                $BAM1 $BAM2 ... $BAMN

This will use BisSNP to perform multi-sample SNP and CpG calling to create
`out.cpg.vcf` and `out.snp.vcf` as well as a BED file of methylation for
each input BAM.

The `--trim` removes the first and last 3 bases from reads to avoid bias at
read ends.
