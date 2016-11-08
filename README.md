bwa-meth
========

Fast and accurante alignment of BS-Seq reads.

## NOTE!!!

As of 2016-08-18, bwa-meth now outputs sam to stdout. It is up to the user to
convert to bam. This means that the --prefix and --calmd flags are gone.

## Update 2016

`bwa-meth` is still among (if not *the*) best aligners for BS-Seq.
While it is fairly stable, I will continue to support the alignment
part of `bwa-meth`--fixing any bugs or updating as needed.

There are now several (likely better) alternatives for tabulation
and SNP calling than provided here so I will not develop those further.

For tabulation, bias, and plotting, use [PileOMeth](https://github.com/dpryan79/PileOMeth)

For SNP calling (a more modern BisSNP), use [biscuit](https://github.com/zwdzwd/biscuit)


## Intro

This works for single-end reads and for **paired-end reads from the
directional protocol** (most common).

Uses the method employed by methylcoder and Bismark of *in silico*
conversion of all C's to T's in both reference and reads.

Recovers the original read (needed to tabulate methylation) by attaching it
as a comment which **bwa** appends as a tag to the read.

Performs favorably to existing aligners gauged by number of on and off-target reads for a capture method that targets CpG-rich region. Some off-target regions may be enriched, but all aligners are be subject to the same assumptions.
See manuscript: http://arxiv.org/abs/1401.1129 for details.
Optimal alignment is the upper-left corner. Curves are drawn by varying the
mapping quality cutoff for alingers that use it.

This image is on real reads and represents an attempt to find good parameters
for all aligners tested.

![Untrimmed reads comparison](https://gist.githubusercontent.com/brentp/bf7d3c3d3f23cc319ed8/raw/d5f1ebcc53b924a05a5980159bfcb97494ec34f2/real.gif)

Note that *bwa-meth* and *Last* perform well without trimming.

run.sh scripts for each method are here: https://github.com/brentp/bwa-meth/tree/master/compare
I have done my best to have each method perform optimally, but no doubt there
could be improvements.

QuickStart
==========

Without installation, you can use as `python bwameth.py` with install, the
command is `bwameth.py`.

The commands:

    bwameth.py index $REF
    bwameth.py --reference $REF some_R1.fastq.gz some_R2.fastq.gz > some.output.sam

will create `some.output.bam` and `some.output.bam.bai`.
To align single end-reads, specify only 1 file.

See the **full example** at: https://github.com/brentp/bwa-meth/tree/master/example/

Installation
============

The following snippet should work for most systems that have samtools and bwa
installed and the ability to install python packages. (Or, you can send this
to your sys-admin). See the dependencies section below for further instructions: 

```Shell

    # these 4 lines are only needed if you don't have toolshed installed
    wget https://pypi.python.org/packages/source/t/toolshed/toolshed-0.4.0.tar.gz
    tar xzvf toolshed-0.4.0.tar.gz
    cd toolshed-0.4.0
    sudo python setup.py install

    wget https://github.com/brentp/bwa-meth/archive/v0.10.tar.gz
    tar xzvf v0.10.tar.gz
    cd bwa-meth-0.10/
    sudo python setup.py install

```

After this, you should be able to run: `bwameth.py` and see the help.

Dependencies
------------

`bwa-meth` depends on 

 + python 2.7+ (including python3)
   - `toolshed` library. can be installed with: 
      * `easy_install toolshed` or
      * `pip install toolshed`

   - if you don't have root or sudo priviledges, you can run
     `python setup.py install --user` from this directory and the bwameth.py
     executable will be at: ~/.local/bin/bwameth.py

   - if you do have root or sudo run: `[sudo] python setup.py install` from
     this directory

   - users unaccustomed to installing their own python packages should 
     download anaconda: https://store.continuum.io/cshop/anaconda/ and
     then install the toolshed module with pip as described above.

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
         --reference $REFERENCE \
         $FQ1 $FQ2 > some.sam
         
The output will pass will have the reads in the correct location (flipped from G => A reference).

Handles clipped alignments and indels correctly. Fastqs can be gzipped
or not.

The command above will be sent to BWA to do the work as something like:

    bwa mem -L 25 -pCM -t 15  $REFERENCE.c2t.fa \
            '<python bwameth.py c2t $FQ1 $FQ2'

So the converted reads are streamed directly to bwa and **never written
to disk**. The output from that is modified by `bwa-meth` and streamed
straight to a bam file.
