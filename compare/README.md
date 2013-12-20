Setup
=====

The base files for the analysis should be in data/
    real_R1.fastq.gz
    real_R2.fastq.gz

You can then generate the trimmed versions as:

    bash src/trim.sh data/real_R{1,2}.fastq.gz 


Simulated
+++++++++

To generate the simulated data, see `src/gen-simulated.sh`
which should require only minimal changes to point to the
reference fasta for `mm10` on your system.


Once the simulated data is generated, run:

    bash src/trim.sh data/sim_R{1,2}.fastq.gz 


Adjust ./common.sh so that `REF` points the the `mm10` on your system.

Index
=====
You will then need to create the appropriate indexes for each of the
alingers you wish to test.

Bowtie:
    bismark_genome_preparation /path/to/reference/files/

GSNAP:
    python src/gsnap-meth.py index $REF

bwa-meth:
    python ../bwa-meth.py index $REF

last:
    lastdb -w 2 -u /path/to/last-hg/examples/bisulfite_f.seed $REF.last_f $REF 
    lastdb -w 2 -u /path/to/last-hg/examples/bisulfite_r.seed $REF.last_r $REF

bsmap:
    no indexing


Align
=====

Once you have the indexes in place and `REF=` in common.sh in place, you can
use the run-{method} scripts to run each of the aligners.

All of the run-{method} scripts source `common.sh` so make sure that has
what you want.

Assessment
==========

For simulated reads:

    $ python src/sim-roc.py \
        --reads data/sim_R1.fastq.gz \ # this is to get the number of input reads
        results/*-sim.bam

For real reads:

    $ python src/target-roc.py \
        --regions data/mm10.capture-regions.bed.gz \
        --reads data/read_R1.fastq.gz \ # this is to get the number of input reads
        results/*-real.bam

Both real and simulated data can also be assesed after trimming.
BAMs from trimmed input appear in results/trim/

Note
====

If you want to simply align some reads. The entire syntax is:

    python bwa-meth.py index /path/to/ref.fasta
    python bwa-meth.py --reference /path/to/ref.fasta reads_R1.fastq reads_R2.fastq -p myoutput

and the result will appear in myoutput.bam

