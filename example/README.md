These commands can be run from this directory once `bwa-meth` is installed

1. Index The Reference.
2. Align the Reads.

```Shell

# rm ref.fa.bwameth.* 
bwameth.py index ref.fa
bwameth.py --reference ref.fa t_R1.fastq.gz t_R2.fastq.gz -t 12

```

Then check the alignments:

```Shell
samtools flagstat bwa-meth.bam
```

    92791 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 duplicates
    92723 + 0 mapped (99.93%:nan%)
    92791 + 0 paired in sequencing
    46399 + 0 read1
    46392 + 0 read2
    92276 + 0 properly paired (99.44%:nan%)
    92652 + 0 with itself and mate mapped
    71 + 0 singletons (0.08%:nan%)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)


We can create a bias-plot if the `matplotlib` and `seaborn` python modules
are installed.::

    $ python ../bias-plot.py bwa-meth.bam ref.fa
    wrote to bwa-meth.bias.txt
    saving to bwa-meth.bias.png

The text otuput looks like this::

    $ head -6 bwa-meth.bias.txt
    base	read_1	read_2
    1	0.014	0.535
    2	0.344	0.569
    3	0.689	0.523
    4	0.488	0.635
    5	0.516	0.560

We can see a strong bias in the firs couple bases of the first read.
This is even more apparent in the image.

![bias-plot](https://gist.githubusercontent.com/brentp/bf7d3c3d3f23cc319ed8/raw/729f26a3d3c3ef87c026f534979deafd8ab26900/bias-example.png "Bias Plot")

This plot will become smoother in the middle with larger datasets.

One can then tabulate the methylation using the Bis-SNP.jar file at $BISSNP
and taking the bias into account by adjusting the `--trim` argument.

**NOTE**: you will need to run `samtools faidx ref.fa` before this will work as
Bis-SNP assumes there is an index.

```Shell
    bwameth.py tabulate \
             --trim 3,3 \
             --map-q 60 --bissnp $BISSNP \
             --prefix ex \
             -t 12 \
             --reference ref.fa \
             bwa-meth.bam

```
This will create `ex`.cpg.vcf and `ex`.snp.vcf along with a .bed.gz file for each bam:

```Shell
$ head ext_R.cpg.bed
```


    #chrom  start   start   pct cs  ts  ctx
    chrREF  5482    5482    0.0 0   1   CG
    chrREF  5540    5540    100.0   1   0   CG
    chrREF  14606   14606   100.0   1   0   CG
    chrREF  23207   23207   0.0 0   2   CG
    chrREF  23219   23219   0.0 0   2   CG
    chrREF  23260   23260   0.0 0   2   CG
    chrREF  23269   23269   0.0 0   1   CG
    chrREF  23288   23288   0.0 0   2   CG
    chrREF  23323   23323   0.0 0   1   CG

