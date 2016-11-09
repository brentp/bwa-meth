These commands can be run from this directory once `bwa-meth` is installed

1. Index The Reference.
2. Align the Reads.

```Shell

# rm ref.fa.bwameth.* 
bwameth.py index ref.fa
bwameth.py --reference ref.fa t_R1.fastq.gz t_R2.fastq.gz -t 12 | samtools view -b - > bwa-meth.bam

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


From here, it is recommended to use [PileOMeth](https://github.com/dpryan79/PileOMeth) for extraction and tabulation of the methylation.

