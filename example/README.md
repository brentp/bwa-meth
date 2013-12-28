These commands can be run from this directory once `bwa-meth` is installed

1. Index The Reference.
2. Align the Reads.

```Shell

# rm ref.c2t* ref.dict ref.fa.fai
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


Then tabulate the methylation assuming you have the BisSNP.jar file at $BISSNP

```Shell
    bwameth.py tabulate \
             --trim 3,3 \
             --map-q 60 --bissnp $BISSNP \
             --prefix ex \
             -t 12 \
             --reference ref.fa \
             bwa-meth.bam

```
This will create `ex`.cpg.vcf and `ex`.snp.vcf along with a .BED file for each bam:

```Shell
$ head ext_R.cpg.bed 
```

    #CHROM	START	POS	PCT	CS	TS	CTX
    chrREF	507	508	100	1	0	CG
    chrREF	524	525	0	0	1	C
    chrREF	566	567	0	0	1	C
    chrREF	582	583	0	0	1	CG
    chrREF	623	624	0	0	1	CG
    chrREF	1782	1783	100	1	0	CG
    chrREF	1791	1792	0	0	1	CG
    chrREF	1798	1799	0	0	1	C
    chrREF	1808	1809	100	1	0	CG
