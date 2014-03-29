suppressPackageStartupMessages(library(bsseq))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(TRUE)

cov=data.frame(Sample=c(paste("2", 1:8, sep="_"), paste("1", 1:8, sep="_")),
               grp_treat=as.factor(c(rep("tol", 8), rep("nontol", 8))))

#cov = read.delim("data/Mouse_allergen_clinical_data.txt")
#cov = cov[tolower(cov$Methylation) == "yes",]
#cov$grp_treat = paste(cov$Group, cov$Treatment, sep="_")
#cov = cov[cov$grp_treat %in% c('MTHFR-/-_HDM', 'Control_HDM'),]
#cov$grp_treat = relevel(cov$grp_treat, ref="Control_HDM")


rownames(cov) = cov$Sample

read.bs = function(fname, sample){
    dat = as.data.frame(fread(fname))
    chr = as.matrix(dat[,1], drop=FALSE)
    pos = as.matrix(dat[,3], drop=FALSE)

    M = as.matrix(dat$n_same, drop=FALSE)
    C = as.matrix(dat$n_same + dat$n_converted, drop=FALSE)
    BSseq(chr=chr, pos=pos, M=M,
               C=C, sampleNames=sample)
}


bobjs = list()

for(i in 1:nrow(cov)){
    sample = cov[i, "Sample"]
    f = paste(sample, "methylation.txt", sep=".")

    bobjs[i] = read.bs(f, sample)
}

bs = do.call("combine", unlist(bobjs))

bs.smooth = BSmooth(bs, mc.cores=8, verbose=TRUE, maxGap=4000)
bs.cov = getCoverage(bs.smooth)

keep = rowMeans(bs.cov) > 2
keep = rowSums(which(bs.cov > 2)) > 2
keep = rowSums(bs.cov >= 2) >= 2
bs.smooth = bs.smooth[keep,]

bs.tstat = BSmooth.tstat(bs.smooth, group1=as.character(cov$Sample[1:8]),
                                    group2=as.character(cov$Sample[9:16]), 
                                    estimate.var="group2",
                                    local.correct=TRUE, mc.cores=4)

dmrs0 = dmrFinder(bs.tstat, cutoff = c(-5, 5), maxGap=200, qcutoff=c(0.01, 0.99))
dmrs = subset(dmrs0, n >= 5 & abs(meanDiff) >= 0.1)

pData(bs.smooth)$col = ifelse(cov$grp_treat == "tol", "red", "blue")

pdf('regions.pdf')
plotManyRegions(bs.smooth, dmrs, extend=100, addRegions=dmrs, BSseqTstat=bs.tstat, stat.ylim=c(-7, 7), addPoints=TRUE)
dev.off()

