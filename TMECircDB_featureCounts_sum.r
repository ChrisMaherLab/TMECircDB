#!/usr/bin/Rscript

# Run this script after you have finished running TMECircDB_featureCounts.sh

# Pre-requisites
library(ggplot2)
library(reshape2)

# Change this to your working directory
setwd('current_directory')

files = list.files('featureCounts', '.cnt$')

# Aggreate read counts for genes across samples
y = NULL
for (f in files){
    x = read.table(file.path('featureCounts', f), header=T, 
        stringsAsFactors=F, sep='\t', check.names=F)
    x = x[, c(1,6,7)]
    samp = gsub('(^bams.|.bam$)', '', colnames(x)[3])
    colnames(x) = c('gene', 'length', samp)
    if (is.null(y)){y = x}else{y = merge(y,x)}
}

# Read annotation gene meta data
# gene2name.tsv is available at https://github.com/ChrisMaherLab/TMECircDB/sample_files
a = read.table('sample_files/gene2name.tsv', header=T, sep='\t', quote='', stringsAsFactors=F,fill=T)
a = a[, c('gene', 'symbol','biotype')]

# Attach gene meta-data to count matrix
y = merge(a,y)
saveRDS(y, file='count.rds')

# Keep gene meta-data (now with gene length in col 4)
a = y[, 1:4]

# Calculate TPM
cnt = y[, -(1:4)]
rownames(cnt) = y$gene
len = as.matrix(y[, rep('length', ncol(cnt))])
cnt = as.matrix(cnt)
colnames(len) = colnames(cnt)
rpk = cnt*1000/len
m = colSums(rpk)/1000000
m = matrix(m, nrow=1, byrow=T)
m = m[rep(1, nrow(rpk)),]
tpm = rpk/m
tpm = as.data.frame.matrix(tpm)
tpm = cbind(gene=rownames(tpm), tpm)
tpm = merge(a, tpm)

saveRDS(tpm, file='tpm.rds')


#stop()

# Create expression file in tablular format
tpmm = tpm; tpmm$length=NULL
tpmm = melt(tpmm)
colnames(tpmm) = c('gene', 'symbol','biotype', 'sample', 'tpm')
tpmm$gene = as.character(tpmm$gene)
tpmm$sample = as.character(tpmm$sample)
saveRDS(tpmm, file='tpm.tab.rds')

# Summary of featureCounts
files = list.files('featureCounts', '.summary$')
z = NULL
for (f in files){
    x = read.table(file.path('featureCounts', f), header=T, stringsAsFactors=F, sep='\t', check.names=F)
    if (is.null(z)){z = x}else{z = merge(z,x)}
}

colnames(z) = gsub('(^bams/|.bam$)', '', colnames(z))
rownames(z) = z$Status; z$Status = NULL
z = as.matrix(z)
x = melt(z)
colnames(x) = c('metric', 'sample', 'count')
x$metric = as.character(x$metric)
x$sample = as.character(x$sample)

pdf('featureCounts.qc.pdf', width=8, height=11, title='')
ggplot(x[x$count > 0,], aes(x=sample, y=count, fill=metric)) + geom_bar(stat='identity') + coord_flip()
dev.off()

write.table(x, file='featureCounts.summary.tsv', sep='\t', row.names=F, quote=F)
