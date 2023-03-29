#!/usr/bin/Rscript

# Run this script after you have finished running TMECircDB_CIRCexplorer_sum.r and TMECircDB_featureCounts_sum.r

# Change this to your working directory
# This directory should contain the linear gene counts matrix count.rds and the circRNA counts matrix aggregate_circRNA_counts.txt
# For demonstration purposes, count.rds and aggregate_circRNA_counts.txt are available at https://github.com/ChrisMaherLab/TMECircDB/sample_files
setwd('current_directory')

# Pre-requisites
library(grid)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(gplots)
library(edgeR)
library(devtools)
library(RColorBrewer)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

##### Data preparation #####

# Import sample info
# Sample info file should have the following columns: sample_id (our IDs have the format of "group_subject_sample number"),	subject (subject index),	group (N, P, M, etc),	match_status (whether or not there is a matched sample),	mapped (uniquely mapped reads from STAR)
# crc_matched_patients_read_info.txt is an example of our data, and is available at https://github.com/ChrisMaherLab/TMECircDB/sample_files
group.file = 'crc_matched_patients_read_info.txt';
m = read.table(group.file, header=T, stringsAsFactors=F, sep='\t', quote = '');

# Only use normal, primary, and metastasis samples
groups = c('N','P', 'M');
m = m[m$group %in% groups,];
rownames(m) = m$sample_id;

# Import and clean both linear genes and circRNA expression matrices
## 1. Linear gene
count.file = 'featureCounts/count.rds';
cnt0 = readRDS(count.file);

# Preserve metadata columns, adjust the numbers according to your data
num.meta.cols0 = 4;
meta.cols0 = colnames(cnt0)[1:num.meta.cols0];
meta0 = cnt0[, meta.cols0];

# Import gene annotation file
# transcripts.meta.tsv is available at https://github.com/ChrisMaherLab/TMECircDB/sample_files
annotation<-read.table("transcripts.meta.tsv",header = T,sep='\t',stringsAsFactors = F,quote = "");
annotation_biotype<-unique(annotation[,c("gene","gene.biotype")]);
meta0<-merge(x=meta0,y=annotation_biotype,by.x="gene",by.y="gene",all.x=T);
meta0<-meta0[,c(1,2,5)];
colnames(meta0) <- c("transcript","gene","transcript.biotype");

rownames(meta0) = meta0$transcript;
rownames(cnt0) = rownames(meta0);
cnt0 = cnt0[, !(colnames(cnt0) %in% meta.cols0)];
colnames(cnt0) <- data.frame(strsplit(colnames(cnt0), "\\/"))[11,];
colnames(cnt0) <- data.frame(strsplit(colnames(cnt0), "Aligned"))[1,];

cnt0 <- cnt0[,order(colnames(cnt0))]

colnames(cnt0) <- m$sample_id
cnt0 = cnt0[, colnames(cnt0) %in% m_filter$sample_id]
cnt0$host_gene <- rownames(cnt0)

## 2. CircRNA
count.file1 = 'aggregate_circRNA_counts.txt';
cnt1 = read.table(count.file1, header=T, stringsAsFactors=F, sep='\t', quote='',check.names=FALSE);

# Pre-filtering of expression values
# In our study, we required ≥1 backspliced read in 15% of all samples for circRNAs
cnt1 <- cnt1[rowSums(cnt1[colnames(cnt1)%in% m$sample[m$group %in% groups]]>=1) >= 0.15*nrow(m),];

# Preserve metadata columns, adjust the numbers according to your data
num.meta.cols1 = 3;
meta.cols1 = colnames(cnt1)[1:num.meta.cols1];
meta1 = cnt1[, meta.cols1];

meta1[,2] <- "circRNA"
meta1[,3] <- paste0(data.frame(strsplit(meta1$uniq_identifier, "\\|"))[7,])
meta1[,4] <- paste0(data.frame(strsplit(meta1$uniq_identifier, "\\|"))[8,])
meta1<-meta1[,c(1,3,4,2)];
colnames(meta1) <- c("circRNA","gene","transcript","transcript.biotype")

rownames(meta1) = meta1$circRNA
rownames(cnt1) = rownames(meta1)

cnt1 = cnt1[, !(colnames(cnt1) %in% meta.cols1)]
cnt1 <- cnt1[,order(colnames(cnt1))]

colnames(cnt1) <- m$sample_id
cnt1 = cnt1[, colnames(cnt1) %in% m_filter$sample_id]
cnt1$circRNA <- rownames(cnt1)
cnt1$host_gene <- meta1$transcript

# Merge circRNA and linear gene read count matrices based on host gene, adjust column numbers accordingly
colnames(cnt1)[1:27] <- paste0(colnames(cnt1)[1:27],"_circ")
colnames(cnt0)[1:27] <- paste0(colnames(cnt0)[1:27],"_linear")

cnt_merge <- merge(cnt1,cnt0,by="host_gene",all=T)

Circ <- cnt_merge[,c(29,2:28)]
Linear <- cnt_merge[,c(29:56)]

N_cols <- grep("N_",colnames(Circ))
P_cols <- grep("P_",colnames(Circ))
M_cols <- grep("M_",colnames(Circ))

PvN_test <- Circ.test(Circ[,c(1,N_cols,P_cols)], Linear[,c(1,N_cols,P_cols)], group=c(rep(1,8),rep(2,6)), circle_description = 1)
MvN_test <- Circ.test(Circ[,c(1,N_cols,M_cols)], Linear[,c(1,N_cols,M_cols)], group=c(rep(1,8),rep(2,13)), circle_description = 1)
MvP_test <- Circ.test(Circ[,c(1,M_cols,P_cols)], Linear[,c(1,M_cols,P_cols)], group=c(rep(1,6),rep(2,13)), circle_description = 1)

PvN_results <- cbind(Circ[,c(1,N_cols,P_cols)], Linear[,c(1,N_cols,P_cols)])
PvN_results <- cbind(PvN_results,PvN_test$p.val,PvN_test$p.adj)
PvN_results_sig <- subset(PvN_results,PvN_results$`PvN_test$p.val`<=0.05)

MvN_results <- cbind(Circ[,c(1,N_cols,M_cols)], Linear[,c(1,N_cols,M_cols)])
MvN_results <- cbind(MvN_results,MvN_test$p.val,MvN_test$p.adj)
MvN_results_sig <- subset(MvN_results,MvN_results$`MvN_test$p.val`<=0.05)

MvP_results <- cbind(Circ[,c(1,M_cols,P_cols)], Linear[,c(1,M_cols,P_cols)])
MvP_results <- cbind(MvP_results,MvP_test$p.val,MvP_test$p.adj)
MvP_results_sig <- subset(MvP_results,MvP_results$`MvP_test$p.val`<=0.05)

aggregate_results <- cbind(Circ,Linear[2:28],
                           PvN_test$p.val,PvN_test$p.adj,
                           MvN_test$p.val,MvN_test$p.adj,
                           MvP_test$p.val,MvP_test$p.adj)
aggregate_results$PvN_DE <- aggregate_results$`PvN_test$p.val`<=0.05
aggregate_results$MvN_DE <- aggregate_results$`MvN_test$p.val`<=0.05
aggregate_results$MvP_DE <- aggregate_results$`MvP_test$p.val`<=0.05
aggregate_results_sig <- subset(aggregate_results,aggregate_results$`PvN_test$p.val`<=0.05 |
                                  aggregate_results$`MvN_test$p.val` <=0.05 |
                                  aggregate_results$`MvP_test$p.val` <=0.05)

write.table(aggregate_results,"CircTest_results.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(aggregate_results_sig,"CircTest_results_DE_only.txt",sep = "\t",col.names = T,row.names = F,quote = F)
