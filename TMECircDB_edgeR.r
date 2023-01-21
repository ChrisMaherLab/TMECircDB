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
meta1<-meta1[,c(1,3,2)];
colnames(meta1) <- c("transcript","gene","transcript.biotype");

# Here we are changing all biotypes for circRNA transcripts to "circRNA" for the heatmaps
meta1[,3] <- "circRNA";
meta1[,2] <- paste0("circ",data.frame(strsplit(meta1$transcript, "\\|"))[7,]);

rownames(meta1) = meta1$transcript;
rownames(cnt1) = rownames(meta1);
cnt1 = cnt1[, !(colnames(cnt1) %in% meta.cols1)];

# Generate a shorter display name for circRNA in heatmap
meta0$display_name <- paste0(meta0$gene,"|",meta0$transcript);
meta1$display_name <- paste0(meta1$gene,"|",data.frame(strsplit(meta1$transcript, "\\|"))[8,],"|",data.frame(strsplit(meta1$transcript, "\\|"))[9,]);

## 3. Merge linear gene and circRNA
cnt = rbind(cnt0,cnt1);
cnt = cnt[,colnames(cnt) %in% m$sample];
cnt = cnt[,order(colnames(cnt))];
colnames(cnt) = m[order(m$sample),'sample_id'];
cnt = as.matrix(cnt);
meta = rbind(meta0,meta1);

# Paired test or not paired test?
paired = FALSE;

# Prepare output directory
out.dir = 'edgeR_output';
dir.create(out.dir);

# Keep samples shared by group file and count matrix only and label samples based on groups
shared.samples = intersect(rownames(m), colnames(cnt));
cnt = cnt[, shared.samples];
m = m[rownames(m) %in% shared.samples,];

colnames(cnt) = paste0(colnames(cnt), '.counts');
mets.samples = m$sample_id[m$group == 'M'];
tumor.samples = m$sample_id[m$group == 'P'];
normal.samples = m$sample_id[m$group == 'N'];

# Generate average counts for each group
mean.cnt = data.frame(
  N.mean.cnt = rowMeans(cnt[, paste0(normal.samples, '.counts')]),
  P.mean.cnt = rowMeans(cnt[, paste0(tumor.samples, '.counts')]),
  M.mean.cnt = rowMeans(cnt[, paste0(mets.samples, '.counts')])
);

##### EdgeR starting #####


# Create edgeR expression object
y = DGEList(counts=cnt, genes=rownames(cnt));

# Normalize counts
y = calcNormFactors(y);
cpm = cpm(y);

# CPM/TPM filter of linear gene expression values
# In our study, we required ≥5 TPM in 25% of all samples for linear genes
min.cnt.linear = 5; min.cnt.num.samples.linear = 0.25*nrow(m);
cpm.filter.pass.linear = rowSums(cpm[rownames(cpm)%in%rownames(meta0),]>=min.cnt.linear) >= min.cnt.num.samples.linear;

keep = cpm.filter.pass.linear; 
table(keep);
y = y[keep, ];

# After filtering, need to renormalize and recalculate CPM
y <- calcNormFactors(y);
meta = as.data.frame(meta[rownames(meta) %in% rownames(y),]);


# Create grouping factor
grps = factor(as.character(m$group), levels=groups);
subjects = factor(as.character(m$subject));

# Design matrix
if(paired){
  design <- model.matrix(~subjects+grps);
}else{
  design = model.matrix(~grps);
}

rownames(design) <- colnames(y);
print(design);

# Run edgeR test
y = estimateGLMCommonDisp(y, design, verbose=TRUE);
y = estimateGLMTrendedDisp(y, design);
y = estimateGLMTagwiseDisp(y, design, trend=F);
fit = glmFit(y, design);


# Recalculate CPM and edgeR RPM
cat('Generating CPM and edgeR rpm\n');
cpm = cpm(y, log=F);
write.table(cpm, file=paste0(out.dir, '/cpm.tsv'), sep='\t', quote=F);
colnames(cpm) = paste0(colnames(cpm), '.CPM');
edgeR.rpm = cpm(y, log=F);
write.table(edgeR.rpm, file=paste0(out.dir, '/rpm.edgeR.tsv'), sep='\t', quote=F);
colnames(edgeR.rpm) = paste0(colnames(edgeR.rpm), '.edgeR.RPM');
mean.edgeR.rpm = NULL;
for (grp in groups){
  x = data.frame(mean.rpm=rowMeans(edgeR.rpm[, paste0(m$sample_id[m$group == grp], '.counts.edgeR.RPM')]));
  colnames(x) = paste0(grp, '.mean.edgeR.rpm');
  if (is.null(mean.edgeR.rpm)){mean.edgeR.rpm = x}else{mean.edgeR.rpm = cbind(mean.edgeR.rpm, x)};
}

res.all = NULL;
n.groups=length(groups);
for (i in 1:(n.groups - 1)){
  for (j in (i+1):n.groups){
    grp1 = groups[i];
    grp2 = groups[j];
    cmp = paste0(grp2, 'v', grp1);
    grp1 = paste0('grps', grp1);
    grp2 = paste0('grps', grp2);
    cat('*****\t', cmp, '\t*****\n');
    con = rep(0, ncol(design));
    con[grepl(grp1, colnames(design))] = -1;
    con[grepl(grp2, colnames(design))] = 1;
    names(con) = colnames(design);
    cat('contrast:\n'); print(con);
    lrt <- glmLRT(fit, contrast=con);
    res = topTags(lrt, n=nrow(lrt))[[1]];
    res.meta = meta[res$genes,];
    
    cat('Number of test with pvalue = 1:', sum(res$PValue == 1), '. If this is too high, consider pre/post-filtering.\n');

    # Output dataframe
    # In our study, we required FDR adjusted p≤0.15
    # You can adjust the stringency of DE criteria here, including p, adjusted p, and fold change, etc
    res = cbind(res.meta, res, mean.cnt[res$genes,], mean.edgeR.rpm[res$genes,]);
    res$status = 'unchanged';
    res$status[res$PValue <= 0.05 & res$logFC > 0 & res$FDR<=0.15] = 'upregulated';
    res$status[res$PValue <= 0.05 & res$logFC < 0 & res$FDR<=0.15] = 'downregulated';
  
    cmp.cols = c('logFC', 'LR', 'PValue', 'FDR', 'status');
    res = res[, c(setdiff(colnames(res), cmp.cols), cmp.cols)];
    colnames(res)[colnames(res) %in% cmp.cols] = paste0(cmp, '.', colnames(res)[colnames(res) %in% cmp.cols]);
    
    if (is.null(res.all)){res.all = res}else{res.all = merge(res.all, res, all=T)};
    
    fdr = paste0(cmp, '.FDR');
    PVALUE = paste0(cmp, '.PValue');
    stat = paste0(cmp, '.status');
    
    # Volcano plot
    logFC = paste0(cmp, '.logFC');
    pval = paste0(cmp, '.PValue');

    # Tag circRNAs
    res$is.circRNA = FALSE;
    res$is.circRNA[res$transcript.biotype=="circRNA"] = TRUE;
    
    known.idx = res$is.circRNA;
    known.x = res[known.idx, logFC];
    known.p = res[known.idx, pval];
    known.y = -log10(known.p);
    res$logFC = res[[logFC]];
    res$PValue = res[[pval]];
    res$status = res[[stat]];
    p = (ggplot(res, aes(x=logFC, y=-log10(PValue), color=status))
         + geom_point() + theme_bw() + ggtitle(cmp)
         + annotate('text', x=known.x, y=known.y, label='·', color='cyan',size=14)
    );
    ggsave(p, file=paste0(out.dir, '/', cmp, '.volcano.png'),width=10,height=10);
    print(summary(de <- decideTestsDGE(lrt, p=.001)));
  }
}


# Tag lncRNAs & circRNAs
lncrna.patterns = '3prime_overlapping_ncrna|Maherlab_novel|antisense|lincRNA|tucp|miRNA|misc_RNA|processed_transcript|refseq_rna|rRNA|sense_intronic|sense_overlapping|snRNA|snoRNA';
res.all$is.lncRNA = FALSE;
res.all$is.lncRNA[grepl(lncrna.patterns, res.all$transcript.biotype) & !grepl('protein', res.all$transcript.biotype)] = TRUE;
res.all$is.circRNA = FALSE;
res.all$is.circRNA[res.all$transcript.biotype=="circRNA"] = TRUE;

# Write output
## All results
write.table(res.all, file=paste0(out.dir, '/results.tsv'), row.names=F, sep='\t', quote=F);

## Only DE genes and circRNAs
select = rowSums(res.all[,grep("FDR",colnames(res.all))]<=0.15,na.rm=T)>=1 & rowSums(res.all[,grep("status",colnames(res.all))]!="unchanged",na.rm=T)>=1;
res.sig <- res.all[select,];
write.table(res.sig, file=paste0(out.dir, '/results.fdr-0.15.tsv'), row.names=F, sep='\t', quote=F);

## Only DE circRNAs
res.circ <- subset(res.all,res.all$transcript.biotype=="circRNA");
res.circ.sig <- subset(res.sig,res.sig$transcript.biotype=="circRNA");
write.table(res.circ.sig, file=paste0(out.dir, '/results.fdr-0.15-circ.tsv'), row.names=F, sep='\t', quote=F);

##### Heatmap of top DE genes #####
top.de.transcirpts = res.sig;
top.de.circ = subset(top.de.transcirpts,top.de.transcirpts$transcript.biotype=="circRNA");
top.de.linear = subset(top.de.transcirpts,top.de.transcirpts$transcript.biotype!="circRNA");

e = log2(cpm+1);
e=e[order(rownames(e)),order(colnames(e))];
meta = meta[order(rownames(meta)),];
rownames(e) = meta$display_name;

colnames(e) = gsub('.counts.CPM', '', colnames(e));

e_col_N <- grep("N",colnames(e));
e_col_P <- grep("P",colnames(e));
e_col_M <-grep("M",colnames(e));
e <- e[,c(e_col_N,e_col_P,e_col_M)];
 
## 1. All DE genes & circRNAs heatmap
e_all <- e[top.de.transcirpts$display_name,];
groupCols <- ifelse(m[colnames(e_all), 'group'] == 'N', 'purple',
                    ifelse(m[colnames(e_all), 'group'] == 'P','orange',
                           'magenta'));
matchStatusCols <- ifelse(m[colnames(e_all), 'match_status'] %in% c('N','Y_no_normal','Y_no_primary'), 
                          'pink', 
                          'gray');
clab=cbind(groupCols,matchStatusCols);
colnames(clab)=c("Tissue type","Match status");
transcriptTypeRows = t(as.matrix(ifelse(res.all[res.all$display_name %in%rownames(e_all), 'is.lncRNA'] == 'TRUE', 'cyan',
                                        ifelse(res.all[res.all$display_name %in%rownames(e_all), 'is.circRNA'] == 'TRUE', 'yellow',       
                                               'black'))));
rownames(transcriptTypeRows) <- "Gene biotype";


colBreaks = c(-2, seq(-1,1,0.1), 2);
col = colorpanel(length(colBreaks)-1, low='blue', mid='white', high='red');
ngenes = nrow(e_all);


png(paste0(out.dir, '/heatmap_all.png'), res=600, units='in', width=11+ncol(e_all)*0.150, height=18);
h = heatmap.3(e_all, 
              scale='row',
              labRow = FALSE,
              #labCol = FALSE,
              col=col, breaks=colBreaks, 
              trace='none',
              dendrogram="none", cexCol=1, cexRow=0.9, 
              key=T, keysize=2, 
              main=paste0("Top DE Genes (Linear & CircRNA) in Mets vs Primary vs Normal\nn=", ngenes),
              Rowv=T, Colv=F,
              margins=c(10,10),
              hclust=function(x) hclust(x,method="complete"),
              distfun=function(x) as.dist((1-cor(t(x)))/2),
              ColSideColors=clab,
              RowSideColors=transcriptTypeRows,
              ColSideColorsSize=2, 
              RowSideColorsSize=1,
              colsep = c(length(e_col_N),length(e_col_N)+length(e_col_P)),
              sepcolor = "white",
              symbreaks=FALSE,
              symkey=FALSE,
              lhei=c(1,10), lwid=c(2,5)
);

legend("left", legend=c("Normal","Primary","Mets",
                        "","Matched","Unmatched",
                        "","Protein coding","LncRNA","CircRNA"), 
       fill=c("purple","orange","magenta",
              "white","pink","gray",
              "white","black","cyan","yellow"), 
       border=FALSE,
       cex=1, box.lty=0,
       inset=-0.005);

dev.off()

## 2. Only DE linear genes heatmap
e_linear <- e[top.de.linear$display_name,];
groupCols <- ifelse(m[colnames(e_linear), 'group'] == 'N', 'purple',
                    ifelse(m[colnames(e_linear), 'group'] == 'P','orange',
                           'magenta'));
matchStatusCols <- ifelse(m[colnames(e_linear), 'match_status'] %in% c('N','Y_no_normal','Y_no_primary'), 
                          'pink', 
                          'gray');
clab=cbind(groupCols,matchStatusCols);
colnames(clab)=c("Tissue type","Match status");
transcriptTypeRows = t(as.matrix(ifelse(res.all[res.all$display_name %in%rownames(e_linear), 'is.lncRNA'] == 'TRUE', 'cyan',
                                        ifelse(res.all[res.all$display_name %in%rownames(e_linear), 'is.circRNA'] == 'TRUE', 'yellow',       
                                        'black'))));
rownames(transcriptTypeRows) <- "Gene biotype";


colBreaks = c(-2, seq(-1,1,0.1), 2);
col = colorpanel(length(colBreaks)-1, low='blue', mid='white', high='red');
ngenes = nrow(e_linear);


png(paste0(out.dir, '/heatmap_linear.png'), res=600, units='in', width=11+ncol(e_linear)*0.150, height=18);
h = heatmap.3(e_linear, 
              scale='row',
              labRow = FALSE,
              #labCol = FALSE,
              col=col, breaks=colBreaks, 
              trace='none',
              dendrogram="none", cexCol=1, cexRow=0.9, 
              key=T, keysize=2, 
              main=paste0("Top DE Linear Genes in Mets vs Primary vs Normal\nn=", ngenes),
              Rowv=T, Colv=F,
              margins=c(10,10),
              hclust=function(x) hclust(x,method="complete"),
              distfun=function(x) as.dist((1-cor(t(x)))/2),
              ColSideColors=clab,
              RowSideColors=transcriptTypeRows,
              ColSideColorsSize=2, 
              RowSideColorsSize=1,
              colsep = c(length(e_col_N),length(e_col_N)+length(e_col_P)),
              #rowsep = 20,
              sepcolor = "white",
              symbreaks=FALSE,
              symkey=FALSE,
              lhei=c(1,10), lwid=c(2,5)
);

legend("left", legend=c("Normal","Primary","Mets",
                        "","Matched","Unmatched",
                        "","Protein coding","LncRNA"), 
       fill=c("purple","orange","magenta",
              "white","pink","gray",
              "white","black","cyan"), 
       border=FALSE,
       cex=1, box.lty=0,
       inset=-0.005);

dev.off()

## 2. Only DE circRNAs heatmap
e_circ <- e[top.de.circ$display_name,];
groupCols <- ifelse(m[colnames(e_circ), 'group'] == 'N', 'purple',
                    ifelse(m[colnames(e_circ), 'group'] == 'P','orange',
                           'magenta'));
matchStatusCols <- ifelse(m[colnames(e_circ), 'match_status'] %in% c('N','Y_no_normal','Y_no_primary'), 
                          'pink', 
                          'gray');
clab=cbind(groupCols,matchStatusCols);
colnames(clab)=c("Tissue type","Match status");

colBreaks = c(-2, seq(-1,1,0.1), 2);
col = colorpanel(length(colBreaks)-1, low='blue', mid='white', high='red');
ngenes = nrow(e_circ);

png(paste0(out.dir, '/heatmap_circ.png'), res=600, units='in', width=10, height=18);
h = heatmap.3(e_circ, 
              scale='row',
              #labRow = FALSE,
              #labCol = FALSE,
              col=col, breaks=colBreaks, 
              trace='none',
              dendrogram="none", cexCol=1, cexRow=0.3, 
              key=T, keysize=2, 
              main=paste0("Top DE CircRNAs in Mets vs Primary vs Normal\nn=", ngenes),
              Rowv=T, Colv=F,
              margins=c(10,20),
              # hclust=function(x) hclust(x,method="complete"),
              # distfun=function(x) as.dist((1-cor(t(x)))/2),
              distfun = function(x) as.dist(1-cor(t(x))),
              hclustfun = function(x) hclust(x, method="average"),
              ColSideColors=clab,
              #RowSideColors=transcriptTypeRows,
              ColSideColorsSize=2, 
              #RowSideColorsSize=1,
              colsep = c(length(e_col_N),length(e_col_N)+length(e_col_P)),
              #rowsep = 20,
              sepcolor = "white",
              symbreaks=FALSE,
              symkey=FALSE,
              lhei=c(1,10), lwid=c(2,5)
);

legend("left", legend=c("Normal","Primary","Mets",
                        "","Matched","Unmatched"),
       fill=c("purple","orange","magenta",
              "white","pink","gray"),
       border=FALSE,
       cex=2, box.lty=0,
       inset=0);

dev.off()

