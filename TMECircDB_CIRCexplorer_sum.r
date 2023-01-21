#!/usr/bin/Rscript

# Change this to your working directory
setwd('current_directory')

# Import gene annotation file
# transcripts.meta.tsv is available at https://github.com/ChrisMaherLab/TMECircDB/sample_files

annotation<-read.table("transcripts.meta.tsv",header = T,sep='\t',stringsAsFactors = F,quote = "");
annotation_biotype<-unique(annotation[,c("transcript","transcript.biotype")]);

# Import sample info
# Sample info file should have the following columns: sample_id (our IDs have the format of "group_subject_sample number"),	subject (subject index),	group (N, P, M, etc),	match_status (whether or not there is a matched sample),	mapped (uniquely mapped reads from STAR)
# crc_matched_patients_read_info.txt is an example of our data, and is available at https://github.com/ChrisMaherLab/TMECircDB/sample_files
sampleInfo<-read.table("crc_matched_patients_read_info.txt",header=T,sep='\t',stringsAsFactors = F,quote = "");

#Insert your list of samples and loop through them
sample_ids <-c('sample 1', 'sample 2');


for (i in 1:length(sample_ids))
{
  sample_id <- sample_ids[i];
  print(sample_id);
  
  # Read CIRCexplorer out file and add headers (Please see https://circexplorer2.readthedocs.io/en/latest/modules/annotate for more details)
  circRNA<-read.table(paste0(sample_id,"_circexplorer2_output",sample_id,".circularRNA_known.txt"),header=T,sep='\t',stringsAsFactors = F,quote = "");
  colnames(circRNA)<-c("chrom","start","end","name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes","exonOffsets","readNumber","circType","geneName","isoformName","index","flankIntron");
  
  # Filter out alternative haplotypes
  circRNA <- subset(circRNA,circRNA$chrom
                 %in% c("chr1","chr2","chr3","chr4","chr5",
                        "chr6","chr7","chr8","chr9","chr10",
                        "chr11","chr12","chr13","chr14","chr15",
                        "chr16","chr17","chr18","chr19","chr20",
                        "chr21","chr22","chrX","chrY","chrM"));
  
  # Here the "biotype" indicates the biotype of the parental gene of a circRNA
  circRNA<-merge(x=circRNA,y=annotation_biotype,by.x="isoformName",by.y="transcript",all.x=T);
  
  # Calculate normalized backspliced read counts
  circRNA$readNumberNormalized<- circRNA$readNumber/sampleInfo[sampleInfo$sample_id==sample_id,"mapped"]*1000000;
  
  # Generate unique ID for each circRNA
  cols <- c("chrom","start","end","strand","exonCount","circType","geneName","isoformName","index");
  circRNA$uniq_identifier <- apply( circRNA[ , cols ] , 1 , paste , collapse = "|" );
  circRNA$uniq_identifier <- gsub(" ","",circRNA$uniq_identifier,fixed=TRUE);
  
  # Make raw counts dataframe and normalized counts dataframe
  readNum<-data.frame(circRNA$uniq_identifier,circRNA$transcript.biotype,circRNA$readNumber);
  readNumNormalized<-data.frame(circRNA$uniq_identifier,circRNA$transcript.biotype,circRNA$readNumberNormalized);
  colnames(readNum)<-c("uniq_identifier","transcript.biotype",paste0("readNumber_",sample_id));
  colnames(readNumNormalized)<-c("uniq_identifier","transcript.biotype",paste0("readNumberNormalized_",sample_id));
  
  #Assign variable names to this sample's output dataframes
  assign(paste0("circRNA_counts_",sample_id),readNum);
  assign(paste0("circRNA_norm_counts_",sample_id),readNumNormalized);
  
  }

# Merge all individual dataframes into aggregated expression matrices for raw counts and normalized counts
# You will need to edit the list() according to your samples
merged_counts <- Reduce(function(x,y) merge(x = x, y = y, by = c("uniq_identifier","transcript.biotype"), all = T), 
                         list(`circRNA_counts_sample 1`, 
                              `circRNA_counts_sample 2`));
merged_norm_counts <- Reduce(function(x,y) merge(x = x, y = y, by = c("uniq_identifier","transcript.biotype"), all = T), 
                        list(`circRNA_norm_counts_sample 1`, 
                             `circRNA_norm_counts_sample 2`));

# Counte the number of samples with expression for each circRNA
merged_counts$numSamples <- rowSums(!is.na(merged_counts[,grepl('readNumber', colnames(merged_counts))]),na.rm=T);
merged_norm_counts$numSamples <- rowSums(!is.na(merged_norm_counts[,grepl('readNumberNormalized', colnames(merged_norm_counts))]),na.rm=T);


write.table(merged_counts,"aggregate_circRNA_counts.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(merged_norm_counts,"aggregate_circRNA_normalized_counts.txt",sep = "\t",col.names = T,row.names = F,quote = F)
