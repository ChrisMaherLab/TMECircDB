#!/usr/bin/Rscript

# Run this script after you have downloaded or pre-processed scRNA-Seq data

# Pre-requisites
library(dplyr)
library(tidyverse)

# Change this to your working directory
# This directory should contain a cell annotation file and an expression matrix file
# For demonstration purposes, we have included data downloaded from Lee 2020 with GEO accesion code GSE144735: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735
setwd('current_directory')
annotation <- read.table("GSE144735_processed_KUL3_CRC_10X_annotation.txt",
                         header = T, stringsAsFactors=F, sep='\t', quote='');
log_tpm_matrix <- read.table("GSE144735_processed_KUL3_CRC_10X_natural_log_TPM_matrix.txt",
                             header = T, stringsAsFactors=F, sep='\t', quote='');
cell_types <- unique(annotation$Cell_type);


# Reverse the log(tpm+1) transform
tpm_matrix <- log_tpm_matrix;
tpm_matrix[2:dim(tpm_matrix)[2]] <- exp(log_tpm_matrix[2:dim(log_tpm_matrix)[2]])-1;

colnames(tpm_matrix)[2:dim(tpm_matrix)[2]] <- gsub("\\.","-",colnames(tpm_matrix)[2:dim(tpm_matrix)[2]]);
rownames(tpm_matrix) <- tpm_matrix$Index;
new_tpm_matrix <- tpm_matrix[,2:dim(tpm_matrix)[2]];

# Use all patients for scRNA-seq reference
# You can change which patients you want to include
chosen_patients <- unique(annotation$Patient);

# Make scRNA-seq reference matrix
sampled_cells_tpm <- data.frame(log_tpm_matrix$Index);
ground_truth_reference <- data.frame(cell_types);
ground_truth_reference[,chosen_patients] <-0;

for (i in 1:length(chosen_patients))
{
  print(paste0("......Processing for patient ",chosen_patients[i]));
  chosen_patient <- chosen_patients[i];
  chosen_cells_per_patient <- subset(annotation,annotation$Patient==chosen_patient);
  print(paste0("......Total number of cells: ",dim(chosen_cells_per_patient)[1]));
  chosen_cells_per_patient_tpm <- new_tpm_matrix[,colnames(new_tpm_matrix) %in% chosen_cells_per_patient$Index];
  
  for (j in 1:length(cell_types))
  {
    cell_type <- cell_types[j];
    print(paste0("...Processing for cell subtype No. " ,j,": ",cell_type));
    cell_type_annotation <- subset(chosen_cells_per_patient,chosen_cells_per_patient$Cell_type==cell_type);
    print(paste0("Number of cells in subtype: ", dim(cell_type_annotation)[1]));
    
    # We used 10% of all cells for a down-sample
    # Adjust to the proportion you need
    new_sampled_cells <- slice_sample(cell_type_annotation,prop = 0.1);
    
    if (dim(new_sampled_cells)[1]>0)
    {
      if (dim(new_sampled_cells)[1]==1)
      {
        print("only 1 cells sampled");
        new_sampled_cells_tpm <- data.frame(chosen_cells_per_patient_tpm[,colnames(chosen_cells_per_patient_tpm) %in% new_sampled_cells$Index]);
        colnames(new_sampled_cells_tpm) <- new_sampled_cells$Index[1];
      }
      else
      {
        print(paste0(dim(new_sampled_cells)[1]," cells sampled"));
        new_sampled_cells_tpm <- data.frame(chosen_cells_per_patient_tpm[,colnames(chosen_cells_per_patient_tpm) %in% new_sampled_cells$Index]);
      }
      ground_truth_reference[ground_truth_reference$cell_types==cell_type,
                             chosen_patient] <- dim(new_sampled_cells)[1];
      print("tpm extracted for cells sampled");
      sampled_cells_tpm <- cbind(sampled_cells_tpm,new_sampled_cells_tpm);
      print("cbind complete");
    }
    else
    {
      print("0 cells sampled");
    }
  }
}

colnames(sampled_cells_tpm)[1] <- "Gene";
colnames(sampled_cells_tpm)[2:dim(sampled_cells_tpm)[2]] <- gsub("\\.","-",colnames(sampled_cells_tpm)[2:dim(sampled_cells_tpm)[2]]);

# Output the sampled cells with their cell-tag IDs
write.table(sampled_cells_tpm,paste0("scRNA_Seq_reference_with_ID_GSE144735_full.txt"),sep='\t', quote=F,col.names = T,row.names = F)

# Re-write all cell-tag IDs with corresponding cell types
for (i in 2:dim(sampled_cells_tpm)[2] )
{
  cell_index <- colnames(sampled_cells_tpm)[i];
  colnames(sampled_cells_tpm)[i] <- annotation[annotation$Index==cell_index,"Cell_type"];
}

# Calculate ground truth cell numbers
ground_truth_reference[7,2:7] <- colSums(ground_truth_reference[,2:7]);
ground_truth_reference[7,1] <- "total_no_cells";

# Output sampled cells with their cell types
write.table(sampled_cells_tpm,paste0("scRNA_Seq_reference_GSE144735_full.txt"),sep='\t', quote=F,col.names = T,row.names = F)

# Output ground truth sampled cell numbers
write.table(ground_truth_reference,paste0("scRNA_Seq_reference_ground_truth_GSE144735_full.txt"),sep='\t', quote=F,col.names = T,row.names = F)