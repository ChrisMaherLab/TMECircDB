# TMECircDB
TMECircDB (Tumor MicroEnvironment specific CircRNA DataBase, https://www.maherlab.com/tmecircdb-overview) is a integrative computational pipeline that leverages published scRNA-Seq datasets (Gene Expression Omnibus accession `GSE144735` and `GSE81861`) and our own unique matched metastatic colorectal cancer (mCRC) patient cohort (Gene Expression Omnibus accession `GSE221240`) to model the cell-type specific expression of circular RNAs (circRNA) in mCRC. The resource has 11,682 circRNAs predicted to have significant cell-type specific expression, including 667 that are exclusively expressed in one cell type. 361 circRNAs showed dysregulation in mCRC proressgion and are hence termed `Circular RNAs Associated with Metastasis (CRAMS)`. We hope this transcriptomic and statistical analysis will serve as a resource for evaluating the cell-type specificity circRNAs that could aide future mechanistic studies exploring their function in cancer. Here we provide scripts and steps used to generate the data.  

## RNA-Seq alignment and circRNA detection
1. Aligh FASTQ files with STAR

```sh TMECircDB_STAR_align.sh```

2. Parse chimeric STAR alignments to identify backsplice junctions with CIRCexplorer

```sh TMECircDB_CIRCexplorer_parse.sh```

3. Annotate backsplice junctions using known transcripts with CIRCexplorer

(1) Individual samples: ```sh TMECircDB_CIRCexplorer_annotate.sh```

(2) Summarize all samples and make circRNA expression matrix: ```Rscript TMECircDB_CIRCexplorer_sum.r```

4. Quantify linear gene reads with featureCounts

(1) Individual samples: ```sh TMECircDB_featureCounts.sh```

(2) Summarize all samples and make linear gene expression matrix: ```Rscript TMECircDB_featureCounts_sum.r```  

## Differential expression analysis & heatmap

```Rscript TMECircDB_edgeR.r```

## Cell-type specific expression modeling
1. Download or pre-process scRNA-Seq expression matrices

Your working directory should contain a cell annotation file and an expression matrix file.
For demonstration purposes, please see GEO accesion code `GSE144735`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735

2. Down-sample scRNA-Seq data and make CIBERSORT reference

```Rscript TMECircDB_scRNA_cibersort_reference.r```

CIBERSORT has a size limit for the scRNA-Seq reference, so it is advisable to select a random down-sample if your scRNA-Seq data has too many cells. In our manuscript, we used 10% of all cells, while keeping the cell-type proportions consistent with the ground truth.

For demonstration purposes, we have included all outputs from this script and our run of the CIBERSORT deconvolution at https://github.com/ChrisMaherLab/TMECircDB/sample_files.

For all other questions related to CIBERSORT, please refer to CIBERSORT tutorial: https://cibersortx.stanford.edu/tutorial.php.

3. Run CIBERSORT's `Create Signature Matrix` module to identify signature genes

Use scRNA_Seq_reference generated in step 2. The following paremeters were used in our study:

`Single cell reference matrix file: scRNA_Seq_reference_GSE144735_full.txt`
`Disable quantile normalization: true`
`kappa: 999`
`q-value: 0.01`
`No. barcode genes: 300 to 500`
`Min. Expression: 1`
`Replicates: 5`
`Sampling: 0.5`
`Filter non-hematopoietic genes from signature matrix during construction: false`

4. Make CIBERSORT mixture file

Use `tpm.rds` generated from `TMECircDB_featureCounts_sum.r` to make the mixture file according to CIBERSORT's instructions.

5. Run CIBERSORT's `Impute Cell Fractions` module to deconvolve cell type proportions

Use signature matrix generated in step 3. The following paremeters were used in our study:

`Signature matrix file: scRNA_Seq_reference_GSE144735_full_inferred_phenoclasses.inferred_refsample.bm.K999.txt`
`Mixture file: crc_matched_patients_linear_tpm_mixture.txt`
`Disable quantile normalization: true`

6. Run Non-Negative Least Squares (NNLS) model to estimate cell-type specific expression & benchmark NNLS predictions against ground truth for each sample

For demonstration purposes, we have included NNLS predictions and benchmarking results at https://github.com/ChrisMaherLab/TMECircDB/sample_files.

(1) For circRNAs: ```Rscript TMECircDB_NNLS_circ.r```

(2) For linear genes: ```Rscript TMECircDB_NNLS_linear_gene.r```

## How to cite TMECircDB?  
Coming soon...  

## Contact  
sidizhao (at) wustl (dot) edu  
