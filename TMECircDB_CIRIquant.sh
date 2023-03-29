# Run this script after you have finished running TMECircDB_CIRCexplorer_annotate.sh
# Pre-requisites: BWA v0.7.17, HISAT2 v2.2.0, and StringTie v2.1.1

# Make featureCounts output directory
mkdir featureCounts


# Insert your list of samples and loop through them
for s in 'sample 1' 'sample 2'

do
    echo "Quantifying circRNAs in sample "$s"..."
    mkdir "$s"_quantification_output
    #config.yml is set up based on instructions at https://ciri-cookbook.readthedocs.io/en/latest/CIRIquant_2_quantification.html
    #CIRCexplorer2 outputs were generated from TMECircDB_CIRCexplorer_annotate.sh
    CIRIquant -t 16 -1 "$s".R1.fastq.gz -2 "$s".R2.fastq.gz -v --config config.yml -e "$s"_quantification_output/"$s".log -o "$s"_quantification_output -p "$s_" --circ "$s"_circexplorer2_output/"$s".circularRNA_known.txt --tool CIRCexplorer2

done


#Next, run prep_CIRIquant prepDE.py based on instructions at https://ciri-cookbook.readthedocs.io/en/latest/CIRIquant_4_de.html so each comparison group is ready to be analyzed

##NvP

prep_CIRIquant -i CIRIquant_NvP.lst --lib library_info_NvP.csv --circ circRNA_info_NvP.csv --bsj circRNA_bsj_NvP.csv --ratio circRNA_ratio_NvP.csv

prepDE.py -i CIRIquant_DE_NvP.lst -g gene_count_matrix_NvP.csv

CIRI_DE_replicate --lib library_info_NvP.csv --bsj circRNA_bsj_NvP.csv --gene gene_count_matrix_NvP.csv --out circRNA_de_NvP.tsv

##NvM

prep_CIRIquant -i CIRIquant_NvM.lst --lib library_info_NvM.csv --circ circRNA_info_NvM.csv --bsj circRNA_bsj_NvM.csv --ratio circRNA_ratio_NvM.csv

prepDE.py -i CIRIquant_DE_NvM.lst -g gene_count_matrix_NvM.csv

CIRI_DE_replicate --lib library_info_NvM.csv --bsj circRNA_bsj_NvM.csv --gene gene_count_matrix_NvM.csv --out  circRNA_de_NvM.tsv

##PvM

prep_CIRIquant -i CIRIquant_PvM.lst --lib library_info_PvM.csv --circ circRNA_info_PvM.csv --bsj circRNA_bsj_PvM.csv --ratio circRNA_ratio_PvM.csv

prepDE.py -i CIRIquant_DE_PvM.lst -g gene_count_matrix_PvM.csv

CIRI_DE_replicate --lib library_info_PvM.csv --bsj circRNA_bsj_PvM.csv --gene gene_count_matrix_PvM.csv --out circRNA_de_PvM.tsv