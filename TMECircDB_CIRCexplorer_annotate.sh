# Run this script after you have finished running TMECircDB_CIRCexplorer_parse.sh
# There should be a back_spliced_junction.bed intermediate file in the output folder of each of your samples

# Insert your list of samples and loop through them
for s in 'sample 1' 'sample 2'

do
        echo $s

# hg38_ref_all.txt is available at https://github.com/ChrisMaherLab/TMECircDB/sample_files
# all-chrs.fa is the reference genome FATSA you used to build your STAR genome index

        CIRCexplorer2 annotate -r hg38_ref_all.txt -g all-chrs.fa -b "$s"_circexplorer2_output/back_spliced_junction.bed -o "$s"_circexplorer2_output/"$s".circularRNA_known.txt > "$s"_circexplorer2_output/"$s".CIRCexplorer2_annotate.log

done
