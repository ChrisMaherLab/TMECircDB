# Run this script after you have finished running TMECircDB_STAR_align.sh
# Pre-requisite: CIRCexplorer v2.3.8

# Insert your list of samples and loop through them
for s in 'sample 1' 'sample 2'

do

        echo "Parsing chimeric alignments for sample "$s"..."
        mkdir "$s"_circexplorer2_output
        cd "$s"_circexplorer2_output
        CIRCexplorer2 parse -t STAR "$s"_star_output/"$s".Chimeric.out.junction > "$s"_circexplorer2_output/"$s"_CIRCexplorer2_parse.log

done
