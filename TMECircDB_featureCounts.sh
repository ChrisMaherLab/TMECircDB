# Run this script after you have finished running TMECircDB_STAR_align.sh
# Pre-requisite: featureCounts v1.6.0

# Make featureCounts output directory
mkdir featureCounts


# Insert your list of samples and loop through them
for s in 'sample 1' 'sample 2'

do
    echo $s

# transcripts.gtf is available at https://github.com/ChrisMaherLab/TMECircDB/sample_files

    featureCounts --tmpDir /tmp -p -t exon -g gene_id -T 1 -s 0 -M -O -d 40 -D 1000 -Q 1 -a transcripts.gtf -o featureCounts/"$s".cnt "$s"_star_output/"$s".Aligned.sortedByCoord.out.bam"

done
