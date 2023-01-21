# Pre-requisite: STAR v2.7.3a
# You will need to build a STAR genome index prior to running this step

# Insert your list of samples and loop through them
for s in 'sample 1' 'sample 2'

do
	echo "Aligning sample "$s"..."
	mkdir "$s"_star_output

# transcripts.gtf is available at https://github.com/ChrisMaherLab/TMECircDB/sample_files
# star_index is the directory your STAR index is built in
# R1.fastq.gz and R2.fastq.gz are your pair-end RNA-Seq raw data files

	STAR --chimSegmentMin 18 --chimOutType WithinBAM Junctions --sjdbGTFfile transcripts.gtf --outReadsUnmapped Fastx --outFileNamePrefix "$s"_star_output/"$s". --genomeDir star_index --readFilesCommand zcat --readFilesIn "$s".R1.fastq.gz "$s".R2.fastq.gz --runThreadN 18 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --genomeLoad NoSharedMemory  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --twopassMode Basic --alignSoftClipAtReferenceEnds No --outSAMattributes NH HI NM MD AS XS --outSAMstrandField intronMotif

done
