#!/bin/bash

# Define the input and output directories
for i in SRX23809085 SRX23809086 SRX23809075 SRX23809074 SRX23809076 SRX23809077 SRX23809078 SRX23809080 SRX23809082 SRX23809083; do

INPUT_DIR="/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${i}_bam"
OUTPUT_DIR="/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${i}_bam_tagged"
MERGED_BAM="/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${i}_bam_tagged/merged_tagged.bam"

# Create an output directory for tagged BAM files
mkdir -p "$OUTPUT_DIR"

# Function to tag a single BAM file
tag_bam_file() {
    BAM_FILE="$1"
    SRR_NAME=$(basename "$BAM_FILE" | sed -E 's/Aligned.sortedByCoord.out.bam//')
    TAGGED_BAM="$OUTPUT_DIR/${SRR_NAME}_tagged.bam"
    
    echo "Processing $BAM_FILE -> $TAGGED_BAM with CB:Z:$SRR_NAME"
    samtools view -h "$BAM_FILE" | \
    awk -v CB="CB:Z:$SRR_NAME" 'BEGIN {OFS="\t"} {if($1 ~ /^@/) {print $0} else {print $0, CB}}' | \
    samtools view -bo "$TAGGED_BAM"
}

# Export the function and variables for GNU Parallel
export -f tag_bam_file
export OUTPUT_DIR

# Use GNU Parallel to process all BAM files in the directory with multiple threads
find "$INPUT_DIR" -name "*Aligned.sortedByCoord.out.bam" | parallel -j 8 tag_bam_file

# Merge all tagged BAM files into a single BAM file
echo "Merging all tagged BAM files into $MERGED_BAM"
samtools merge "$MERGED_BAM" "$OUTPUT_DIR"/*_tagged.bam

# Index the merged BAM file
echo "Indexing $MERGED_BAM"
samtools index "$MERGED_BAM"

echo "All done! Merged BAM file: $MERGED_BAM"

done

#nohup bash /mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/BAMfile_merge_celltag.sh &