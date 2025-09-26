
# Generate counts
# Base path where the directories are located
base_path="/mnt/raidtmp/Alejandro/Smart-seq_Athero_vanlandewijck/"

# Path to the featureCounts binary
featureCounts_path="/mnt/raidtmp/Alejandro/subread-2.0.6-Linux-x86_64/bin/featureCounts"

# Path to the annotation file
annotation_file="/mnt/raidtmp/Alejandro/GRCh38/gencode.v43.annotation.gtf"

# Iterate through directories ending with _bam
for x in SRX23809079_bam
do
    # Define the output file name
    output_file="${base_path}featureCounts/${x}_counts.txt"
    
    # Run featureCounts for all .bam files in the current directory
    $featureCounts_path -T 4 -a $annotation_file -o $output_file ${base_path}/${x}/*.bam
done

#nohup bash /mnt/raidtmp/Alejandro/Smart-seq_Athero_vanlandewijck/Featurecounts.sh &