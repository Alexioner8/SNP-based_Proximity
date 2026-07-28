#!/bin/bash

# Define the sample IDs
samples=("carotid_4_Symptomatic" "carotid_6_Asymptomatic" "carotid_9_Symptomatic" "carotid_13_Symptomatic" "carotid_14_Symptomatic" "carotid_15_Asymptomatic" "carotid_16_Asymptomatic" "carotid_17_Asymptomatic")

# Set the path
path1="/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Monopogen"
export LD_LIBRARY_PATH=/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Monopogen/apps:$LD_LIBRARY_PATH

for sample in "${samples[@]}"; do

# Define sample-specific paths
bam_path="/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/${sample}/possorted_genome_bam.bam"
cell_reads_path="/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/${sample}_cell_reads.csv"
out_path="out_${sample}"

# Create bam.somatic.lst file
echo "${sample},${bam_path}" > /mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Monopogen_out/bam.somatic.lst

# Run Monopogen preProcess
    python3.13 ${path1}/src/Monopogen.py preProcess -b /mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Monopogen_out/bam.somatic.lst \
     -o ${out_path} -a ${path1}/apps -t 6

# Run Monopogen germline
    python3.13 ${path1}/src/Monopogen.py germline \
        -a ${path1}/apps -t 6 -r ${path1}/resource/GRCh38.region.somatic_call.lst \
        -p /mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/1KG3_reference/ \
        -g /mnt/raidexttmp/Alejandro/GRCh38/GRCh38.primary_assembly.genome.fa -m 3 -s all -o ${out_path}

# Run Monopogen somatic (featureInfo)
    python3.13 ${path1}/src/Monopogen.py somatic \
        -a ${path1}/apps -t 11 -r ${path1}/resource/GRCh38.region.somatic_call.lst \
        -i ${out_path} -l ${cell_reads_path} -s featureInfo \
        -g /mnt/raidexttmp/Alejandro/GRCh38/GRCh38.primary_assembly.genome.fa
    
   

# Run Monopogen somatic (cellScan)
python3.13 ${path1}/src/Monopogen.py somatic \
            -a ${path1}/apps -r ${path1}/resource/GRCh38.region.somatic_call.lst -t 6 \
            -i ${out_path} -l ${cell_reads_path} -s cellScan \
            -g /mnt/raidexttmp/Alejandro/GRCh38/GRCh38.primary_assembly.genome.fa

# Run Monopogen somatic (LDrefinement)
python3.13 ${path1}/src/Monopogen.py somatic \
            -a ${path1}/apps -r ${path1}/resource/GRCh38.region.somatic_call.lst -t 6 \
            -i ${out_path} -l ${cell_reads_path} -s LDrefinement \
            -g /mnt/raidexttmp/Alejandro/GRCh38/GRCh38.primary_assembly.genome.fa

done





#nohup bash /mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/monopogen_process_nohup.sh &

