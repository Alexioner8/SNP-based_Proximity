######
source("/mnt/raidbio/tmp/Alejandro/functions.R")
#saveRDS(obj.in, file = "/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
obj.in=readRDS("/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")


# Define the sample IDs
samples <- c("SRX23809074", "SRX23809075", "SRX23809076", "SRX23809077", "SRX23809078", 
             "SRX23809080", "SRX23809082", "SRX23809083", "SRX23809085", "SRX23809086")

# Loop over each sample
for (x in samples) {
  
  # Identify the correct prefix ("Asymptomatic_" or "Symptomatic_")
  prefix <- ifelse(any(grepl(paste0("Asymptomatic_", x), obj.in$library)), "Asymptomatic_", "Symptomatic_")
  
  # Subset the data based on the identified prefix
  sample_name <- paste0(prefix, x)
  sample_data <- obj.in[, obj.in$library == sample_name]
  
  # Extract barcode and read count
  cell_barcodes <- colnames(sample_data)
  num_reads <- sample_data$nCount_RNA
  
  # Remove prefix and "-1" from barcodes
  cleaned_barcodes <- gsub(paste0("^", sample_name, "_"), "", cell_barcodes)
  
  # Create output data frame
  output_data <- data.frame(cell = cleaned_barcodes, id = num_reads)
  
  # Define output file path
  output_file <- paste0("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/", x, "_cell_reads.csv")
  
  # Write to CSV
  write.csv(output_data, file = output_file, row.names = FALSE, quote = FALSE)
}



###########################
##########################################
#!/bin/bash

# Define the sample IDs
samples=("SRX23809085" "SRX23809086" "SRX23809075" "SRX23809076" "SRX23809077" "SRX23809078" "SRX23809082" "SRX23809083")

# Set the path
path1="/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path1}/apps

# Loop over each sample
for sample in "${samples[@]}"; do

    # Determine prefix (Asymptomatic or Symptomatic)
    if [[ -f "/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${sample}_bam_tagged/merged_tagged.bam" ]]; then
        prefix="Asymptomatic"
    else
        prefix="Symptomatic"
    fi

    # Define sample-specific paths
    bam_path="/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/${sample}_bam_tagged/merged_tagged.bam"
    cell_reads_path="/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/${sample}_cell_reads.csv"
    out_path="out_${sample}"

    # Create bam.somatic.lst file
    echo "${prefix}_${sample},${bam_path}" > /mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/bam.somatic.lst

    # Run Monopogen preProcess
    python3.12 ${path1}/src/Monopogen.py preProcess -b bam.somatic.lst -o ${out_path} -a ${path1}/apps -t 6

    # Run Monopogen germline
    python3.12 ${path1}/src/Monopogen.py germline \
        -a ${path1}/apps -t 6 -r ${path1}/resource/GRCh38.region.somatic_call.lst \
        -p /mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/1KG3_reference/ \
        -g /mnt/raidexttmp/Alejandro/GRCh38/GRCh38.primary_assembly.genome.fa -m 3 -s all -o ${out_path}

    # Run Monopogen somatic (featureInfo)
    python3.12 ${path1}/src/Monopogen.py somatic \
        -a ${path1}/apps -t 11 -r ${path1}/resource/GRCh38.region.somatic_call.lst \
        -i ${out_path} -l ${cell_reads_path} -s featureInfo \
        -g /mnt/raidexttmp/Alejandro/GRCh38/GRCh38.primary_assembly.genome.fa

    # Loop through chromosomes
    for chr in {1..22}; do
        echo "chr${chr}" >> ${path1}/resource/GRCh38.region.somatic_cellScan.lst
        
        # Run Monopogen somatic (cellScan)
        python3.12 ${path1}/src/Monopogen.py somatic \
            -a ${path1}/apps -r ${path1}/resource/GRCh38.region.somatic_cellScan.lst -t 6 \
            -i ${out_path} -l ${cell_reads_path} -s cellScan \
            -g /mnt/raidexttmp/Alejandro/GRCh38/GRCh38.primary_assembly.genome.fa

        # Run Monopogen somatic (LDrefinement)
        python3.12 ${path1}/src/Monopogen.py somatic \
            -a ${path1}/apps -r ${path1}/resource/GRCh38.region.somatic_cellScan.lst -t 6 \
            -i ${out_path} -l ${cell_reads_path} -s LDrefinement \
            -g /mnt/raidexttmp/Alejandro/GRCh38/GRCh38.primary_assembly.genome.fa

        # Clean up temporary chromosome file
        rm ${path1}/resource/GRCh38.region.somatic_cellScan.lst
    done

    # Remove bam.somatic.lst after processing each sample
    rm /mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/bam.somatic.lst

done
