
######

source("/mnt/raidexttmp/Alejandro/functions.R")
setwd("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR")
library(SoupX)

#
# this searches for all input matrices
#
files=list()
for (x in c("HV1", "HV2"))
{
file<- Sys.glob(paste(paste("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/", x, sep=""), "/outs/raw_feature_bc_matrix/matrix.mtx.gz", sep=""))

files[x]= file
}




#
# here we iterate over all input matrices and load them into allfiles (gene expression) and allABs (antibody capture)
#
allfiles.raw = list()
for (file in files)
{
  samplename = sub(".*matrix","",str_split(dirname(file), "/")[[1]][6]) #estos numeros indican el nombre de la carpeta
  foldername = dirname(file)
  
  print(paste(samplename, foldername))
  
  h5file = Read10X(foldername,unique.features = TRUE)

  if (is.null(names(h5file)))
  {
      print(paste("WITHOUT AB", samplename))
    allfiles.raw[[samplename]] = h5file
  } else {
      print(paste("WITH AB", samplename))
    allfiles.raw[[samplename]] = h5file$`Gene Expression`
    allABs.raw[[samplename]] = h5file$`Antibody Capture`
  }

  print(paste(samplename, nrow(allfiles.raw[[samplename]]), "x", ncol(allfiles.raw[[samplename]]), "genes x cells"))
}

length(allfiles.raw)
names(allfiles.raw)


#
# here we create a list of seurat object. each entry corresponds to an input matrix from above
#

objlist = list()
for (x in names(allfiles.raw))
{

    matrix = allfiles.raw[[x]]
    
    # this creates a Seurat object from the count matrix. it sets the object's project to x and prepends the sample name to all cells
    # the patternlist.mouse contains patterns for mt and RP-genes
    filteredObj = makeSeuratObj(matrix, x, patternList.human)
    
    # this creates log-normalized count matrices in RNA assay
    filteredObj <- NormalizeData(filteredObj, verbose = FALSE)
    # this calculates the most (2000) variable features per data set. variable features are features which show a high variance between all cells of a sample
    filteredObj <- FindVariableFeatures(filteredObj, verbose = FALSE)
    
    objlist[[x]] = filteredObj

    print(x)
    print(filteredObj)
    
    
}

names(objlist)

objlist.raw = objlist

objlist <- lapply(X = objlist.raw, FUN = function(obj) {
  # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
  print(paste("Seurat obj project", obj@project.name))
  print(obj)
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA > 300)
  obj <- subset(obj, subset = percent.mt < 15)
  print(obj)
  
  return(obj)
})


for (name in names(objlist))
{
  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  save_plot(p, paste("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/QC", paste(name, "filtered_violins_qc", sep="_"), sep="/"), fig.width=10, fig.height=6)

  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine=F)
  p[[1]] = p[[1]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[2]] = p[[2]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[3]] = p[[3]] + scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,5))
  p = combine_plot_grid_list(plotlist=p, ncol=3)
  save_plot(p, paste("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/QC", paste(name, "filtered_violins_detail_qc", sep="_"), sep="/"), fig.width=18, fig.height=6)
  
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/QC", paste(name, "filtered_scatter_ncount_mt", sep="_"), sep="/"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/QC", paste(name, "filtered_scatter_ncount_rp", sep="_"), sep="/"), fig.width=10, fig.height=6)
}




####Integration

prepareFinalList = function(finalList)
{



print("cells per experiment")
print(mapply(sum, lapply(finalList, function(x) {dim(x)[2]})))
print("total cells")
print(sum(mapply(sum, lapply(finalList, function(x) {dim(x)[2]}))))


objlist = list()
for (objname in names(finalList))
{

    x = finalList[[objname]]  
    print(objname)
    print(paste("Seurat obj project", x@project.name))
    
    Project(x) = objname
    print(paste("Seurat obj project", x@project.name))

    DefaultAssay(x) = "RNA"
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)

    x$library = objname

    objlist[[objname]] = x
}

features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 2000)

objlist <- lapply(X = objlist, FUN = function(x) {

    print(paste("Seurat obj project", x@project.name))
    print(x)



    features_use <- intersect(features, rownames(x))
    x <- ScaleData(x, features = features_use, verbose = FALSE, assay="RNA", vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA"))
    x <- RunPCA(x, verbose = FALSE, reduction.name="pca", assay="RNA")
    #x <- suppressWarnings(SCTransform(x, verbose = FALSE,vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score")))
    
    plot1 <- ElbowPlot(x)
    save_plot(plot1, paste("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/QC", paste(x@project.name, "elbowplot_dimensionality", sep="_"), sep="/"), fig.width=10, fig.height=6)


    x$project = x@project.name

    return(x)
})

return(objlist)

}



analyseFinalList = function(objlist, intname) 
{

dir.create(intname, recursive = TRUE)


#
# integrate based on RNA/GEX assay
#
objSamples = objlist
print(objSamples)

objSamples = lapply(objSamples, function(x) {
  DefaultAssay(x) <- 'RNA'
  x <- RunPCA(x, verbose = FALSE, reduction.name="pca",  assay="RNA")
#  DefaultAssay(x) <- 'SCT'
#  print(x@reductions$pca)
  return(x)
})
print("GEX integration features")
print(objSamples)

features_gex <- SelectIntegrationFeatures(object.list = objSamples, nfeatures = 2000)
objlist.anchors <- FindIntegrationAnchors(object.list = objSamples,  reduction = "cca", dims = 1:30, anchor.features = features_gex) 
obj.list.integrated <- IntegrateData(new.assay.name = "integratedgex", anchorset = objlist.anchors, dims = 1:30, verbose=T, k.weight = 30) 
print("GEX integration done")

#
# integrated GEX viz
#
obj.list.integrated = ScaleData(obj.list.integrated, assay="integratedgex")
obj.list.integrated <- RunPCA(obj.list.integrated, reduction.name="igpca", assay="integratedgex")
obj.list.integrated <- RunUMAP(obj.list.integrated, reduction = "igpca", dims = 1:30, reduction.key = "UMAPig_",)
p=DimPlot(obj.list.integrated, reduction="umap", shuffle = T, seed = 1, group.by= "orig.ident")
save_plot(p, paste(intname, "wnn_ig_dimplot", sep="/"), 12, 8)

p=DimPlot(obj.list.integrated, reduction="igpca", group.by= "orig.ident")
save_plot(p, paste(intname, "wnn_pca_ig_dimplot", sep="/"), 12, 8)

p=DimPlot(obj.list.integrated, shuffle = T, seed = 1, group.by= "orig.ident")
save_plot(p, "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/orig.ident_dimplot", 12, 8)





return(obj.list.integrated)

}


finalList_sample = prepareFinalList(objlist)
integratedList = analyseFinalList(finalList_sample, "HV_int")


source("/mnt/raidexttmp/Alejandro/functions.R")
setwd("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR")
#saveRDS(integratedList, file = "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/HV.rds")
#integratedList=readRDS("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/HV.rds")




##############

obj.in <- FindVariableFeatures(integratedList, selection.method = "vst")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj.in), 30)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj.in) 
p= LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(p, "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/HV_feat", 10, 10)


#cluster identification
DefaultAssay(integratedList) = "integratedgex"
obj.in <- FindNeighbors(integratedList, reduction="igpca", dims = 1:30)
DefaultAssay(obj.in) <- "integratedgex"
obj.in <- FindClusters(obj.in, resolution = 0.2, algorithm = 4)

p=DimPlot(obj.in, pt.size = 2, label=T, reduction = "umap")
save_plot(p, "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/Dimplot_umap_ident", fig.width=10, fig.height=8)


############
# Define the sample IDs
samples <- c("HV1", "HV2")

# Loop over each sample
for (x in samples) {
  

  sample_data <- obj.in
  
  # Extract barcode and read count
  cell_barcodes <- colnames(sample_data)
  num_reads <- sample_data$nCount_RNA
  
  # Remove prefix and "-1" from barcodes
  cleaned_barcodes <- gsub(paste0("^", x, "_"), "", cell_barcodes)
  
  # Create output data frame
  output_data <- data.frame(cell = cleaned_barcodes, id = num_reads)
  
  # Define output file path
  output_file <- paste0("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/", x, "_cell_reads.csv")
  
  # Write to CSV
  write.csv(output_data, file = output_file, row.names = FALSE, quote = FALSE)
}


################


# List of CSVs and their sample prefixes
tcr_files <- c(
  HV1 = "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/TCR_seq/HV1_vdj.csv",
  HV2 = "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/TCR_seq/HV2_vdj.csv"
)

# Initialize empty metadata
all_meta <- data.frame(
  cell = colnames(obj.in),
  clonotype = NA,
  stringsAsFactors = FALSE
)

# Loop through samples
for (sample_name in names(tcr_files)) {

  # Load CSV
  df <- read.csv(tcr_files[sample_name], stringsAsFactors = FALSE)

  # Summarize clonotype per barcode
  df_meta <- df %>%
    group_by(barcode) %>%
    summarise(raw_clonotype_id = first(raw_clonotype_id), .groups = "drop")

  # Select Seurat barcodes from this sample only
  sample_cells <- colnames(obj.in)[grepl(paste0("^", sample_name, "_"), colnames(obj.in))]
  sample_cells_clean <- gsub(paste0("^", sample_name, "_"), "", sample_cells)

  # Build temporary metadata
  temp_meta <- data.frame(
    cell = sample_cells,
    clonotype = df_meta$raw_clonotype_id[match(sample_cells_clean, df_meta$barcode)],
    stringsAsFactors = FALSE
  )

  # Prepend sample name to make clonotypes unique
  temp_meta$clonotype <- ifelse(!is.na(temp_meta$clonotype),
                                paste0(sample_name, " ", temp_meta$clonotype),
                                NA)

  # Filter clonotypes with >=3 cells
  clonotype_counts <- table(temp_meta$clonotype)
  valid_clonotypes <- names(clonotype_counts[clonotype_counts > 3])
  temp_meta$clonotype_filtered <- ifelse(temp_meta$clonotype %in% valid_clonotypes,
                                         temp_meta$clonotype,
                                         NA)

  # Update main metadata only for these sample cells
  all_meta$clonotype[match(temp_meta$cell, all_meta$cell)] <- temp_meta$clonotype_filtered
}

# Add the filtered clonotype metadata to Seurat object
obj.in <- AddMetaData(obj.in, metadata = all_meta$clonotype, col.name = "clonotype")

obj= obj.in[,obj.in$clonotype!="NA"]
# Plot TSNE colored by filtered clonotype
p <- DimPlot(obj, pt.size = 1, reduction = "umap", group.by = "clonotype", 
        label = FALSE, shuffle = TRUE) +
  labs(
    title = "HV T-cell clonotypes (>3 cells)",
    x = "UMAP 1",
    y = "UMAP 2"
  )

save_plot(p, "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/Dimplot_clonotypes_filtered3cells", fig.width = 8, fig.height = 6)



###########Monopogen SNPs

# Set the directory containing your CSV files
setwd("/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/Monopogen_run/out_HV1/somatic")  # Change this to your actual path

# List all CSV files that match the pattern "chr*.putativeSNVs.csv"
file_list <- list.files(pattern = "chr[0-9]+\\.putativeSNVs\\.csv")

# Read and merge all files
library(dplyr)  # For data manipulation

merged_data <- file_list %>%
  lapply(read.csv) %>%
  bind_rows()

# Save the merged file
write.csv(merged_data, "merged.putativeSNVs.csv", row.names = FALSE)

rds_list <- list.files(pattern = "chr[0-9]+\\.SNV_mat\\.RDS")

merged_rds <- rds_list %>%
  lapply(readRDS)%>%
  bind_rows()

colnames(merged_rds) <- c(colnames(merged_rds[,1:18]), paste0("HV1_", colnames(merged_rds[,19:ncol(merged_rds)])))


saveRDS(merged_rds, file= "merged.SNV_mat.RDS")

#################

# Define sample names and paths
samples <- c(HV1 = "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/Monopogen_run/out_HV1/somatic/",
             HV2 = "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/Monopogen_run/out_HV2/somatic/")

# Initialize empty matrix list
all_mat <- list()

for (sample_name in names(samples)) {
  
  cat("Processing sample:", sample_name, "\n")
  
  # Load meta and matrix files
  meta_file <- file.path(samples[sample_name], "merged.putativeSNVs.csv")
  mat_file  <- file.path(samples[sample_name], "merged.SNV_mat.RDS")
  
  meta <- read.csv(meta_file)
  mat  <- readRDS(mat_file)
  
  # Filter variants based on depth and SVM score
  meta_filter <- meta[meta$Depth_ref > 5 & meta$Depth_alt > 5, ]
  meta_filter <- meta_filter[meta_filter$SVM_pos_score > 0.1, ]
  snv_pass <- paste0(meta_filter$chr, ":", meta_filter$pos, ":", meta_filter$Ref_allele, ":", meta_filter$Alt_allele)
  
  mat <- mat[snv_pass, , drop = FALSE]
  
  # Remove columns with any NA
  mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]
  
  # Only keep cells that overlap with the Seurat object
  overlap <- intersect(colnames(mat), colnames(obj.in))
  mat <- mat[, overlap, drop = FALSE]
  
  # Convert data.frame to matrix for numeric manipulation
  mat <- as.matrix(mat)
  
  # Extract numeric values
  split_values <- do.call(rbind, strsplit(mat, "/", fixed = TRUE))
  left_values  <- as.numeric(split_values[, 1])
  right_values <- as.numeric(split_values[, 2])
  
  # Adapt SNV matrix
  mat[left_values == 0] <- 1
  mat[right_values == 0] <- -1
  mat[right_values == 0 & left_values == 0] <- 0
  mat[right_values > 0] <- 1
  mat[left_values > (right_values * 2)] <- -1
  
  # Convert back to data.frame (optional)
  mat <- as.data.frame(mat)
  
  # Store in list
  all_mat[[sample_name]] <- mat
}



# Create Seurat assay
Geno1 <- CreateAssayObject(counts = all_mat[["HV1"]])
obj.in[["HV1SNP"]] <- Geno1
Geno2 <- CreateAssayObject(counts = all_mat[["HV2"]])
obj.in[["HV2SNP"]] <- Geno2

# Keep only cells with clonotype info
obj <- obj.in[, !is.na(obj.in$clonotype)]


#########################################################################
library(Matrix)
library(data.table)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)

# Get patient IDs
patients <- unique(obj$orig.ident)


# Store results
all_alluvial_data <- list()

for (patient in patients) {
  cat("Processing patient:", patient, "\n")

  patient_cells <- colnames(obj)[obj$orig.ident == patient]
  if (length(patient_cells) < 2) next
  
  # Extract SNP assay matrix
  snp_matrix <- obj[[paste0(patient, "SNP")]]@counts


  snp_presence <- snp_matrix

  cell_types <- obj$clonotype



  # Compute SNP sharing matrix
  shared_snp_mat <- as.matrix(crossprod(snp_presence))
  diag(shared_snp_mat) <- 0

  # Max shared SNPs per cell
  max_shared_per_cell <- apply(shared_snp_mat, 1, max)

  # Compute 99th percentile threshold
  threshold_99 <- quantile(max_shared_per_cell, 0.99)

  # Select cells above threshold
  filtered_cells <- names(max_shared_per_cell[max_shared_per_cell >= threshold_99])
  if (length(filtered_cells) < 2) next

  filtered_mat <- shared_snp_mat[filtered_cells, filtered_cells]

  # Nearest neighbors
  closest_df <- data.table(
    source_cell = rownames(filtered_mat),
    target_cell = apply(filtered_mat, 1, function(x) names(x)[which.max(x)])
  )

  closest_df[, source_cluster := cell_types[source_cell]]
  closest_df[, target_cluster := cell_types[target_cell]]

  # Transition summary
  alluvial_data <- closest_df[, .N, by = .(source_cluster, target_cluster)]
  alluvial_data[, patient := patient]
  all_alluvial_data[[patient]] <- alluvial_data
}


# Combine and plot as before
combined_alluvial <- rbindlist(all_alluvial_data)
combined_alluvial

# Compute global counts
TP <- sum(combined_alluvial$N[combined_alluvial$source_cluster == combined_alluvial$target_cluster])
FP <- sum(combined_alluvial$N[combined_alluvial$source_cluster != combined_alluvial$target_cluster])

FDR <- FP / (TP + FP)
precision <- TP / (TP + FP)


#
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
 
color_list <- ggplotColours(n=5)

# Summarize total transitions
label_data <- combined_alluvial[, .(total_N = sum(N)), by = .(source_cluster, target_cluster)]
label_data[, label := total_N * 2]

# Compute max y to place labels nicely above the flows
max_y <- max(label_data$total_N)

p <- ggplot(label_data,
            aes(axis1 = source_cluster, axis2 = target_cluster, y = total_N)) +
  geom_alluvium(aes(fill = source_cluster), width = 1/12, alpha = 0.9) +
  geom_stratum(width = 1/12, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 4, fontface = "bold", color = "black") +
  # Label only unique transitions with total N*2
  geom_text(
    data = label_data,
    aes(label = label, group = interaction(source_cluster, target_cluster)),
    stat = "flow", hjust = -5,
    size = 4, fontface = "bold", color = "black"
  ) +
  scale_fill_manual(values = color_list) +
  labs(
    title = "SNV-based NN-cells TCR-clonotype control",
    y = "Number of Cells", x = NULL
  ) +
  # Add annotations (adjust x, y positions as needed)
  annotate("text", x = 2.5, y = max_y * 1.1, 
           label = "FDR 95 percentile = 0.028",
           hjust = 1, size = 5, fontface = "italic") +
  annotate("text", x = 2.5, y = max_y * 1.05, 
           label = "FDR 99 percentile = 0",
           hjust = 1, size = 5, fontface = "italic") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )


save_plot(p, "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/Alluvial_TCR_Filtered", 8, 8)




########################
library(Matrix)
library(data.table)
library(ggplot2)

# Define thresholds (99% to 90%)
percentiles <- seq(0.99, 0.80, by = -0.02)

# Store results
results <- list()

for (q in percentiles) {
  cat("Processing percentile:", q, "\n")

  TP_all <- 0
  FP_all <- 0
  
  for (patient in patients) {
    patient_cells <- colnames(obj)[obj$orig.ident == patient]
    if (length(patient_cells) < 2) next
    
    # Extract SNP assay matrix
    snp_matrix <- obj[[paste0(patient, "SNP")]]@counts
    snp_presence <- snp_matrix
    cell_types <- obj$clonotype
    
    # Compute SNP sharing
    shared_snp_mat <- as.matrix(crossprod(snp_presence))
    diag(shared_snp_mat) <- 0
    
    max_shared_per_cell <- apply(shared_snp_mat, 1, max)
    
    # Threshold for this quantile
    threshold <- quantile(max_shared_per_cell, q)
    
    filtered_cells <- names(max_shared_per_cell[max_shared_per_cell >= threshold])
    if (length(filtered_cells) < 2) next
    
    filtered_mat <- shared_snp_mat[filtered_cells, filtered_cells]
    
    # Nearest neighbors
    closest_df <- data.table(
      source_cell = rownames(filtered_mat),
      target_cell = apply(filtered_mat, 1, function(x) names(x)[which.max(x)])
    )
    
    closest_df[, source_cluster := cell_types[source_cell]]
    closest_df[, target_cluster := cell_types[target_cell]]
    
    # Count true/false positives
    TP <- sum(closest_df$source_cluster == closest_df$target_cluster)
    FP <- sum(closest_df$source_cluster != closest_df$target_cluster)
    
    TP_all <- TP_all + TP
    FP_all <- FP_all + FP
  }
  
  precision <- ifelse((TP_all + FP_all) > 0, TP_all / (TP_all + FP_all), NA)
  FDR <- ifelse((TP_all + FP_all) > 0, FP_all / (TP_all + FP_all), NA)
  
  results[[as.character(q)]] <- data.table(
    percentile = q,
    precision = precision,
    FDR = FDR
  )
}

# Combine results
results_dt <- rbindlist(results)

# Convert to percentage for plotting
results_dt[, precision_pct := precision * 100]

library(ggplot2)
library(scales)

# Improved barplot
p <- ggplot(results_dt, aes(x = factor(percentile), y = precision_pct, fill = percentile)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("FDR = ", round(FDR, 3))),
            color = "white", fontface = "bold", size = 4, angle = 90,
            position = position_stack(vjust = 0.5)) +
  scale_fill_gradient(low = "#8dbde1ff", high = "#042241ff", guide = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::percent_format(scale = 1)) +
  labs(
    title = "SNP-based NN-cell Assignment vs TCR Clonotype",
    x = "Percentile threshold",
    y = "Properly assigned T cells (%)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_plot(p, "/mnt/raidexttmp/Alejandro/Benchmark_MonopogenSNP_and_TCR/Barplot_TCR_Filtered", 8, 7)

