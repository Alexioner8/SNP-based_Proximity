

source("/mnt/raidbio/tmp/Alejandro/functions.R")
#saveRDS(obj.in, file = "/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
obj.in=readRDS("/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
new.cluster.ids <- c("SMCs", "MPs","ECs", "T-cells","PCs")
names(new.cluster.ids) <- levels(obj.in)
obj.in <- RenameIdents(obj.in, new.cluster.ids)
obj.in$CellType <- Idents(obj.in)



# Loop through each unique library in obj.in
exclude_libraries <- c("Asymptomatic_SRX23809081", "Symptomatic_SRX23809084")
filtered_libraries <- setdiff(unique(obj.in$orig.ident), exclude_libraries)

obj.list = list()

for (x in filtered_libraries) {
  # Subset the obj.in Seurat object based on the current library
  sample <- obj.in[, obj.in$orig.ident == x]
  
  # Get the directory name dynamically
  library_name <- unique(sample$orig.ident)  # Extract the library name
  file_path <- paste0("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/out_", gsub(".*_", "", library_name), "/somatic/merged.SNV_mat.RDS")
  
  # Read the corresponding merged SNV matrix
  mat <- readRDS(file = file_path)
  
  meta<- read.csv(file=paste0("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/out_", gsub(".*_", "", library_name), "/somatic/merged.putativeSNVs.csv"))

  meta_filter <-meta[meta$Depth_ref>5 & meta$Depth_alt>5,]
  meta_filter <- meta_filter[meta_filter$SVM_pos_score>0.1,]
  snv_pass <- paste0(meta_filter$chr,":", meta_filter$pos, ":", meta_filter$Ref_allele,":", meta_filter$Alt_allele)
  mat <- mat[snv_pass,]

  
  # Find overlap between the sample and the matrix columns
  overlap <- intersect(colnames(mat), colnames(sample))
  
  if (length(overlap) > 0) {
    # Subset the matrix and sample based on the overlapping columns
    mat <- mat[, overlap, drop = FALSE]
    sample <- subset(sample, cells = overlap)
    
    # Convert data.frame to matrix for easier manipulation
    mat <- as.matrix(mat)
    
    # Extract numbers using regular expressions
    split_values <- do.call(rbind, strsplit(mat, "/", fixed = TRUE))
    
    # Convert extracted values to numeric
    left_values <- as.numeric(split_values[, 1])
    right_values <- as.numeric(split_values[, 2])
    
    # Adapt SNV matrix
    mat[left_values == 0] <- 1
    mat[right_values == 0] <- (-1)
    mat[right_values == 0 & left_values == 0] <- 0
    mat[right_values > 0] <- 1
    
    # Convert back to data.frame if needed
    mat <- as.data.frame(mat)
    
    # Identify all cell types present in this sample
    sample_celltypes <- unique(sample$CellType)

    # Function to get SNPs detected in >80% of cells per cluster
    get_high_conf_snps <- function(cells, threshold = 0.8) {
    submat <- mat[, cells, drop = FALSE]
    detected_fraction <- rowMeans(submat != 0)
    rownames(submat)[detected_fraction >= threshold]
    }

    # For each cluster, get high-confidence SNPs
    cluster_snps <- lapply(sample_celltypes, function(ct) {
    cells_in_ct <- colnames(sample)[sample$CellType == ct]
    get_high_conf_snps(cells_in_ct, threshold = 0.8)
    })

    # Keep only SNPs shared across all clusters
    common_snps <- Reduce(intersect, cluster_snps)

    # Subset mat to keep only common SNPs
    mat_filtered <- mat[common_snps, , drop = FALSE]

    # Add the filtered SNP matrix to the list
    obj.list[[x]] <- mat_filtered
  }
}

library(Matrix)

# Convert obj.list elements to sparse matrices
sparse_list <- lapply(obj.list, function(df) {
  # Convert to numeric matrix first
  mat <- as.matrix(df)
  mat <- apply(mat, 2, as.numeric)  # Ensure values are numeric
  # Convert to sparse matrix
  sparse_mat <- Matrix(mat, sparse = TRUE)
  rownames(sparse_mat) <- rownames(df)  # Copy row names
  return(sparse_mat)
})

# Remove the original obj.list to free up memory
rm(obj.list)
gc()  # Run garbage collection


# Step 1: Collect all unique SNP positions from row names
all_snp_positions <- unique(unlist(lapply(sparse_list, function(mat) rownames(mat))))

# Step 2: Create an empty sparse matrix to store merged data
merged_matrix <- Matrix(0, nrow = length(all_snp_positions), 
                        ncol = sum(sapply(sparse_list, ncol)), sparse = TRUE)

# Assign row names
rownames(merged_matrix) <- all_snp_positions

# Store column names
col_names <- c()
col_start <- 1

# Step 3: Align each matrix and fill merged matrix
for (mat_name in names(sparse_list)) {
  mat <- sparse_list[[mat_name]]
  mat_rownames <- rownames(mat)
  
  # Match SNP positions to unified set
  match_idx <- match(mat_rownames, all_snp_positions)
  
  # Determine column indices for placement
  col_end <- col_start + ncol(mat) - 1
  merged_matrix[match_idx, col_start:col_end] <- mat
  col_start <- col_end + 1
  
  # Store column names
  col_names <- c(col_names, colnames(mat))
}

# Assign column names
colnames(merged_matrix) <- col_names


# Save merged sparse matrix
saveRDS(merged_matrix, file = "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Allpatients.SNV_sparse_mat_merged.RDS")





# Convert sparse matrix to a dense format
source("/mnt/raidbio/tmp/Alejandro/functions.R")
library(Matrix)
#saveRDS(obj.in, file = "/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
obj.in=readRDS("/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
merged_matrix=readRDS("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Allpatients.SNV_sparse_mat_merged.RDS")
new.cluster.ids <- c("SMCs", "MPs","ECs", "T-cells","PCs")
names(new.cluster.ids) <- levels(obj.in)
obj.in <- RenameIdents(obj.in, new.cluster.ids)
obj.in$CellType <- Idents(obj.in)



###################alluvial plot showing connections from source cell type to high confidence nearest-neighbor cell type
library(Matrix)
library(data.table)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)

source("/mnt/raidbio/tmp/Alejandro/functions.R")

# Load data
obj.in <- readRDS("/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
merged_matrix <- readRDS("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Allpatients.SNV_sparse_mat_merged.RDS")

# Rename clusters
new.cluster.ids <- c("SMCs", "MPs","ECs", "T-cells","PCs")
names(new.cluster.ids) <- levels(obj.in)
obj.in <- RenameIdents(obj.in, new.cluster.ids)
obj.in$CellType <- Idents(obj.in)


# Keep only overlapping cells
overlap <- intersect(colnames(merged_matrix), colnames(obj.in))
obj <- subset(obj.in, cells = overlap)

# Add SNP matrix to Seurat
Geno <- CreateAssayObject(counts = merged_matrix)
obj[["Geno"]] <- Geno

# Feature vs SNPs plot
mut_snp_per_cell <- colSums(merged_matrix == 1)
obj[["totalSNPs"]] <- mut_snp_per_cell

#Extract data
df_plot <- data.frame(
  Counts = obj$nCount_RNA,
  SNPs = obj[["totalSNPs"]][, 1]
)

# Scatter plot with linear model and correlation coefficient
# Perform Pearson correlation test
cor_test <- cor.test(df_plot$Counts, df_plot$SNPs, method = "pearson")

# Extract r and p
cor_val <- cor_test$estimate     # the Pearson r
p_val <- cor_test$p.value        # the p-value

# Format correlation label
cor_label <- paste0(
  "r = ", round(cor_val, 3), 
  ", p", format.pval(p_val, digits = 3, scientific = TRUE)
)


# Get line fit values to estimate annotation position
fit <- lm(SNPs ~ Counts, data = df_plot)
x_pos <- quantile(df_plot$Counts, 0.90)
y_pos <- predict(fit, newdata = data.frame(Counts = x_pos))*3.5

# Plot
p_scatter <- ggplot(df_plot, aes(x = Counts, y = SNPs)) +
  geom_point(alpha = 0.5, size = 1.8, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#F1948A", size = 1) +
  annotate("text", x = x_pos, y = y_pos,
           label = cor_label,
           hjust = 0, vjust = -1.5,
           size = 6, color = "#E74C3C", fontface = "bold") +
  theme_minimal(base_size = 16) +
  labs(
    x = "RNA Counts per cell",
    y = "Total SNV Counts per cell"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Save plot
ggsave(
  filename = "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/counts_vs_snp_scatter_lm.png",
  plot = p_scatter,
  width = 7,
  height = 6
)

#Extract data
df_plot <- data.frame(
  Counts = obj$nFeature_RNA,
  SNPs = obj[["totalSNPs"]][, 1]
)

# Scatter plot with linear model and correlation coefficient
# Perform Pearson correlation test
cor_test <- cor.test(df_plot$Counts, df_plot$SNPs, method = "pearson")

# Extract r and p
cor_val <- cor_test$estimate     # the Pearson r
p_val <- cor_test$p.value        # the p-value

# Format correlation label
cor_label <- paste0(
  "r = ", round(cor_val, 3), 
  ", p", format.pval(p_val, digits = 3, scientific = TRUE)
)


# Get line fit values to estimate annotation position
fit <- lm(SNPs ~ Counts, data = df_plot)
x_pos <- quantile(df_plot$Counts, 0.90)
y_pos <- predict(fit, newdata = data.frame(Counts = x_pos))*3

# Plot
p_scatter <- ggplot(df_plot, aes(x = Counts, y = SNPs)) +
  geom_point(alpha = 0.5, size = 1.8, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#F1948A", size = 1) +
  annotate("text", x = x_pos, y = y_pos,
           label = cor_label,
           hjust = 0, vjust = -1.5,
           size = 6, color = "#E74C3C", fontface = "bold") +
  theme_minimal(base_size = 16) +
  labs(
    x = "Feature Counts per cell",
    y = "Total SNV Counts per cell"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Save plot
ggsave(
  filename = "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/features_vs_snp_scatter_lm.png",
  plot = p_scatter,
  width = 7,
  height = 6
)

# Get patient IDs
patients <- unique(obj$orig.ident)[unique(obj$orig.ident)!="Asymptomatic_SRX23809079"]


get_high_conf_snps <- function(cells, threshold = 0.8) {
  mat <- merged_matrix[, cells, drop = FALSE]
  detected_fraction <- rowMeans(mat != 0)
  rownames(mat)[detected_fraction >= threshold]
}


# Store results
# Store results
all_alluvial_data <- list()

for (patient in patients) {
  cat("Processing patient:", patient, "\n")

  patient_cells <- colnames(obj)[obj$orig.ident == patient]
  if (length(patient_cells) < 2) next

  common = get_high_conf_snps(patient_cells)
  snp_presence <- merged_matrix[common, patient_cells] == 1
  cell_types <- obj$CellType[patient_cells]
  names(cell_types) <- patient_cells

  # Compute SNP sharing matrix
  shared_snp_mat <- as.matrix(crossprod(snp_presence))
  diag(shared_snp_mat) <- 0

  # Compute z-scores of max shared SNPs
  max_shared_per_cell <- apply(shared_snp_mat, 1, max)
  mean_shared <- mean(max_shared_per_cell)
  sd_shared <- sd(max_shared_per_cell)
  z_scores <- (max_shared_per_cell - mean_shared) / sd_shared

  # Select cells with z-score > 2
  filtered_cells <- names(z_scores[z_scores > 2])
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

#
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
 
color_list <- ggplotColours(n=5)
# Plot
p <- ggplot(combined_alluvial,
            aes(axis1 = source_cluster, axis2 = target_cluster, y = N)) +
  geom_alluvium(aes(fill = source_cluster), width = 1/12, alpha = 0.9) +
  geom_stratum(width = 1/12, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 4, fontface = "bold", color = "black") +
  scale_fill_manual(values = color_list) +
  labs(
    title = "High-confidence SNV-based NN-cell",
    y = "Number of Cells", x = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

save_plot(p, "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Alluvial_80common_PatientWise_Filtered", 8, 10)



# Store patient-wise SNP pair distributions
all_shared_snp_values <- list()

for (patient in patients) {
  cat("Processing patient:", patient, "\n")

  patient_cells <- colnames(obj)[obj$orig.ident == patient]
  if (length(patient_cells) < 2) next

  common = get_high_conf_snps(patient_cells)
  snp_presence <- merged_matrix[common, patient_cells] == 1
  cell_types <- obj$CellType[patient_cells]
  names(cell_types) <- patient_cells

  # Compute SNP sharing matrix
  shared_snp_mat <- as.matrix(crossprod(snp_presence))
  diag(shared_snp_mat) <- 0

  # Store upper triangle values (i.e., unique cell pairs)
  shared_values <- shared_snp_mat[upper.tri(shared_snp_mat)]
  all_shared_snp_values[[patient]] <- data.table(
    shared_snp = shared_values,
    patient = patient
  )

  # (Optional: continue with max SNP filtering / alluvial logic...)
}

# Combine all into one data.table
shared_snp_dt <- rbindlist(all_shared_snp_values)

# Compute per-patient mean and sd
stats_dt <- shared_snp_dt[, .(
  mean_shared = mean(shared_snp),
  sd_shared = sd(shared_snp)
), by = patient]

# Merge with plotting data
shared_snp_dt <- merge(shared_snp_dt, stats_dt, by = "patient")

# Plot with mean + SD line
dist_plot <- ggplot(shared_snp_dt, aes(x = shared_snp)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", boundary = 0) +
  facet_wrap(~ patient, scales = "free_y") +
  geom_vline(aes(xintercept = mean_shared + sd_shared), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Shared SNPs per Cell Pair (Mean + SD threshold per patient)",
       x = "Number of Shared SNPs", y = "Count of Cell Pairs") +
  theme_minimal(base_size = 14)

# Save it
ggsave(
  filename = "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/SNVpair_distribution_by_patient_facet_mean_sd.png",
  plot = dist_plot, width = 12, height = 6, dpi = 300
)



################### Annotate all NN-cells for heatmap

# Store results
all_alluvial_data <- list()
all_celllevel_data <- list()

for (patient in patients) {
  cat("Processing patient:", patient, "\n")

  patient_cells <- colnames(obj)[obj$orig.ident == patient]
  if (length(patient_cells) < 2) next

  common <- get_high_conf_snps(patient_cells)
  snp_presence <- merged_matrix[common, patient_cells] == 1
  cell_types <- obj$CellType[patient_cells]
  names(cell_types) <- patient_cells

  # Compute SNP sharing matrix
  shared_snp_mat <- as.matrix(crossprod(snp_presence))
  diag(shared_snp_mat) <- 0

  # Max shared SNPs per cell
  max_shared_per_cell <- apply(shared_snp_mat, 1, max)

  # Compute 99th percentile threshold
  threshold_99 <- quantile(max_shared_per_cell, 0.7)

  # Select cells above threshold
  filtered_cells <- names(max_shared_per_cell[max_shared_per_cell > threshold_99])
  if (length(filtered_cells) < 2) next

  filtered_mat <- shared_snp_mat[filtered_cells, filtered_cells]


  closest_df <- data.table(
    source_cell = rownames(filtered_mat),
    target_cell = apply(filtered_mat, 1, function(x) names(x)[which.max(x)])
  )

  closest_df[, source_cluster := cell_types[source_cell]]
  closest_df[, target_cluster := cell_types[target_cell]]

  # Remove intra-cluster
  #closest_df <- closest_df[source_cluster != target_cluster]
  #if (nrow(closest_df) == 0) next

  # Summary for alluvial
  alluvial_data <- closest_df[, .N, by = .(source_cluster, target_cluster)]
  alluvial_data[, patient := patient]
  all_alluvial_data[[patient]] <- alluvial_data

  # Store barcode-level NN data
  closest_df[, patient := patient]
  all_celllevel_data[[patient]] <- closest_df
}


cellint_data <- rbindlist(all_celllevel_data)

library(data.table)

# Ensure data.table
setDT(cellint_data)

# STEP 1: Intra-cluster transitions only
intra <- cellint_data[source_cluster == target_cluster]
intra[, label := paste0("NN-", source_cluster, "_", source_cluster)]

intra_labels <- rbind(
  intra[, .(cell = source_cell, label)],
  intra[, .(cell = target_cell, label)]
)

# Keep one label per cell
intra_labels_unique <- intra_labels[, .N, by = .(cell, label)][order(-N), .SD[1], by = cell]

# STEP 2: Inter-cluster transitions only
inter <- cellint_data[source_cluster != target_cluster]
inter[, source_label := paste0("NN-", target_cluster, "_", source_cluster)]
inter[, target_label := paste0("NN-", source_cluster, "_", target_cluster)]

inter_labels <- rbind(
  inter[, .(cell = source_cell, label = source_label)],
  inter[, .(cell = target_cell, label = target_label)]
)

# Keep one label per cell
inter_labels_unique <- inter_labels[, .N, by = .(cell, label)][order(-N), .SD[1], by = cell]

# STEP 3: Combine, giving priority to intra-cluster labels
all_labels_unique <- rbindlist(list(
  intra_labels_unique,
  inter_labels_unique[!cell %in% intra_labels_unique$cell]
), use.names = TRUE)

# STEP 4: Add to Seurat metadata
obj$NN_transition <- NA_character_
matching_cells <- intersect(all_labels_unique$cell, colnames(obj))
obj$NN_transition[matching_cells] <- all_labels_unique[match(matching_cells, cell)]$label
table(obj$NN_transition)



objt = obj[, obj$NN_transition != "NA"]


library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gtools) 

library(ComplexHeatmap)
library(circlize)
library(gtools)

for (patient in patients) {
  cat("Processing patient:", patient, "\n")

  patient_cells <- colnames(objt)[objt$orig.ident == patient]
  if (length(patient_cells) < 2) next

  common <- get_high_conf_snps(patient_cells)
  if (length(common) < 2) next

  snp_matrix <- merged_matrix[common, patient_cells, drop = FALSE]
  snp_names <- rownames(snp_matrix)
  cell_names <- colnames(snp_matrix)
  
  celltypes_vec <- objt$NN_transition[patient_cells]
  celltypes <- factor(celltypes_vec, levels = unique(celltypes_vec))  # ensure factor
  
  # Extract chromosome from SNP names
  snp_chrom <- sub(":.*", "", snp_names)
  names(snp_chrom) <- snp_names

  # Hierarchical ordering by cell type
  ordered_cells <- unlist(lapply(levels(celltypes), function(ct) {
    ct_cells <- names(celltypes)[celltypes == ct]
    if (length(ct_cells) > 2) {
      hc <- hclust(dist(t(snp_matrix[, ct_cells, drop = FALSE])))
      ct_cells[hc$order]
    } else {
      ct_cells
    }
  }))

  snp_matrix <- snp_matrix[, ordered_cells]
  celltypes <- celltypes[ordered_cells]

  # Convert to dense matrix
  dense_mat <- as.matrix(snp_matrix)
  dense_mat <- t(dense_mat)  # Cells x SNPs

  # Extract SNP info
  snp_ids <- colnames(dense_mat)
  snp_info <- do.call(rbind, strsplit(snp_ids, "[:-]"))
  snp_df <- data.frame(
    chr = snp_info[,1],
    pos = as.numeric(snp_info[,2]),
    row.names = snp_ids
  )

  # Filter and order chromosomes
  valid_chrs <- paste0("chr", 1:22)
  snp_df <- snp_df[snp_df$chr %in% valid_chrs, ]
  snp_df$chr <- factor(snp_df$chr, levels = mixedsort(unique(snp_df$chr)))
  snp_df <- snp_df[order(snp_df$chr, snp_df$pos), ]

  # Reorder matrix
  dense_mat_ordered <- dense_mat[, rownames(snp_df), drop = FALSE]

  # Chromosome and cell type splits
  chrom_split <- factor(snp_df$chr, levels = levels(snp_df$chr))
  cell_split <- factor(celltypes, levels = levels(celltypes))

  # Cell annotation as data.frame
  cell_anno_df <- data.frame(CellType = celltypes)
  col_anno_colors <- list(CellType = setNames(
    scales::hue_pal()(length(levels(celltypes))), levels(celltypes)
  ))

  cell_anno <- rowAnnotation(df = cell_anno_df, col = col_anno_colors)

  # Heatmap
  ht <- Heatmap(
    dense_mat_ordered,
    name = "SNP",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_split = chrom_split,
    row_split = cell_split,
    left_annotation = cell_anno,
    use_raster = TRUE,
    raster_device = "png",
    bottom_annotation = HeatmapAnnotation(
      Chr = anno_block(
        gp = gpar(fill = NA),
        labels = levels(chrom_split),
        labels_gp = gpar(fontsize = 10)
      )
    ),
    gap = unit(1, "mm"),
    column_gap = unit(1, "mm"),
    row_gap = unit(1, "mm"),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 1)
  )

  # Save to PNG
  png_filename <- paste0("SNP_heatmap_", patient, ".png")
  png(png_filename, width = 5000, height = 2000, res = 400, bg = "black")
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}


