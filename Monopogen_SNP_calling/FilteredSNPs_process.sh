source("/mnt/raidbio/tmp/Alejandro/functions.R")
#saveRDS(obj.in, file = "/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
obj.in=readRDS("/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
new.cluster.ids <- c("SMCs", "MPs","ECs", "T-cells","PCs")
names(new.cluster.ids) <- levels(obj.in)
obj.in <- RenameIdents(obj.in, new.cluster.ids)
obj.in$CellType <- Idents(obj.in)



SRX23809079=obj.in[,obj.in$orig.ident=="Asymptomatic_SRX23809079"]

# Set the directory containing your CSV files
setwd("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/out_SRX23809079/somatic/")  # Change this to your actual path

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

colnames(merged_rds) <- c(colnames(merged_rds[,1:18]), paste0("Asymptomatic_SRX23809079_", colnames(merged_rds[,19:ncol(merged_rds)])))


saveRDS(merged_rds, file= "merged.SNV_mat.RDS")



setwd("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/")  # Change this to your actual path
source("/mnt/raidbio/tmp/Alejandro/functions.R")
#saveRDS(obj.in, file = "/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
obj.in=readRDS("/mnt/raidbio/tmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
new.cluster.ids <- c("SMCs", "MPs","ECs", "T-cells","PCs")
names(new.cluster.ids) <- levels(obj.in)
obj.in <- RenameIdents(obj.in, new.cluster.ids)
obj.in$CellType <- Idents(obj.in)



SRX23809079=obj.in[,obj.in$orig.ident=="Asymptomatic_SRX23809079"]
SRX23809079[,SRX23809079$CellType=="T-cells"]$orig.ident
table(SRX23809079$CellType)

meta<- read.csv(file="/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/out_SRX23809079/somatic/merged.putativeSNVs.csv")
head(meta)

mat <- readRDS(file="/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/out_SRX23809079/somatic/merged.SNV_mat.RDS")

meta_filter <-meta[meta$Depth_ref>5 & meta$Depth_alt>5,]
meta_filter <- meta_filter[meta_filter$SVM_pos_score>0.1,]
snv_pass <- paste0(meta_filter$chr,":", meta_filter$pos, ":", meta_filter$Ref_allele,":", meta_filter$Alt_allele)
mat <- mat[snv_pass,]


overlap <- intersect(colnames(mat),colnames(SRX23809079))
obj <- subset(SRX23809079,cells=overlap)
mat <- mat[,colnames(obj)]
# Convert data.frame to matrix for easier manipulation
mat <- as.matrix(mat)

# Extract numbers using regular expressions
split_values <- do.call(rbind, strsplit(mat, "/", fixed = TRUE))

# Convert extracted values to numeric
left_values <- as.numeric(split_values[, 1])
right_values <- as.numeric(split_values[, 2])

#Adapt SNV matrix 
mat[left_values== 0] <- 1
mat[right_values== 0] <- (-1)
mat[right_values== 0 & left_values== 0] <- 0
mat[right_values > 0] <- 1


# Convert back to data.frame if needed
mat <- as.data.frame(mat)

# Print the first few rows to check
head(mat)

Geno= CreateAssayObject(counts = mat)
obj[["Geno"]] <- Geno


####################################################################### Single SNP discovery

#####
target_cluster <- "T-cells"  # Replace with your cluster ID

# Extract SNP assay matrix
snp_matrix <- obj[["Geno"]]@counts


# Get cluster assignments
cluster_assignments <- obj$CellType

# Get cells belonging to the target cluster
cells_in_cluster <- colnames(obj)[cluster_assignments == target_cluster]

# Subset SNP matrix for those cells
snp_matrix_cluster <- snp_matrix[, cells_in_cluster, drop = FALSE]

# Select rows that do not contain any "1"
snps_1_in_cluster <- rownames(snp_matrix_cluster[apply(snp_matrix_cluster, 1, function(x) all(x != 1)), ])
snps_2_in_cluster <- rownames(snp_matrix_cluster[snps_1_in_cluster, ][rowSums(snp_matrix_cluster[snps_1_in_cluster, ]) < -15, ])

# Count occurrences of 1 in each row
counts <- rowSums(snp_matrix[snps_2_in_cluster, ] == 1)

# Sort row names by count in decreasing order
top_snps <- names(sort(counts, decreasing = TRUE))[1:20]

#Plot heatmap
library(viridis)
p=DoHeatmap(obj, features = top_snps, assay="Geno", group.by = "CellType", slot="counts")+ scale_fill_viridis()
save_plot(p, "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/SNP_heatmap_AtheroSRX23809079_Monopogen", 25, 10)


SNV_ex <- top_snps[5] 
geno_mat <- GetAssayData(obj, assay = "Geno")
obj@meta.data$Mutec <- geno_mat[SNV_ex, ]
cell_ref <- colnames(obj)[obj@meta.data$Mutec%in%c(-1)]
cell_alt <- colnames(obj)[obj@meta.data$Mutec%in%c(1)]

p1 <- DimPlot(obj,group.by = "CellType", pt.size=1, label=T) + NoLegend()
p2 <- DimPlot(obj ,cells.highlight = list("Ref"= cell_ref, "Alt"= cell_alt),  
    cols.highlight = c("blue", "red"), cols = "gray", pt.size = 1, sizes.highlight = 1) + 
    #annotate("text", x = 5, y = -12, label = paste("p-value: ", rest[rest$SNV == SNV_ex,]$pval, sep=""), color = "black", size = 6) + 
    ggtitle(paste0("Patient SRX23809079, SNV ", SNV_ex)) +
        theme(plot.title = element_text(hjust = 0.5))


p <- ggarrange(p1,p2)
save_plot(p, "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/MonopogenSNV_SRX23809079_4", 16, 6)

#obj[,obj$CellType=="Tcells"]


############## All single somatic SNVs
library(patchwork)
geno_mat <- GetAssayData(obj, assay = "Geno")

snv_list <- top_snps[c(1,2,3,5)] 

plot_list <- lapply(snv_list, function(SNV_ex) {
  if (!SNV_ex %in% rownames(geno_mat)) {
    warning(paste("SNV not found:", SNV_ex))
    return(NULL)
  }

  obj@meta.data$Mutec <- as.vector(geno_mat[SNV_ex, ])

  cell_ref <- colnames(obj)[obj@meta.data$Mutec %in% c(-1)]
  cell_alt <- colnames(obj)[obj@meta.data$Mutec %in% c(1)]

  DimPlot(
    obj,
    cells.highlight = list("Ref" = cell_ref, "Alt" = cell_alt),
    cols.highlight = c("blue", "red"), cols = "gray",
    pt.size = 1.5, sizes.highlight = 2
  ) +
    ggtitle(paste0("Patient SRX23809079, SNV ", SNV_ex)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
})

# Remove failed plots
plot_list <- Filter(Negate(is.null), plot_list)

# Combine and save
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 2)
ggsave("SNV_DimPlots_SRX23809079.png", combined_plot, width = 8*2, height = 6*2, dpi = 600)








############Permutation test (null distribution)
# Your observed count of interesting SNPs
num_detections=sort(counts, decreasing = TRUE)[20]
observed_count <- sum(rowSums(snp_matrix[snps_2_in_cluster, ] == 1) >= num_detections)


set.seed(123)  # for reproducibility
n_permutations <- 1000
perm_counts <- numeric(n_permutations)

all_cells <- colnames(obj)
true_t_cells <- colnames(obj)[obj$CellType == "T-cells"]
n_tcells <- length(true_t_cells)

for (i in 1:n_permutations) {
  # Randomly sample the same number of T-cells
  random_cells <- sample(all_cells, n_tcells)
  
  # Subset SNP matrix
  snp_cluster <- snp_matrix[, random_cells, drop = FALSE]
  
  # SNPs not carrying '1' (mutation) in the shuffled T-cells
  snps_without_1 <- rownames(snp_cluster[apply(snp_cluster, 1, function(x) all(x != 1)), ])
  
  # Among those, select SNPs with reference support in those cells
  snps_with_high_ref <- rownames(snp_cluster[snps_without_1, ][rowSums(snp_cluster[snps_without_1, ]) < -10, ])
  
  # Now count how many of these SNPs are altered in other cells
  count <- rowSums(snp_matrix[snps_with_high_ref, ] == 1)
  
  # Keep count of SNPs that pass
  perm_counts[i] <- sum(count >= num_detections)
}

# Empirical p-value
empirical_p <- mean(perm_counts >= observed_count)
cat("Empirical p-value:", empirical_p, "\n")

# Convert permutation count to a data frame
perm_df <- data.frame(SNP_count = perm_counts)


# Histogram
png("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Histogram_permutation_test_SRX23809079.png", 
    width = 2000, height = 1800, res = 300)

ggplot(perm_df, aes(x = SNP_count)) +
  geom_histogram(bins = 50, fill = "#d0d0d0", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = observed_count, color = "#d73027", linetype = "dashed", size = 1) +
  annotate("text",
           x = observed_count,
           y = max(table(cut(perm_counts, breaks = 50))) * 0.95,
           label = paste0("Empirical p = ", round(empirical_p, 4)),
           hjust = 1.2, color = "#d73027", size = 5, fontface = "bold") +
  labs(
    title = "Permutation Test: Null Distribution of SNP Counts",
    subtitle = "Random sampling of T-cells vs. non-T cell SNP distribution",
    x = "SNPs detected",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

dev.off()




######################################################################Pairwise Jaccard index
#########################################################

snp_matrix <- obj[["Geno"]]@counts
cell_types <- obj$CellType
clusters <- unique(as.character(cell_types))  # Exclude these clusters

# Function to get SNPs detected in >80% of cells for a given set of cells
get_high_conf_snps <- function(cells, threshold = 0.8) {
  mat <- snp_matrix[, cells, drop = FALSE]
  detected_fraction <- rowMeans(mat != 0)
  rownames(mat)[detected_fraction >= threshold]
}

# Function to get SNPs present (value == 1) in at least 5 cells
get_present_snps <- function(cells, snp_subset) {
  mat <- snp_matrix[snp_subset, cells, drop = FALSE]
  present <- rowSums(mat == 1) >= 5
  rownames(mat)[present]
}

# Initialize matrices
jaccard_matrix <- matrix(0, nrow = length(clusters), ncol = length(clusters),
                         dimnames = list(clusters, clusters))
label_matrix <- matrix("", nrow = length(clusters), ncol = length(clusters),
                       dimnames = list(clusters, clusters))

# Pairwise Jaccard calculation with annotation
for (i in clusters) {
  cells_i <- colnames(snp_matrix)[cell_types == i]
  for (j in clusters) {
    cells_j <- colnames(snp_matrix)[cell_types == j]

    # High-confidence SNPs
    snps_i <- get_high_conf_snps(cells_i)
    snps_j <- get_high_conf_snps(cells_j)
    snps_shared <- intersect(snps_i, snps_j)

    # SNPs present (>= 5 cells with value 1)
    snps_present_i <- get_present_snps(cells_i, snps_shared)
    snps_present_j <- get_present_snps(cells_j, snps_shared)

    inter <- length(intersect(snps_present_i, snps_present_j))
    union <- length(union(snps_present_i, snps_present_j))

    jaccard_matrix[i, j] <- ifelse(union > 0, inter / union, NA)
    label_matrix[i, j] <- ifelse(union > 0, paste0(inter, "/", union), "")
  }
}


# Plot with annotations
png("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Jaccard_SNP_similarity_index_SRX23809079.png", 
    width = 2600, height = 1200, res = 300)
pheatmap(jaccard_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = label_matrix,
         main = "Jaccard Similarity Index for Patient SRX23809079")
dev.off()







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
patients <- unique(obj$orig.ident)


get_high_conf_snps <- function(cells, threshold = 0.8) {
  mat <- merged_matrix[, cells, drop = FALSE]
  detected_fraction <- rowMeans(mat != 0)
  rownames(mat)[detected_fraction >= threshold]
}



# Store results
all_alluvial_data <- list()

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

#
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
   if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
 
color_list <- ggplotColours(n=5)

# Summarize total transitions
label_data <- combined_alluvial[, .(total_N = sum(N)), by = .(source_cluster, target_cluster)]
label_data[, label := total_N * 2]

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
    title = "High-confidence SNV-based NN-cells (99th percentile) ",
    y = "Number of Cells", x = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

save_plot(p, "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Alluvial_80common_PatientWise_Filtered", 8, 8)



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
stats_dt <- shared_snp_dt[, .(percentile99 = quantile(shared_snp, probs = 0.99)), by = patient]



# Plot with mean + SD line
# Corrected percentile plot
dist_plot <- ggplot(shared_snp_dt, aes(x = shared_snp)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", boundary = 0) +
  facet_wrap(~ patient, scales = "free_y", ncol = 3) +
  geom_vline(data = stats_dt, aes(xintercept = percentile99), 
             color = "red", linetype = "dashed", size = 1) +
  labs(title = "Shared SNPs per Cell Pair (99th percentile threshold per patient)",
       x = "Number of Shared SNPs", y = "Count of Cell Pairs") +
  theme_minimal(base_size = 14)


# Save it
ggsave(
  filename = "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/SNVpair_distribution_by_patient_mean_sd.png",
  plot = dist_plot, width = 9, height = 6, dpi = 300
)






################Check macrophage ontogeny differences

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
  threshold_99 <- quantile(max_shared_per_cell, 0.95)

  # Select cells above threshold
  filtered_cells <- names(max_shared_per_cell[max_shared_per_cell >= threshold_99])
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

color_list <- ggplotColours(n=4)

p <- DimPlot(
  obj,
  cells.highlight = list(
    "NN-SMCs_MPs" = colnames(obj[, obj$NN_transition == "NN-SMCs_MPs"]),
    "NN-MPs_MPs"   = colnames(obj[, obj$NN_transition == "NN-MPs_MPs"])
  ),
  cols.highlight = color_list,
  cols = "gray",
  shuffle = FALSE,
  pt.size = 1,
  sizes.highlight = 1
) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  annotate("text", x = 5, y = 10, label = "NN-MPs_MPs: 92 cells", size = 5) +
  annotate("text", x = 5, y = 9, label = "NN-SMCs_MPs: 8 cells", size = 5)

save_plot(
  p,
  "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos/DEmarkers/NN-macrophages_UMAP",
  fig.width = 8, fig.height = 6
)

####################

# Step 1: Tabulate and filter NN_transition categories with >3 cells
nn_counts <- table(obj$NN_transition)
valid_transitions <- names(nn_counts[nn_counts > 3])

# Step 2: Get cells belonging to those transitions
cells_to_keep <- WhichCells(obj, expression = NN_transition %in% valid_transitions)

# Step 3: Subset the Seurat object
obj_filtered <- subset(obj, cells = cells_to_keep)

table(obj_filtered$NN_transition)


##############
deResTT = makeDEResults(obj_filtered, group.by="NN_transition", assay="RNA", test="t")
write.table(deResTT,"/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos/DEmarkers/expr_test_clusters_t.tsv", sep="\t", row.names=F, quote = F)
write_xlsx(deResTT, "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos/DEmarkers/expr_test_clusters_t.xlsx")


#deResTT<-read_tsv("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos/DEmarkers/expr_test_clusters_t.tsv")
DefaultAssay(obj_filtered) <- "RNA"
markers.use.tt= subset(deResTT , avg_log2FC>0&p_val_adj<0.05&!startsWith(gene, "mt-")&!startsWith(gene, "rp"))
finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:10) #subset(clusterID==4)
finalMarkers.use.tt


data_dupli= finalMarkers.use.tt[!duplicated(finalMarkers.use.tt[ , "gene"]), ]


events= data_dupli %>% count(clusterID)
inser=cumsum(events$n)+0.5-events$n
insert=replace(inser, inser==0.5, 0)

xmi<- insert
xmin<- xmi[c(FALSE, TRUE)]
xma<- insert+events$n
xmax<- xma[c(FALSE, TRUE)]
ymi<- 0*events$n
ymin<- ymi[c(FALSE, TRUE)]
yma<- rep(length(events$n)+0.5, each=length(events$n))
ymax<- yma[c(FALSE, TRUE)]



p_dp_genes_idents = DotPlot(obj_filtered, features = data_dupli$gene, assay="RNA", dot.scale = 5, group.by="NN_transition")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
    annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin , ymax=ymax, alpha=0.2, fill="blue") #rep(c("blue", "grey"), times= length(events$n)/2)
save_plot(p_dp_genes_idents, "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos/DEmarkers/dotplot_cluster_genes_colored", 35, 8)























####################
deGroup = "NN_transition"


allConds = unique(obj$NN_transition)[unique(obj$NN_transition) %in% c("NN-MPs_MPs","NN-SMCs_MPs")]


print("all conds")
print(allConds)

all.cells = list()

for (cond in allConds)
{

  all.cells[[cond]] = cellIDForClusters(obj, deGroup, c(cond))

}

retList = list()
retListW= list()
for (i in 1:(length(allConds)-1))
{

  for (j in (i+1):length(allConds))
  {
    print(paste(i,"<->",j))

    condI = allConds[i]
    condJ = allConds[j]

    condNameI = tolower(str_replace_all(condI, " ", "_"))
    condNameJ = tolower(str_replace_all(condJ, " ", "_"))

    print(paste(condJ, condNameJ))
    print(paste(condI, condNameI))

    DEmarkers = compareClusters(obj, all.cells[[condJ]], all.cells[[condI]], condNameJ, condNameI,
                                                outfolder=paste("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos", "de", sep="/"), fcCutoff=0.1)


    DEmarkersW = compareClusters(obj, all.cells[[condJ]], all.cells[[condI]], condNameJ, condNameI,
                                                    outfolder=paste("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos", "dewilcox", sep="/"), test="wilcox", fcCutoff=0.1)

    compName = paste(condNameJ, condNameI, sep="_")
        
    retList[[compName]] = DEmarkers
    retListW[[compName]] = DEmarkersW
  }
}




makeVolcanos(retList, deGroup, paste("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos", "de_volcano", sep="/"),
                turnExpression=F, FCcutoff=0.1, pCutoff = 0.05, width = 8, height = 7)

makeVolcanos(retListW, deGroup, paste("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos", "dewilcox_volcano",  sep="/"),
                turnExpression=F, FCcutoff=0.1, pCutoff = 0.05, width = 8, height = 7)





##############
obj.M= obj[,obj$NN_transition %in% c("NN-MPs_MPs","NN-SMCs_MPs")]

markers.use.tt <- retList[[1]] %>%
  filter(
    (avg_log2FC > 0.1 | avg_log2FC < -0.1) &
    p_val_adj < 0.05
)

# Rank by significance and effect size
finalMarkers.use.tt <- markers.use.tt %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) 
finalMarkers.use.tt$gene

top_genes <- finalMarkers.use.tt[!grepl("^ENS", finalMarkers.use.tt$gene), ][1:25, ]$gene

p_dp_genes_idents <- DotPlot(
  obj.M,
  features = top_genes,
  assay = "RNA",
  dot.scale = 6,
  col.min = 0,
  group.by = "NN_transition"
) +
  coord_flip() +
  scale_color_gradient(low = "lightgrey", high = "firebrick") +  # Better contrast
  scale_size(range = c(1, 6)) +  # Control dot size range
  labs(
    title = "SMC-adjacent MPs vs MP-adjacent MPs",
    x = NULL, y = "Gene Expression"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "right"
  )

# Save plot
save_plot(
  p_dp_genes_idents,
  "/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/volcanos/DEmarkers/dotplot_NN-MPs",
  fig.width = 5, fig.height = 7
)

