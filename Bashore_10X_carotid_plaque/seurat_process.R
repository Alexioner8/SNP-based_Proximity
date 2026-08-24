
######

source("/mnt/raidexttmp/Alejandro/functions.R")
setwd("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs")
library(SoupX)

#
# this searches for all input matrices
#
files=list()
for (x in c("carotid_4_Symptomatic", "carotid_6_Asymptomatic", "carotid_9_Symptomatic", "carotid_13_Symptomatic", "carotid_14_Symptomatic", "carotid_15_Asymptomatic", "carotid_16_Asymptomatic", "carotid_17_Asymptomatic"))
{
file<- Sys.glob(paste(paste("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/", x, sep=""), "/raw_feature_bc_matrix.h5", sep=""))

files[x]= file
}



#
# here we iterate over all input matrices and load them into allfiles (gene expression) and allABs (antibody capture)
#
allfiles.raw = list()
for (file in files)
{
  samplename = str_split(dirname(file), "/")[[1]][6] #estos numeros indican el nombre de la carpeta
  
  print(samplename)
  
  h5file = Read10X_h5(file,unique.features = TRUE)

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
    # the patternlist.human contains patterns for MT and RP-genes
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





objlist <- lapply(X = objlist, FUN = function(obj) {
  # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
  print(paste("Seurat obj project", obj@project.name))
  print(obj)
  obj <- subset(obj, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA > 500)
  obj <- subset(obj, subset = percent.mt < 15)
  print(obj)
  
  return(obj)
})



for (name in names(objlist))
{
  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  save_plot(p, paste("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/QC", paste(name, "filtered_violins_qc", sep="_"), sep="/"), fig.width=10, fig.height=6)

  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine=F)
  p[[1]] = p[[1]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[2]] = p[[2]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[3]] = p[[3]] + scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,5))
  p = combine_plot_grid_list(plotlist=p, ncol=3)
  save_plot(p, paste("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/QC", paste(name, "filtered_violins_detail_qc", sep="_"), sep="/"), fig.width=18, fig.height=6)
  
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/QC", paste(name, "filtered_scatter_ncount_mt", sep="_"), sep="/"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/QC", paste(name, "filtered_scatter_ncount_rp", sep="_"), sep="/"), fig.width=10, fig.height=6)
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



    x <- ScaleData(x, features = features, verbose = FALSE, assay="RNA", vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA"))
    x <- RunPCA(x, verbose = FALSE, reduction.name="pca", assay="RNA")
    #x <- suppressWarnings(SCTransform(x, verbose = FALSE,vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score")))
    
    plot1 <- ElbowPlot(x)
    save_plot(plot1, paste("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/QC", paste(x@project.name, "elbowplot_dimensionality", sep="_"), sep="/"), fig.width=10, fig.height=6)


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

features_gex <- SelectIntegrationFeatures(object.list = objSamples, nfeatures = 2000)#, assay=rep("RNA", length(objSamples)))
#objSamples <- PrepSCTIntegration(object.list = objSamples, anchor.features = features_gex)
objlist.anchors <- FindIntegrationAnchors(object.list = objSamples,  reduction = "rpca", dims = 1:20, anchor.features = features_gex) #normalization.method = "SCT",
obj.list.integrated <- IntegrateData(new.assay.name = "integratedgex", anchorset = objlist.anchors, dims = 1:20, verbose=T) #normalization.method = "SCT",
print("GEX integration done")

#
# integrated GEX viz
#
obj.list.integrated = ScaleData(obj.list.integrated, assay="integratedgex")
obj.list.integrated <- RunPCA(obj.list.integrated, reduction.name="igpca", assay="integratedgex")
obj.list.integrated <- RunUMAP(obj.list.integrated, reduction = "igpca", dims = 1:20, reduction.key = "UMAP",)
p=DimPlot(obj.list.integrated, reduction="umap", shuffle = T, seed = 1, group.by= "orig.ident")
save_plot(p, paste(intname, "wnn_ig_dimplot", sep="/"), 12, 8)

p=DimPlot(obj.list.integrated, reduction="igpca", group.by= "orig.ident")
save_plot(p, paste(intname, "wnn_pca_ig_dimplot", sep="/"), 12, 8)





return(obj.list.integrated)

}


finalList_sample = prepareFinalList(objlist)
integratedList_plaque = analyseFinalList(finalList_sample, "Seurat/HumanPlaque_int")



saveRDS(integratedList_plaque, file = "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/bashoreint.rds")
#integratedList_plaque=readRDS("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/bashoreint.rds")


obj.in <- FindVariableFeatures(integratedList_plaque, selection.method = "vst")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj.in), 30)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj.in) 
p= LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(p, "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/HV_feat", 10, 10)


#cluster identification
DefaultAssay(integratedList_plaque) = "integratedgex"
obj.in <- FindNeighbors(integratedList_plaque, reduction="igpca", dims = 1:20)
DefaultAssay(obj.in) <- "integratedgex"
obj.in <- FindClusters(obj.in, resolution = 0.2, algorithm = 4)

p=DimPlot(obj.in, pt.size = 1.5, label=T, reduction = "umap")
save_plot(p, "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/Dimplot_umap_ident", fig.width=10, fig.height=8)


############
# Define the sample IDs
samples <- c("carotid_4_Symptomatic", "carotid_6_Asymptomatic", "carotid_9_Symptomatic", "carotid_13_Symptomatic", "carotid_14_Symptomatic", "carotid_15_Asymptomatic", "carotid_16_Asymptomatic", "carotid_17_Asymptomatic")

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
  output_file <- paste0("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/", x, "_cell_reads.csv")
  
  # Write to CSV
  write.csv(output_data, file = output_file, row.names = FALSE, quote = FALSE)
}


DefaultAssay(obj.in) <- "RNA"
obj.in = JoinLayers(obj.in)
###########
deResTT = makeDEResults(obj.in, group.by="seurat_clusters", assay="RNA", test="wilcox")
write.table(deResTT,"/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/expr_test_clusters_t.tsv", sep="\t", row.names=F, quote = F)
write_xlsx(deResTT, "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/expr_test_clusters_t.xlsx")


#deResTT<-read_tsv("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/expr_test_clusters_t.tsv")
DefaultAssay(obj.in) <- "RNA"
markers.use.tt= subset(deResTT , avg_log2FC>0&p_val_adj<0.05&!startsWith(gene, "mt-")&!startsWith(gene, "rp"))

# Ensure clusters are ordered numerically
# Convert clusterID to numeric-ordered factor
markers.use.tt$clusterID <- factor(
  markers.use.tt$clusterID,
  levels = sort(as.numeric(unique(markers.use.tt$clusterID)))
)

finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:10) 
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

 


p_dp_genes_idents = DotPlot(obj.in, features = data_dupli$gene, assay="RNA", dot.scale = 5, group.by="seurat_clusters")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
    annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin , ymax=ymax, alpha=0.2, fill="blue") #rep(c("blue", "grey"), times= length(events$n)/2)
save_plot(p_dp_genes_idents, "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/dotplot_cluster_genes_colored", 30, 8)


## Extract the list of genes for clusterID
genes <- finalMarkers.use.tt[finalMarkers.use.tt$clusterID == 14, ]$gene

# Collapse the gene list into a comma-separated string
genes_comma_separated <- paste(genes, collapse = ",")

# Write the string to a file
write(genes_comma_separated, file = "genes_clusterID.txt")


#############
cellList = colnames(obj.in)

featVec <- vector(mode="character", length=length(cellList))
featVec = as.character(obj.in$seurat_clusters)


featVec[featVec == "1"] = "Fibroblasts 1"
featVec[featVec == "2"] = "ECs"
featVec[featVec == "3"] = "MPs"
featVec[featVec == "4"] = "SMCs"
featVec[featVec == "5"] = "T-cells"
featVec[featVec == "6"] = "B-cells"
featVec[featVec == "7"] = "SMC foam cells"
featVec[featVec == "8"] = "modulated SMCs"
featVec[featVec == "9"] = "Peripheral nerve cells"
featVec[featVec == "10"] = "NK cells"
featVec[featVec == "11"] = "Plasma cells"
featVec[featVec == "12"] = "Monocytes"
featVec[featVec == "13"] = "Fibroblast 2"
featVec[featVec == "14"] = "Mast cells"


obj.in$label=featVec



#~/R.sh

source("/mnt/raidexttmp/Alejandro/functions.R")
setwd("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat")
#saveRDS(obj.in, file = "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/bashoreint.rds")
#obj.in=readRDS("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/bashoreint.rds")



p=DimPlot(obj.in, pt.size = 0.5, group.by="label", reduction = "umap", label=T)
save_plot(p, paste("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat", "Dimplot_umap_label", sep="/"), 12, 8)

samples <- unique(obj.in$orig.ident)

for (s in samples) {
  
  # Subset object
  obj.sub <- subset(obj.in, subset = orig.ident == s)
  
  # Create plot
  p <- DimPlot(
    obj.sub,
    pt.size = 0.5,
    group.by = "label",
    reduction = "umap",
    label = TRUE
  ) + ggtitle(s)
  
  # Define filename (safe formatting)
  fname <- file.path(
    "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat",
    paste0("Dimplot_umap_label_", s, ".png")
  )
  
  # Save
  ggsave(filename = fname, plot = p, width = 12, height = 8)
}


##########################################
######################################### Somatic SNVs Monopogen


samples <- unique(obj.in$orig.ident)

base_dir <- "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Monopogen_out"

for (s in samples) {
  
  message("Processing: ", s)
  
  # Define input directory
  in_dir <- file.path(base_dir, paste0("out_", s, "somatic"))
  in_dir <- gsub("out_", "out_", in_dir)  # safe placeholder if needed
  
  # FIXED: correct path construction (your original missed a slash)
  in_dir <- file.path(base_dir, paste0("out_", s), "somatic")
  
  setwd(in_dir)
  
  #### -------------------------
  #### 1. Merge CSV SNVs
  #### -------------------------
  
  file_list <- list.files(pattern = "chr[0-9]+\\.putativeSNVs\\.csv")
  
  merged_data <- file_list %>%
    lapply(read.csv) %>%
    bind_rows()
  
  write.csv(
    merged_data,
    file = file.path(in_dir, "merged.putativeSNVs.csv"),
    row.names = FALSE
  )
  
  #### -------------------------
  #### 2. Merge RDS SNV matrices
  #### -------------------------
  
  rds_list <- list.files(pattern = "chr[0-9]+\\.SNV_mat\\.RDS")
  
  merged_rds <- rds_list %>%
    lapply(readRDS) %>%
    bind_rows()
  
  # rename columns (sample-specific prefix)
  colnames(merged_rds) <- c(
    colnames(merged_rds[, 1:18]),
    paste0(s, "_", colnames(merged_rds[, 19:ncol(merged_rds)]))
  )
  
  saveRDS(
    merged_rds,
    file = file.path(in_dir, "merged.SNV_mat.RDS")
  )
}


###############


carotid_6_Asymptomatic=obj.in[,obj.in$orig.ident == "carotid_6_Asymptomatic"]

table(carotid_6_Asymptomatic$label)

meta<- read.csv(file="/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Monopogen_out/out_carotid_6_Asymptomatic/somatic/merged.putativeSNVs.csv")
head(meta)

mat <- readRDS(file="/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Monopogen_out/out_carotid_6_Asymptomatic/somatic/merged.SNV_mat.RDS")

meta_filter <-meta[meta$Depth_ref>5 & meta$Depth_alt>5,]
meta_filter <- meta_filter[meta_filter$SVM_pos_score>0.1,]
snv_pass <- paste0(meta_filter$chr,":", meta_filter$pos, ":", meta_filter$Ref_allele,":", meta_filter$Alt_allele)
mat <- mat[snv_pass,]


overlap <- intersect(colnames(mat),colnames(carotid_6_Asymptomatic))
obj <- subset(carotid_6_Asymptomatic,cells=overlap)
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
cluster_assignments <- obj$label

# Get cells belonging to the target cluster
cells_in_cluster <- colnames(obj)[cluster_assignments == target_cluster]

# Subset SNP matrix for those cells
snp_matrix_cluster <- snp_matrix[, cells_in_cluster, drop = FALSE]

# Select rows that do not contain any "1"
keep <- apply(snp_matrix_cluster, 1, function(x) all(x != 1))
keep[is.na(keep)] <- FALSE
snp_subset <- snp_matrix_cluster[keep, , drop = FALSE]
snp_subset1 <- snp_subset[rowSums(snp_subset) < -30, ]

# Count occurrences of 1 in each row
counts <- rowSums(snp_matrix[rownames(snp_matrix) %in% rownames(snp_subset1),] == 1)

# Sort row names by count in decreasing order
top_snps <- names(sort(counts, decreasing = TRUE))[1:20]

#Plot heatmap
library(viridis)
p=DoHeatmap(obj, features = top_snps, assay="Geno", group.by = "label", slot="counts")+ scale_fill_viridis()
save_plot(p, "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/SNVs_bashore/SNP_heatmap_Atherocarotid_6_Asymptomatic", 25, 10)


SNV_ex <- top_snps[1] 
geno_mat <- GetAssayData(obj, assay = "Geno")
obj@meta.data$Mutec <- geno_mat[SNV_ex, ]
cell_ref <- colnames(obj)[obj@meta.data$Mutec%in%c(-1)]
cell_alt <- colnames(obj)[obj@meta.data$Mutec%in%c(1)]

p1 <- DimPlot(obj ,cells.highlight = list("Ref"= cell_ref, "Alt"= cell_alt),  
    cols.highlight = c("blue", "red"), cols = "gray", pt.size = 1, sizes.highlight = 1) + 
    #annotate("text", x = 5, y = -12, label = paste("p-value: ", rest[rest$SNV == SNV_ex,]$pval, sep=""), color = "black", size = 6) + 
    ggtitle(paste0("Patient_6_Asymptomatic, SNV ", SNV_ex)) +
        theme(plot.title = element_text(hjust = 0.5))


p2 <- DimPlot(obj, pt.size = 0.5, group.by="label", reduction = "umap", label=T)+
ggtitle("Bashore 3' Next GEM 10X Genomics, Patient_6_Asymptomatic") +
        theme(plot.title = element_text(hjust = 0.5))



p <- ggarrange(p1,p2)
save_plot(p, "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/SNVs_bashore/MonopogenSNV_patient_6_Asymptomatic", 12*2, 8)


##################### Load all SNV matrices into seurat object

libraries <- unique(obj.in$orig.ident)

snv_list <- list()

for (library_name in libraries) {

  message("Processing: ", library_name)

  # --- Load files ---
  mat_path <- paste0(
    "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Monopogen_out/out_",
    library_name,
    "/somatic/merged.SNV_mat.RDS"
  )

  meta_path <- paste0(
    "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Monopogen_out/out_",
    library_name,
    "/somatic/merged.putativeSNVs.csv"
  )

  if (!file.exists(mat_path) | !file.exists(meta_path)) {
    message("Skipping (missing files): ", library_name)
    next()
  }

  mat <- readRDS(mat_path)
  meta <- read.csv(meta_path)

  # --- Filter SNVs ---
  meta_filter <- meta[meta$Depth_ref > 5 & meta$Depth_alt > 5, ]
  meta_filter <- meta_filter[meta_filter$SVM_pos_score > 0.1, ]

  snv_pass <- paste0(meta_filter$chr, ":", meta_filter$pos, ":", 
                     meta_filter$Ref_allele, ":", meta_filter$Alt_allele)

  mat <- mat[snv_pass, , drop = FALSE]

  # --- Keep only overlapping cells ---
  overlap <- intersect(colnames(mat), colnames(obj.in))
  if (length(overlap) == 0) next

  mat <- mat[, overlap, drop = FALSE]

  # --- Convert genotype encoding ---
  mat <- as.matrix(mat)

  split_values <- do.call(rbind, strsplit(mat, "/", fixed = TRUE))

  left_values  <- as.numeric(split_values[, 1])
  right_values <- as.numeric(split_values[, 2])

  mat[left_values == 0] <- 1
  mat[right_values == 0] <- -1
  mat[right_values == 0 & left_values == 0] <- 0
  mat[right_values > 0] <- 1

  mat <- matrix(as.numeric(mat), nrow = nrow(mat))
  rownames(mat) <- snv_pass
  colnames(mat) <- overlap

  # --- Store ---
  library(Matrix)

  mat <- as.matrix(mat)
  mode(mat) <- "numeric"

  sparse_mat <- Matrix(mat, sparse = TRUE)

  snv_list[[library_name]] <- sparse_mat
}


# Collect all SNVs
all_snvs <- unique(unlist(lapply(snv_list, rownames)))

# Map SNVs to indices
snv_index <- setNames(seq_along(all_snvs), all_snvs)

# Preallocate lists for sparse construction
i_list <- list()
j_list <- list()
x_list <- list()

col_offset <- 0
cell_names <- c()

for (library_name in names(snv_list)) {

  mat <- snv_list[[library_name]]

  # Get indices of non-zero entries
  nz <- summary(mat)   # gives i, j, x

  # Map row indices to global SNV index
  global_i <- snv_index[rownames(mat)[nz$i]]

  # Shift column indices
  global_j <- nz$j + col_offset

  i_list[[library_name]] <- global_i
  j_list[[library_name]] <- global_j
  x_list[[library_name]] <- nz$x

  col_offset <- col_offset + ncol(mat)
  cell_names <- c(cell_names, colnames(mat))
}

# Combine all entries
i <- unlist(i_list)
j <- unlist(j_list)
x <- unlist(x_list)

# Create sparse matrix
global_mat <- sparseMatrix(
  i = i,
  j = j,
  x = x,
  dims = c(length(all_snvs), length(cell_names))
)

rownames(global_mat) <- all_snvs
colnames(global_mat) <- cell_names

for (library_name in names(snv_list)) {

  mat <- snv_list[[library_name]]

  global_mat[rownames(mat), colnames(mat)] <- mat
}


# Feature vs SNPs plot
mut_snp_per_cell <- colSums(global_mat == 1)
obj.in[["totalSNPs"]] <- mut_snp_per_cell


Geno <- CreateAssayObject(counts = global_mat)
obj.in[["Geno"]] <- Geno


saveRDS(obj.in, file = "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/bashoreint_withSNVs.rds")



##########


#~/R.sh

source("/mnt/raidexttmp/Alejandro/functions.R")
setwd("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat")

#obj.in=readRDS("/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/bashoreint_withSNVs.rds")
#smart_matrix <- readRDS("/mnt/raidexttmp/Alejandro/Monopogen_Athero_Smartseq2/Allpatients.SNV_sparse_mat_merged.RDS")


library(Matrix)

geno <- GetAssayData(obj.in, assay = "Geno", slot = "counts")

get_snps_in_cells <- function(mat, cells, threshold = 0.8) {

  submat <- mat[, cells, drop = FALSE]

  frac <- Matrix::rowSums(submat != 0) / ncol(submat)

  frac[is.na(frac)] <- 0

  rownames(submat)[frac >= threshold]
}

patients <- unique(obj.in$orig.ident)

patient_housekeeping_snps <- list()

for (patient in patients) {

  cat("Processing:", patient, "\n")

  # Cells from this patient
  patient_cells <- colnames(obj.in)[obj.in$orig.ident == patient]

  if (length(patient_cells) < 2) next

  # Cell types in this patient
  celltypes <- unique(obj.in$label[patient_cells])

  # Store SNPs per cluster
  cluster_snps <- list()

  for (ct in celltypes) {

    cells_ct <- patient_cells[obj.in$label[patient_cells] == ct]

    if (length(cells_ct) < 2) next

    snps_ct <- get_snps_in_cells(geno, cells_ct, threshold = 0.8)

    if (length(snps_ct) > 0) {
      cluster_snps[[ct]] <- snps_ct
    }
  }

  # 🔥 Key step: intersect across clusters
  if (length(cluster_snps) > 1) {
    common_snps <- Reduce(intersect, cluster_snps)
  } else if (length(cluster_snps) == 1) {
    common_snps <- cluster_snps[[1]]
  } else {
    next
  }

  patient_housekeeping_snps[[patient]] <- common_snps

  cat("  -> SNVs retained:", length(common_snps), "\n")
}


patient_housekeeping_snps

# Smart-seq2
smartseq2_snvs <- nrow(smart_matrix)/length(unique(obj.in$orig.ident))

# 10X (sum across patients)
tenx_snvs <- sum(sapply(patient_housekeeping_snps, length))

plot_df <- data.frame(
  Method = c("Smart-seq2_Mocci", "10X_Bashore"),
  SNVs = c(smartseq2_snvs, tenx_snvs)
)

p=ggplot(plot_df, aes(x = Method, y = SNVs)) +
  geom_bar(stat = "identity", fill = "#2C7FB8") +
  geom_text(aes(label = SNVs), vjust = -0.3, size = 5) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  ) +
  ylab("Average of HLF-sSNV loci per patient") +
  xlab("") +
  ggtitle("Housekeeping-like gene SNV loci present in 80% cells")
save_plot(p, "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/SNV_comparison_smartseq2_10x", 6, 6)

#################
library(Matrix)
library(ggplot2)

geno <- GetAssayData(obj.in, assay = "Geno", slot = "counts")

patients <- unique(obj.in$orig.ident)

thresholds <- seq(0.8, 0.2, by = -0.1)

get_snps_in_cells <- function(mat, cells, threshold) {

  submat <- mat[, cells, drop = FALSE]

  frac <- Matrix::rowSums(submat != 0) / ncol(submat)
  frac[is.na(frac)] <- 0

  rownames(submat)[frac >= threshold]
}


results <- data.frame()

for (th in thresholds) {

  cat("Processing threshold:", th, "\n")

  total_snps <- 0

  for (patient in patients) {

    patient_cells <- colnames(obj.in)[obj.in$orig.ident == patient]
    if (length(patient_cells) < 2) next

    celltypes <- unique(obj.in$label[patient_cells])

    cluster_snps <- list()

    for (ct in celltypes) {

      cells_ct <- patient_cells[obj.in$label[patient_cells] == ct]

      if (length(cells_ct) < 2) next

      snps_ct <- get_snps_in_cells(geno, cells_ct, threshold = th)

      if (length(snps_ct) > 0) {
        cluster_snps[[ct]] <- snps_ct
      }
    }

    # Intersect across clusters
    if (length(cluster_snps) > 1) {
      common_snps <- Reduce(intersect, cluster_snps)
    } else if (length(cluster_snps) == 1) {
      common_snps <- cluster_snps[[1]]
    } else {
      next
    }

    total_snps <- total_snps + length(common_snps)
  }

  results <- rbind(results, data.frame(
    threshold = th,
    total_snps = total_snps
  ))
}


p=ggplot(results, aes(x = factor(threshold), y = total_snps)) +
  geom_bar(stat = "identity", fill = "#2C7FB8") +  # nice scientific blue
  theme_classic()+
  ggtitle("Somatic SNV loci detected at different cutoffs") +
  xlab("Detection threshold (fraction of cells cluster-wise)") +
  ylab("Number of somatic SNV loci") 
save_plot(p, "/mnt/raidexttmp/Alejandro/bashore_human_plaque_SNVs/Seurat/housekeeping_snvs_threshold_sensitivity", 6, 6)


###############

