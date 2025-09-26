
######
source("/mnt/raidexttmp/Alejandro/functions.R")
# Read the counts file

#
# this searches for all input matrixrices
#
filesasym= list()
for (file in c("SRX23809072_bam_counts.txt","SRX23809073_bam_counts.txt","SRX23809074_bam_counts.txt", "SRX23809079_bam_counts.txt", "SRX23809086_bam_counts.txt","SRX23809085_bam_counts.txt","SRX23809081_bam_counts.txt"))
{
filesasym[file]<- Sys.glob(paste("/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/featureCounts/", file, sep=""))
}

filessym= list()
for (file in c("SRX23809075_bam_counts.txt", "SRX23809078_bam_counts.txt", "SRX23809076_bam_counts.txt", "SRX23809077_bam_counts.txt","SRX23809083_bam_counts.txt", "SRX23809080_bam_counts.txt", "SRX23809082_bam_counts.txt", "SRX23809084_bam_counts.txt"))
{
filessym[file]<- Sys.glob(paste("/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/featureCounts/", file, sep=""))
}


#
# here we iterate over all input matrixrices and load them into allfiles (gene expression) and allABs (antibody capture)
#
allfiles.raw = list()
for (file in filesasym)
{
  samplename = paste("Asymptomatic_",sub(".*(SRX[0-9]+).*", "\\1", file), sep="")
  foldername = dirname(file)
  
  print(paste(samplename, foldername))

  counts <- read.delim(file, comment.char = "#", row.names = 1)
  
  counts <- counts[ , -(1:5)]
  colnames(counts)=str_extract(colnames(counts), "SRR[0-9]+")
  
  geneinfo=read.table("/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/geneInfo copy.tab", sep = "\t", header = F)
  
  geneinfo$V4 <- ave(geneinfo$V2, geneinfo$V2, FUN = function(x) {
  if(length(x) > 1) {
    paste0(x, ".", seq_along(x))
  } else {
    x
  }
  })
  
  counts=counts[rownames(counts) %in% geneinfo$V1,]
  symbols= geneinfo[geneinfo$V1 %in% rownames(counts),]$V4
  rownames(counts)=symbols

  allfiles.raw[[samplename]] = counts
  
  print(paste(samplename, nrow(allfiles.raw[[samplename]]), "x", ncol(allfiles.raw[[samplename]]), "genes x cells"))
}

for (file in filessym)
{
  samplename = paste("Symptomatic_",sub(".*(SRX[0-9]+).*", "\\1", file), sep="")
  foldername = dirname(file)
  
  print(paste(samplename, foldername))

  counts <- read.delim(file, comment.char = "#", row.names = 1)
  
  counts <- counts[ , -(1:5)]
  colnames(counts)=str_extract(colnames(counts), "SRR[0-9]+")
  
  geneinfo=read.table("/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/geneInfo copy.tab", sep = "\t", header = F)
  
  geneinfo$V4 <- ave(geneinfo$V2, geneinfo$V2, FUN = function(x) {
  if(length(x) > 1) {
    paste0(x, ".", seq_along(x))
  } else {
    x
  }
  })
  
  counts=counts[rownames(counts) %in% geneinfo$V1,]
  symbols= geneinfo[geneinfo$V1 %in% rownames(counts),]$V4
  rownames(counts)=symbols

  allfiles.raw[[samplename]] = counts
  
  print(paste(samplename, nrow(allfiles.raw[[samplename]]), "x", ncol(allfiles.raw[[samplename]]), "genes x cells"))
}

length(allfiles.raw)
names(allfiles.raw)
str(allfiles.raw)


objlist = list()
for (x in names(allfiles.raw))
{
  df_unique <- allfiles.raw[[x]][, !duplicated(colnames(allfiles.raw[[x]]))]
    matrix = df_unique
    
    
    # this creates a Seurat object from the count matrixrix. it sets the object's project to x and prepends the sample name to all cells
    # the patternlist.mouse contains patterns for MT and RP-genes
    filteredObj = makeSeuratObj(matrix, x, patternList.human)
    
    # this creates log-normalized count matrixrices in RNA assay
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
  obj <- subset(obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & nCount_RNA > 10000 & nCount_RNA < 1000000)
  obj <- subset(obj, subset = percent.mt < 15)
  print(obj)
  
  return(obj)
})


for (name in names(objlist))
{
  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  save_plot(p, paste("Seurat/QC",paste(name, "filtered_violins_qc", sep="_"), sep="/"), fig.width=10, fig.height=6)

  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine=F)
  p[[1]] = p[[1]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[2]] = p[[2]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[3]] = p[[3]] + scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,5))
  p = combine_plot_grid_list(plotlist=p, ncol=3)
  save_plot(p, paste("Seurat/QC", paste(name, "filtered_violins_detail_qc", sep="_"), sep="/"), fig.width=18, fig.height=6)
  
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("Seurat/QC", paste(name, "filtered_scatter_ncount_mt", sep="_"), sep="/"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("Seurat/QC", paste(name, "filtered_scatter_ncount_rp", sep="_"), sep="/"), fig.width=10, fig.height=6)

}



###### now we merge all samples into one object
integratedList_plaque <- merge(x = objlist[[1]], y = objlist[2:length(objlist)], project = "human_atherosclerosis")
integratedList_plaque <- JoinLayers(integratedList_plaque)
# 1. Scale data using only variable genes (not all genes)
integratedList_plaque <- ScaleData(integratedList_plaque, features = VariableFeatures(integratedList_plaque))

# 2. Run PCA using variable genes
integratedList_plaque <- RunPCA(integratedList_plaque, features = VariableFeatures(integratedList_plaque))

# 3. Examine PCA elbow plot to choose number of PCs
p=ElbowPlot(integratedList_plaque, ndims = 50)
save_plot(p, "/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/elbow_plot", fig.width = 10, fig.height = 8)


# 4. Use more PCs and try a higher clustering resolution for better separation
integratedList_plaque <- FindNeighbors(integratedList_plaque, dims = 1:40)
integratedList_plaque <- FindClusters(integratedList_plaque, resolution = 0.1)

# 5. Run UMAP for visualization
integratedList_plaque <- RunUMAP(integratedList_plaque, dims = 1:40)

# 6. Assign to obj.in for downstream compatibility
obj.in <- integratedList_plaque

# 7. Visualize clusters
p = DimPlot(obj.in, reduction = "umap", label = TRUE, pt.size = 1.5)
save_plot(p, "/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/Dimplot_umap_ident", fig.width = 10, fig.height = 8)

p=DimPlot(obj.in, reduction = "umap", pt.size=1.5, group.by= "orig.ident")
save_plot(p, "/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/Dimplot_umap_orig.ident", fig.width=10, fig.height=8)

DefaultAssay(obj.in) <- "RNA"
Feat= FeaturePlot(obj.in, features = "MYH11", max.cutoff= 4)
save_plot(Feat, "Seurat/Featureplot_Myh11", fig.width=15, fig.height=10)

DefaultAssay(obj.in) <- "RNA"
Feat= FeaturePlot(obj.in, features = "PECAM1", max.cutoff= 4)
save_plot(Feat, "Seurat/Featureplot_CD31", fig.width=15, fig.height=10)

DefaultAssay(obj.in) <- "RNA"
Feat= FeaturePlot(obj.in, features = "CD3E", max.cutoff= 4)
save_plot(Feat, "Seurat/Featureplot_TCR", fig.width=15, fig.height=10)

Feat= FeaturePlot(obj.in, features = ("ACTA2"))
save_plot(Feat, "Seurat/Featureplot_Acta2", fig.width=15, fig.height=10)

Feat= FeaturePlot(obj.in, features = ("CSPG4"))
save_plot(Feat, "Seurat/Featureplot_NG2", fig.width=15, fig.height=10)

Feat= FeaturePlot(obj.in, features = c("ANXA2", "S100A10"), ncol=2, max.cutoff= 3)
save_plot(Feat, "Seurat/Featureplot_S100a10_ANXA2", fig.width=30, fig.height=10)


obj.in=AddModuleScore(obj.in, features= c("ANXA2", "S100A10"), nbin=1, ctrl= 5, name="Anxa2S100a10_Enriched")
p= FeaturePlot(obj.in,features = "Anxa2S100a10_Enriched1", label = FALSE, repel = TRUE) +
            scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
save_plot(p, "Seurat/Module_Score", fig.width=8, fig.height=5)

###
cellList = colnames(obj.in)

cellList[startsWith(cellList, "Asymptomatic_SRX23809072_")] = "Asymptomatic"
cellList[startsWith(cellList, "Asymptomatic_SRX23809073_")] = "Asymptomatic"
cellList[startsWith(cellList, "Asymptomatic_SRX23809074_")] = "Asymptomatic"
cellList[startsWith(cellList, "Asymptomatic_SRX23809079_")] = "Asymptomatic"
cellList[startsWith(cellList, "Asymptomatic_SRX23809086_")] = "Asymptomatic"
cellList[startsWith(cellList, "Asymptomatic_SRX23809085_")] = "Asymptomatic"
cellList[startsWith(cellList, "Asymptomatic_SRX23809081_")] = "Asymptomatic"
cellList[startsWith(cellList, "Symptomatic_SRX23809083_")] = "Symptomatic"
cellList[startsWith(cellList, "Symptomatic_SRX23809076_")] = "Symptomatic"
cellList[startsWith(cellList, "Symptomatic_SRX23809082_")] = "Symptomatic"
cellList[startsWith(cellList, "Symptomatic_SRX23809075_")] = "Symptomatic"
cellList[startsWith(cellList, "Symptomatic_SRX23809080_")] = "Symptomatic"
cellList[startsWith(cellList, "Symptomatic_SRX23809077_")] = "Symptomatic"
cellList[startsWith(cellList, "Symptomatic_SRX23809078_")] = "Symptomatic"
cellList[startsWith(cellList, "Symptomatic_SRX23809084_")] = "Symptomatic"

obj.in$project=cellList

p=DimPlot(obj.in, label=F, reduction = "umap", group.by= "project", pt.size=1.1, shuffle=T) 
save_plot(p, paste("Seurat", "Dimplot_umap_condition", sep="/"), fig.width=10, fig.height=6)


##############
obj.in = JoinLayers(obj.in)
deResTT = makeDEResults(obj.in, group.by="seurat_clusters", assay="RNA", test="wilcox")
write.table(deResTT,"Seurat/expr_test_clusters_t.tsv", sep="\t", row.names=F, quote = F)
write_xlsx(deResTT, "Seurat/expr_test_clusters_t.xlsx")


#deResTT<-read_tsv("/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/expr_test_clusters_t.tsv")
DefaultAssay(obj.in) <- "RNA"
markers.use.tt= subset(deResTT , avg_log2FC>0&p_val_adj<0.05&!startsWith(gene, "mt-")&!startsWith(gene, "rp"))
finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:20) #subset(clusterID==4)
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



p_dp_genes_idents = DotPlot(obj.in, features = data_dupli$gene, assay="RNA", dot.scale = 5, group.by="CellType")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
    annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin , ymax=ymax, alpha=0.2, fill="blue") #rep(c("blue", "grey"), times= length(events$n)/2)
save_plot(p_dp_genes_idents, "/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/dotplot_cluster_genes_colored", 35, 8)

#write.table(t(data_dupli[data_dupli$clusterID==7,]$gene), file = "cell_markers.tsv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

setwd("/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/")
source("/mnt/raidexttmp/Alejandro/functions.R")
#saveRDS(obj.in, file = "/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")
#obj.in=readRDS("/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat/finalobj.in.rds")



new.cluster.ids <- c("SMCs", "MPs","ECs", "T-cells","PCs")
names(new.cluster.ids) <- levels(obj.in)
obj.in <- RenameIdents(obj.in, new.cluster.ids)
obj.in$CellType <- Idents(obj.in)

p=DimPlot(obj.in, reduction = "umap", label=T, pt.size=1.1, group.by= "CellType", shuffle=F) +
labs(x = "UMAP 1", y = "UMAP 2")
save_plot(p, paste("/mnt/raidexttmp/Alejandro/Smart-seq_Athero_vanlandewijck/Seurat", "Dimplot_umap_label", sep="/"), fig.width=8, fig.height=6)


