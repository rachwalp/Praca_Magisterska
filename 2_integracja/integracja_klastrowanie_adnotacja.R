# Biblioteki
library(Seurat)
library(ggplot2)
library(mclust)
library(harmony)
library('glmGamPoi')
library(clustree)
library(reshape2)
library(Azimuth)


# Zaladowanie probek
sample_26ctrl <- readRDS("26ctrl.rds")
sample_153ctri <- readRDS("153ctri.rds")
sample_173ctrl <- readRDS("173ctrl.rds")
sample_180ctrl <- readRDS("180ctrl.rds")
sample_A00103 <- readRDS("A00103.rds")
sample_A03211 <- readRDS("A03211.rds")
sample_A05234 <- readRDS("A05234.rds")
sample_Y06102 <- readRDS("Y06102.rds")
sample_1 <- readRDS("sample1_neg.rds")
sample_2 <- readRDS("sample2.rds")
sample_3 <- readRDS("sample3.rds")
sample_4 <- readRDS("sample4_neg.rds")
sample_95 <- readRDS("sample95_neg.rds")
sample_99 <- readRDS("sample99_neg.rds")
sample_BM102 <- readRDS("BM102.rds")
sample_BM106 <- readRDS("BM106.rds")
sample_BM107 <- readRDS("BM107.rds")
sample_BM108 <- readRDS("BM108.rds")
sample_BM114 <- readRDS("BM114.rds")
sample_BM115 <- readRDS("BM115.rds")
sample_BM119 <- readRDS("BM119.rds")
sample_BM120 <- readRDS("BM120.rds")



# Funkcja pomocnicza do wyliczania statystyk.
print_seurat_stats <- function(seurat_object) {
  if (is.null(seurat_object) || length(seurat_object) == 0) {
    cat("The Seurat object is empty or NULL.\n")
    return(NULL)
  }
  
  remaining_cells <- ncol(seurat_object)
  remaining_genes <- nrow(seurat_object)
  cat("Number of cells:", remaining_cells, "\n")
  cat("Number of genes:", remaining_genes, "\n")
}


seurat_objects <- list(sample_26ctrl, sample_153ctri,sample_173ctrl,sample_180ctrl,sample_A00103,sample_A03211,sample_A05234,sample_Y06102,sample_1, sample_2, sample_3, sample_4,sample_95,sample_99, sample_BM102, sample_BM106, sample_BM107, sample_BM108, sample_BM114, sample_BM115, sample_BM119, sample_BM120)

names(seurat_objects) <- c("26ctrl","153ctri","173ctrl","180ctrl","A00103","A03211","A05234","Y06102","sample1", "sample2", "sample3", "sample4","sample95","sample99", "sample_BM102", "sample_BM106", "sample_BM107", "sample_BM108", "sample_BM114", "sample_BM115", "sample_BM119", "sample_BM120")

for (obj_name in names(seurat_objects)) {
  seurat_object <- seurat_objects[[obj_name]]
  cat("\nStatistics for", obj_name, ":\n")
  print_seurat_stats(seurat_object)
}

# Sprawdzenie duplikatow
cells_sample26ctrl <- colnames(sample_26ctrl)
cells_sample1153ctri <- colnames(sample_153ctri)
cells_sample173ctrl <- colnames(sample_173ctrl)
cells_sample180ctrl <- colnames(sample_180ctrl)
cells_sampleA00103 <- colnames(sample_A00103)
cells_sampleA03211 <- colnames(sample_A03211)
cells_sampleA05234 <- colnames(sample_A05234)
cells_sampleY06102 <- colnames(sample_Y06102)
cells_sample1 <- colnames(sample_1)
cells_sample2 <- colnames(sample_2)
cells_sample3 <- colnames(sample_3)
cells_sample4 <- colnames(sample_4)
cells_sample95 <- colnames(sample_95)
cells_sample99 <- colnames(sample_99)
cells_sampleBM102 <- colnames(sample_BM102)
cells_sampleBM106 <- colnames(sample_BM106)
cells_sampleBM107 <- colnames(sample_BM107)
cells_sampleBM108 <- colnames(sample_BM108)
cells_sampleBM114 <- colnames(sample_BM114)
cells_sampleBM115 <- colnames(sample_BM115)
cells_sampleBM119 <- colnames(sample_BM119)
cells_sampleBM120 <- colnames(sample_BM120)
all_cells <- c(cells_sample1, cells_sample2, cells_sample3, cells_sample4, cells_sample95, cells_sample99,cells_sample26ctrl, cells_sample1153ctri, cells_sample173ctrl, cells_sample180ctrl,cells_sampleA00103,cells_sampleA03211, cells_sampleA05234, cells_sampleY06102, cells_sampleBM102, cells_sampleBM106, cells_sampleBM107, cells_sampleBM108, cells_sampleBM114, cells_sampleBM115, cells_sampleBM119, cells_sampleBM120)

cell_counts <- table(all_cells)

duplicate_cells <- names(cell_counts[cell_counts > 1])

length(duplicate_cells)

# Zmiana nazw komorek w zwiazku z obecnoscia duplikatow 
sample_26ctrl   <- RenameCells(sample_26ctrl, add.cell.id = "26ctrl")
sample_153ctri  <- RenameCells(sample_153ctri, add.cell.id = "153ctri")
sample_173ctrl  <- RenameCells(sample_173ctrl, add.cell.id = "173ctrl")
sample_180ctrl  <- RenameCells(sample_180ctrl, add.cell.id = "180ctrl")
sample_A00103   <- RenameCells(sample_A00103, add.cell.id = "A00103")
sample_A03211   <- RenameCells(sample_A03211, add.cell.id = "A03211")
sample_A05234   <- RenameCells(sample_A05234, add.cell.id = "A05234")
sample_Y06102   <- RenameCells(sample_Y06102, add.cell.id = "Y06102")
sample_1        <- RenameCells(sample_1, add.cell.id = "1")
sample_2        <- RenameCells(sample_2, add.cell.id = "2")
sample_3        <- RenameCells(sample_3, add.cell.id = "3")
sample_4        <- RenameCells(sample_4, add.cell.id = "4")
sample_95       <- RenameCells(sample_95, add.cell.id = "95")
sample_99       <- RenameCells(sample_99, add.cell.id = "99")
sample_BM102    <- RenameCells(sample_BM102, add.cell.id = "BM102")
sample_BM106    <- RenameCells(sample_BM106, add.cell.id = "BM106")
sample_BM107    <- RenameCells(sample_BM107, add.cell.id = "BM107")
sample_BM108    <- RenameCells(sample_BM108, add.cell.id = "BM108")
sample_BM114    <- RenameCells(sample_BM114, add.cell.id = "BM114")
sample_BM115    <- RenameCells(sample_BM115, add.cell.id = "BM115")
sample_BM119    <- RenameCells(sample_BM119, add.cell.id = "BM119")
sample_BM120    <- RenameCells(sample_BM120, add.cell.id = "BM120")

samples_list <- list(
  sample_26ctrl, sample_153ctri, sample_173ctrl, sample_180ctrl,
  sample_A00103, sample_A03211, sample_A05234, sample_Y06102,
  sample_1, sample_2, sample_3, sample_4, sample_95, sample_99,
  sample_BM102, sample_BM106, sample_BM107, sample_BM108,
  sample_BM114, sample_BM115, sample_BM119, sample_BM120
)

# Zlaczenie probek 
merged <- merge(
  x = samples_list[[1]],
  y = samples_list[-1],
  merge.data = FALSE
)

cat("Number of cells in merged object:", ncol(merged), "\n")

merged$sample <- sapply(strsplit(colnames(merged), "_"), `[`, 1)
table(merged$sample)

merged[["percent.mt"]] <- PercentageFeatureSet(
  merged,
  pattern = "^MT-"
)

DefaultAssay(merged)
VlnPlot(merged, features = "nFeature_RNA", group.by = "sample")


# SCTransform
merged <- SCTransform(
  merged,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)


# PCA
merged <- RunPCA(merged)


# Integracja
merged <- RunHarmony(
  object = merged,
  group.by.vars = "sample",
  assay.use = "SCT",      # or "RNA" if not using SCT
  reduction.use = "pca",
  dims.use = 1:30
)


merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:15)
merged <- FindClusters(merged, resolution = 1.2)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:15)

merged <- RunUMAP(merged, reduction = "pca", dims = 1:15, reduction.name = "umap_pca")
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:15, reduction.name = "umap_harmony")

DimPlot(merged, reduction = "umap_pca", group.by = "sample")
DimPlot(merged, reduction = "umap_harmony", group.by = "sample")


# Klastrowanie

startRes = 0.2
endRes = 2.8 
step = 0.2

resolutions <- seq(startRes, endRes, by = step)
for (r in resolutions) {
  merged <- FindClusters(merged, resolution = r)
}

clustree(merged, prefix = "SCT_snn_res.")

ARImatrix=matrix(nrow=length(seq(startRes,endRes,step)),ncol=length(seq(startRes,endRes,step)))
for(i in 1:length(seq(startRes,endRes,step))) {
  for(j in 1:length(seq(startRes,endRes,step))) {
    if(i!=j) {
      vecA=merged@meta.data[,paste("SCT_snn_res.",seq(startRes,endRes,step)[i],sep="")]
      vecB=merged@meta.data[,paste("SCT_snn_res.",seq(startRes,endRes,step)[j],sep="")]
      ARImatrix[i,j]=adjustedRandIndex(vecA,vecB)
    }
  }
}

ARImean=NULL
ARIsd=NULL
ARInmbClusters=NULL
for(i in 1:length(seq(startRes,endRes,step))) {
  ARImean=c(ARImean,mean(ARImatrix[i,-i]))
  ARIsd=c(ARIsd,sd(ARImatrix[i,-i]))
  ARInmbClusters=c(ARInmbClusters,length(unique(merged@meta.data[,paste("SCT_snn_res.",seq(startRes,endRes,step)[i],sep="")])))
}
optimalRes=seq(startRes,endRes,step)[which(ARImean==max(ARImean))[1]]
ARIdf=data.frame(resolution=seq(startRes,endRes,step), ARI_mean=ARImean, ARI_sd=ARIsd, number_clusters=ARInmbClusters)
ARIdf=melt(ARIdf,id='resolution')

# Wykres ARI
print(ggplot(ARIdf, aes(x = resolution, y = value, fill = variable)) + geom_line(aes(group=variable,color=variable), size = 1.2) + labs(y=NULL) + labs(x = "resolution") + ggtitle("Adjusted Rand Index") + facet_grid(variable ~ . , scales = "free") + theme_bw() + geom_vline(xintercept=optimalRes,color="orange",linetype="dotted"))

# Wykres klastrowania 
DimPlot(
  merged,
  group.by = "SCT_snn_res.1.2",
  label = TRUE
)

# Adnotacje
# Adnotacje wykonanu przy urzyciu pakietu Azimuth w terminalu Ubuntu.
# Adnotacja zostala wykonana automatycznie przy pomocy funkcji RunAzimuth()

azimuth_adnotation <-  RunAzimuth(merged, reference = "bonemarrowref")

DimPlot(azimuth_analysis, group.by = "predicted.celltype.l2", label = TRUE)

# Adnotacje wykonane przez Azimuth nie zawieraly komorek nowotworowych dlatego adnotacja do dalszych etapow analiz zostala wykonana recznie.
