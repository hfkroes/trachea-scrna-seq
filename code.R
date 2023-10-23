library(Seurat)
library(BPCells)
library(DoubletFinder)
library(glue)

options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")

setwd("//wsl.localhost/Ubuntu/home/krois/projects/gags")

FormatObject <- function(object) {
  for (old_name in names(object@meta.data)) {
    if (grepl("^DF", old_name)) {
      new_name <- 'doubletfinder_category'
      eval(parse(text = sprintf('object@meta.data$%s <- object@meta.data$%s', new_name, old_name)))
    }
  }
  return(object)
}

LoadAndTreatSamples <- function(sample) {
  
  group <- substr(sample, 1, nchar(sample) - 3)
  
  options(Seurat.object.assay.version = "v3")
    
  object <- Read10X(data.dir = glue('./{sample}/outs/filtered_feature_bc_matrix'))
  object <- CreateSeuratObject(object)
  object@meta.data$percent.mt <- PercentageFeatureSet(object, pattern = "^mt-")
  object <- subset(object, subset = nFeature_RNA > 200 & percent.mt < 5)
  object <- SCTransform(object, method = "glmGamPoi", min_cells = 3, vars.to.regress = 'percent.mt')
  object <- RunPCA(object)
  object <- RunUMAP(object, dims = 1:10)
  sweep.res.list <- paramSweep_v3(object, PCs = 1:10, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimized_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic.prop <- modelHomotypic(object@meta.data$ClusteringResults)  
  nExp_poi <- round(0.075*nrow(object@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  object <- doubletFinder_v3(object, PCs = 1:10, pN = 0.25, pK = optimized_pK, nExp = nExp_poi.adj, sct = TRUE)
  object <- FormatObject(object)
  fobject <- subset(object, subset = doubletfinder_category == "Singlet")
  raw_st <- fobject@assays$RNA@counts
  
  options(Seurat.object.assay.version = "v5")
  
  raw_feature_matrix <- convert_matrix_type(raw_st, type = "uint32_t")
  write_matrix_dir(raw_feature_matrix, dir = glue("./temp/{sample})"))
  raw_feature_matrix <- open_matrix_dir(dir = glue("./temp/{sample})"))
  object <- CreateSeuratObject(raw_feature_matrix, project = sample)
  object@meta.data$group <- group
  
  return(object)

}

samples = c('Decell_1m_M1', 'Decell_1m_M2', 'Decell_1m_M3', 'Decell_2wk_M1', 'Decell_2wk_M2', 'Native_U_M1', 'Native_U_M2', 'Native_U_M3', 'Syn_1m_M1', 'Syn_1m_M2', 'Syn_2wk_M1', 'Syn_2wk_M2')

for (sample in samples) {
  assign(sample, LoadAndTreatSamples(sample))
}

"""
#Load all 10X data
Decell_1m_M1.data <- Read10X(data.dir = "./Decell_1m_M1/outs/filtered_feature_bc_matrix")
Decell_1m_M1 <- CreateSeuratObject(counts = Decell_1m_M1.data, project = "Decell_1m_M1")
Decell_1m_M1@meta.data$group <- 'Decell_1m'
rm(Decell_1m_M1.data)
Decell_1m_M1[["percent.mt"]] <- PercentageFeatureSet(Decell_1m_M1, pattern = "^mt-")
Decell_1m_M1 <- subset(Decell_1m_M1, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Decell_1m_M1@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Decell_1m_M2.data <- Read10X(data.dir = "./Decell_1m_M2/outs/filtered_feature_bc_matrix")
Decell_1m_M2 <- CreateSeuratObject(counts = Decell_1m_M2.data, project = "Decell_1m_M2")
Decell_1m_M2@meta.data$group <- 'Decell_1m'
rm(Decell_1m_M2.data)
Decell_1m_M2[["percent.mt"]] <- PercentageFeatureSet(Decell_1m_M2, pattern = "^mt-")
Decell_1m_M2 <- subset(Decell_1m_M2, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Decell_1m_M2@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Decell_1m_M3.data <- Read10X(data.dir = "./Decell_1m_M3/outs/filtered_feature_bc_matrix")
Decell_1m_M3 <- CreateSeuratObject(counts = Decell_1m_M3.data, project = "Decell_1m_M3")
Decell_1m_M3@meta.data$group <- 'Decell_1m'
rm(Decell_1m_M3.data)
Decell_1m_M3[["percent.mt"]] <- PercentageFeatureSet(Decell_1m_M3, pattern = "^mt-")
Decell_1m_M3 <- subset(Decell_1m_M3, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Decell_1m_M3@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Decell_2wk_M1.data <- Read10X(data.dir = "./Decell_2wk_M1/outs/filtered_feature_bc_matrix")
Decell_2wk_M1 <- CreateSeuratObject(counts = Decell_2wk_M1.data, project = "Decell_2wk_M1")
Decell_2wk_M1@meta.data$group <- 'Decell_2wk'
rm(Decell_2wk_M1.data)
Decell_2wk_M1[["percent.mt"]] <- PercentageFeatureSet(Decell_2wk_M1, pattern = "^mt-")
Decell_2wk_M1 <- subset(Decell_2wk_M1, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Decell_2wk_M1@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Decell_2wk_M2.data <- Read10X(data.dir = "./Decell_2wk_M2/outs/filtered_feature_bc_matrix")
Decell_2wk_M2 <- CreateSeuratObject(counts = Decell_2wk_M2.data, project = "Decell_2wk_M2")
Decell_2wk_M2@meta.data$group <- 'Decell_2wk'
rm(Decell_2wk_M2.data)
Decell_2wk_M2[["percent.mt"]] <- PercentageFeatureSet(Decell_2wk_M2, pattern = "^mt-")
Decell_2wk_M2 <- subset(Decell_2wk_M2, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Decell_2wk_M2@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Native_U_M1.data <- Read10X(data.dir = "./Native_U_M1/outs/filtered_feature_bc_matrix")
Native_U_M1 <- CreateSeuratObject(counts = Native_U_M1.data, project = "Native_U_M1")
Native_U_M1@meta.data$group <- 'Native_U'
rm(Native_U_M1.data)
Native_U_M1[["percent.mt"]] <- PercentageFeatureSet(Native_U_M1, pattern = "^mt-")
Native_U_M1 <- subset(Native_U_M1, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Native_U_M1@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Native_U_M2.data <- Read10X(data.dir = "./Native_U_M2/outs/filtered_feature_bc_matrix")
Native_U_M2 <- CreateSeuratObject(counts = Native_U_M2.data, project = "Native_U_M2")
Native_U_M2@meta.data$group <- 'Native_U'
rm(Native_U_M2.data)
Native_U_M2[["percent.mt"]] <- PercentageFeatureSet(Native_U_M2, pattern = "^mt-")
Native_U_M2 <- subset(Native_U_M2, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Native_U_M2@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Native_U_M3.data <- Read10X(data.dir = "./Native_U_M3/outs/filtered_feature_bc_matrix")
Native_U_M3 <- CreateSeuratObject(counts = Native_U_M3.data, project = "Native_U_M3")
Native_U_M3@meta.data$group <- 'Native_U'
rm(Native_U_M3.data)
Native_U_M3[["percent.mt"]] <- PercentageFeatureSet(Native_U_M3, pattern = "^mt-")
Native_U_M3 <- subset(Native_U_M3, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Native_U_M3@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Syn_1m_M1.data <- Read10X(data.dir = "./Syn_1m_M1/outs/filtered_feature_bc_matrix")
Syn_1m_M1 <- CreateSeuratObject(counts = Syn_1m_M1.data, project = "Syn_1m_M1")
Syn_1m_M1@meta.data$group <- 'Syn_1m'
rm(Syn_1m_M1.data)
Syn_1m_M1[["percent.mt"]] <- PercentageFeatureSet(Syn_1m_M1, pattern = "^mt-")
Syn_1m_M1 <- subset(Syn_1m_M1, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Syn_1m_M1@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Syn_1m_M2.data <- Read10X(data.dir = "./Syn_1m_M2/outs/filtered_feature_bc_matrix")
Syn_1m_M2 <- CreateSeuratObject(counts = Syn_1m_M2.data, project = "Syn_1m_M2")
Syn_1m_M2@meta.data$group <- 'Syn_1m'
rm(Syn_1m_M2.data)
Syn_1m_M2[["percent.mt"]] <- PercentageFeatureSet(Syn_1m_M2, pattern = "^mt-")
Syn_1m_M2 <- subset(Syn_1m_M2, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Syn_1m_M2@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Syn_2wk_M1.data <- Read10X(data.dir = "./Syn_2wk_M1/outs/filtered_feature_bc_matrix")
Syn_2wk_M1 <- CreateSeuratObject(counts = Syn_2wk_M1.data, project = "Syn_2wk_M1")
Syn_2wk_M1@meta.data$group <- 'Syn_2wk'
rm(Syn_2wk_M1.data)
Syn_2wk_M1[["percent.mt"]] <- PercentageFeatureSet(Syn_2wk_M1, pattern = "^mt-")
Syn_2wk_M1 <- subset(Syn_2wk_M1, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Syn_2wk_M1@meta.data$nCount_RNA, 0.93) & percent.mt < 5)

Syn_2wk_M2.data <- Read10X(data.dir = "./Syn_2wk_M2/outs/filtered_feature_bc_matrix")
Syn_2wk_M2 <- CreateSeuratObject(counts = Syn_2wk_M2.data, project = "Syn_2wk_M2")
Syn_2wk_M2@meta.data$group <- 'Syn_2wk'
rm(Syn_2wk_M2.data)
Syn_2wk_M2[["percent.mt"]] <- PercentageFeatureSet(Syn_2wk_M2, pattern = "^mt-")
Syn_2wk_M2 <- subset(Syn_2wk_M2, subset = nFeature_RNA > 200 & nCount_RNA < quantile(Syn_2wk_M2@meta.data$nCount_RNA, 0.93) & percent.mt < 5)
"""

aggregated_data <- merge(Decell_1m_M1, y = c(Decell_1m_M2, Decell_1m_M3, Decell_2wk_M1, Decell_2wk_M2, Native_U_M1, Native_U_M2, Native_U_M3, Syn_1m_M1, Syn_1m_M2, Syn_2wk_M1, Syn_2wk_M2), add.cell.ids = c('Decell_1m_M1', 'Decell_1m_M2', 'Decell_1m_M3', 'Decell_2wk_M1', 'Decell_2wk_M2', 'Native_U_M1', 'Native_U_M2', 'Native_U_M3', 'Syn_1m_M1', 'Syn_1m_M2', 'Syn_2wk_M1', 'Syn_2wk_M2'), project = "aggregated_data")
rm(list = c('Decell_1m_M1', 'Decell_1m_M2', 'Decell_1m_M3', 'Decell_2wk_M1', 'Decell_2wk_M2', 'Native_U_M1', 'Native_U_M2', 'Native_U_M3', 'Syn_1m_M1', 'Syn_1m_M2', 'Syn_2wk_M1', 'Syn_2wk_M2'))

aggregated_data <- SCTransform(aggregated_data, method = "glmGamPoi", vars.to.regress = "percent.mt", min_cells = 3)
aggregated_data <- RunPCA(aggregated_data)
pct <- aggregated_data[["pca"]]@stdev / sum(aggregated_data[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pc_min_1 <- which(cumu > 90 & pct < 5)[1]
pc_min_2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
pcs <- min(pc_min_1, pc_min_2)
aggregated_data <- RunUMAP(aggregated_data, dims = 1:pcs)
aggregated_data <- FindNeighbors(aggregated_data, dims = 1:pcs)
aggregated_data <- FindClusters(aggregated_data)

VizDimLoadings(aggregated_data, dims = 1:2, reduction = "pca")
DimHeatmap(aggregated_data, dims = 1:15, cells = 5000, balanced = TRUE)
ElbowPlot(aggregated_data, ndims = pcs+5)
DimPlot(aggregated_data, reduction = "pca")
DimPlot(aggregated_data, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(aggregated_data, reduction = "umap", split.by = 'group')

aggregated_data <- PrepSCTFindMarkers(object = aggregated_data)
all.markers <- FindAllMarkers(aggregated_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
