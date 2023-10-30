library(ExperimentHub)
library(zellkonverter)
library(DoubletFinder)
library(tidyverse)
library(openxlsx)
library(SingleR)
library(BPCells)
library(Seurat)
library(readxl)
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
  write_matrix_dir(raw_feature_matrix, dir = glue("./temp/{sample}"))
  raw_feature_matrix <- open_matrix_dir(dir = glue("./temp/{sample}"))
  object <- CreateSeuratObject(raw_feature_matrix, project = sample)
  object@meta.data$group <- group
  
  return(object)

}

FormatGenes <- function(input_string) {
  words <- unlist(strsplit(input_string, ', '))
  transformed_words <- lapply(words, function(word) {
    if (nchar(word) > 1) {
      paste(toupper(substr(word, 1, 1)), tolower(substr(word, 2, nchar(word))), sep = "")
    } else {
      word
    }
  })
  result_vector <- unlist(transformed_words)
  return(result_vector)
}

PlotGeneExpression <- function(data=aggregated_data, ref=canonical_markers_travaglini, cell_type) {
  
  markers <- FormatGenes(canonical_markers_travaglini$markers[,canonical_markers_travaglini$cell_type==cell_type])
  FeaturePlot(aggregated_data, features = markers)
  
}

samples = c('Decell_1m_M1', 'Decell_1m_M2', 'Decell_1m_M3', 'Decell_2wk_M1', 'Decell_2wk_M2', 'Native_U_M1', 'Native_U_M2', 'Native_U_M3', 'Syn_1m_M1', 'Syn_1m_M2', 'Syn_2wk_M1', 'Syn_2wk_M2')

for (sample in samples) {
  assign(sample, LoadAndTreatSamples(sample))
}


#Sample aggregation

aggregated_data <- merge(Decell_1m_M1, y = c(Decell_1m_M2, Decell_1m_M3, Decell_2wk_M1, Decell_2wk_M2, Native_U_M1, Native_U_M2, Native_U_M3, Syn_1m_M1, Syn_1m_M2, Syn_2wk_M1, Syn_2wk_M2), add.cell.ids = c('Decell_1m_M1', 'Decell_1m_M2', 'Decell_1m_M3', 'Decell_2wk_M1', 'Decell_2wk_M2', 'Native_U_M1', 'Native_U_M2', 'Native_U_M3', 'Syn_1m_M1', 'Syn_1m_M2', 'Syn_2wk_M1', 'Syn_2wk_M2'), project = "aggregated_data")
rm(list = samples)


#Cell clustering

aggregated_data <- RunPCA(aggregated_data)
pct <- aggregated_data[["pca"]]@stdev / sum(aggregated_data[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pc_min_1 <- which(cumu > 90 & pct < 5)[1]
pc_min_2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
pcs <- min(pc_min_1, pc_min_2)

aggregated_data <- RunUMAP(aggregated_data, dims = 1:pcs)
aggregated_data <- FindNeighbors(aggregated_data, dims = 1:pcs)
aggregated_data <- FindClusters(aggregated_data, dims = 1:pcs)

VizDimLoadings(aggregated_data, dims = 1:2, reduction = "pca")
DimHeatmap(aggregated_data, dims = 1:15, cells = 5000, balanced = TRUE)
ElbowPlot(aggregated_data, ndims = pcs+5)
DimPlot(aggregated_data, reduction = "pca")
DimPlot(aggregated_data, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(aggregated_data, reduction = "umap", split.by = 'group')


#Cluster marker identification

aggregated_data <- PrepSCTFindMarkers(object = aggregated_data)
all_markers <- FindAllMarkers(aggregated_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(all.markers, 'markers.xlsx')


#SingleR Annotation from Experiment Hub Tabula Muris data

eh <- ExperimentHub()
tm_mm_trachea_ref <- eh[['EH1617']]
tm_mm_trachea_ref <- tm_mm_trachea_ref[,tm_mm_trachea_ref$tissue == 'Trachea']
tm_mm_trachea_ref <- tm_mm_trachea_ref[,!is.na(tm_mm_trachea_ref$cell_ontology_class)]
tm_mm_trachea_ref <- as.Seurat(tm_mm_trachea_ref, data = NULL)
tm_mm_trachea_ref@meta.data$percent.mt <- PercentageFeatureSet(tm_mm_trachea_ref, pattern = "^mt-")
tm_mm_trachea_ref <- SCTransform(tm_mm_trachea_ref, assay = "originalexp", method = "glmGamPoi", min_cells = 3, vars.to.regress = 'percent.mt')
tm_mm_trachea_ref <- as.SingleCellExperiment(tm_mm_trachea_ref)

tm_labels <- SingleR(test = as.SingleCellExperiment(aggregated_data), ref = tm_mm_trachea_ref, labels = tm_mm_trachea_ref@colData$cell_ontology_class, num.threads = 26)
aggregated_data@meta.data$tm_cell_type <- tm_labels$labels
DimPlot(aggregated_data, reduction = "umap", group.by = 'tm_cell_type', label = TRUE)

setwd("D:/Data/IC Ari")


#SingleR Annotation from CellXGene Tabula Muris data

cxg_mm_trachea_ref <- readH5AD('Trachea_droplet.h5ad')
cxg_mm_trachea_ref@assays@data@listData[["counts"]] <- cxg_mm_trachea_ref@assays@data@listData[["X"]]
cxg_mm_trachea_ref <- as.Seurat(cxg_mm_trachea_ref, data = NULL)
cxg_mm_trachea_ref@meta.data$percent.mt <- PercentageFeatureSet(cxg_mm_trachea_ref, pattern = "^mt-")
cxg_mm_trachea_ref <- SCTransform(cxg_mm_trachea_ref, assay = "originalexp", method = "glmGamPoi", min_cells = 3, vars.to.regress = 'percent.mt')
cxg_mm_trachea_ref <- as.SingleCellExperiment(cxg_mm_trachea_ref)

cxg_labels <- SingleR(test = as.SingleCellExperiment(aggregated_data), ref = cxg_mm_trachea_ref, labels = cxg_mm_trachea_ref@colData$cell_ontology_class, num.threads = 26)
aggregated_data@meta.data$cxg_cell_type <- cxg_labels$labels
DimPlot(aggregated_data, reduction = "umap", group.by = 'cxg_cell_type', label = TRUE)


#SingleR Annotation from Human Cell Atlas Lung data

hca_hs_lung_ref <- ntiss10x.P1.anno@data
hca_hs_lung_ref <- convert_matrix_type(hca_hs_lung_ref, type = "uint32_t")
write_matrix_dir(hca_hs_lung_ref, dir = glue("./temp/hca"))
hca_hs_lung_ref <- open_matrix_dir(dir = glue("./temp/hca"))
hca_hs_lung_ref <- CreateSeuratObject(hca_hs_lung_ref, project = 'hca')
hca_hs_lung_ref@meta.data <- ntiss10x.P1.anno@meta.data
hca_hs_lung_ref <- SCTransform(hca_hs_lung_ref, method = "glmGamPoi", min_cells = 3)
hca_hs_lung_ref <- as.SingleCellExperiment(hca_hs_lung_ref)

hca_labels <- SingleR(test = as.SingleCellExperiment(aggregated_data), ref = hca_hs_lung_ref, labels = hca_hs_lung_ref@colData$free_annotation, num.threads = 26)
aggregated_data@meta.data$hca_cell_type <- hca_labels$labels
DimPlot(aggregated_data, reduction = "umap", group.by = 'hca_cell_type', label = TRUE)

#Plotting of cell type canonical markers distribution

canonical_markers_travaglini <- read_excel('canonical_markers_travaglini.xlsx')
canonical_markers_travaglini <- canonical_markers_travaglini[, c(1, 4)]
colnames(canonical_markers_travaglini) <- c('cell_type', 'markers')
canonical_markers_travaglini <- subset(canonical_markers_travaglini, !cell_type %in% c("Capillary Cell", "Eosinophil", "Nonclassical Monocyte"))
canonical_markers_travaglini <- canonical_markers_travaglini[complete.cases(canonical_markers_travaglini), ]

for (ct in canonical_markers_travaglini$cell_type) {
  print(ct)
  p <- PlotGeneExpression(cell_type=ct)
  ggsave(glue('ct_ge/{ct}.png'), p)
}

#FeaturePlot(aggregated_data, features = "Krt5")
#RidgePlot(aggregated_data, features = "Krt5")
