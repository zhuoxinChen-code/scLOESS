### Use Seurat for clustering

library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(colorspace)
library(magrittr)
library(reshape2)

### Quality control
cas9.data <- Read10X(data.dir = "filtered_feature_bc_matrix")
cas9_larv <- CreateSeuratObject(counts = cas9.data, project = "larval.fish", min.cells = 3, min.features = 200)

cas9_larv[["percent.mt"]] <- PercentageFeatureSet(cas9_larv, pattern = "^mt-")
VlnPlot(cas9_larv, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Normalization
cas9_larv <- NormalizeData(cas9_larv, normalization.method = "LogNormalize", scale.factor = 10000)

### Feature selection
cas9_larv <- FindVariableFeatures(cas9_larv, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cas9_larv), 10)
#plot1 <- VariableFeaturePlot(cas9_larv)
#plot2 <- LabelPoints(plot = plot1, points = top10)
#CombinePlots(plots = list(plot1, plot2))

### Scaling the data
all.genes <- rownames(cas9_larv)
cas9_larv <- ScaleData(cas9_larv, features = all.genes)
#cas9_larv[["RNA"]]@scale.data

### Perform linear dimensional reduction
cas9_larv <- RunPCA(cas9_larv, features = VariableFeatures(object = cas9_larv))

### Determine the dimensionality
cas9_larv <- JackStraw(cas9_larv, num.replicate = 100)
cas9_larv <- ScoreJackStraw(cas9_larv, dims = 1:20)
#JackStrawPlot(cas9_larv, dims = 1:20)   
#ElbowPlot(cas9_larv,50)   

### Clustering
cas9_larv <- FindNeighbors(cas9_larv, dims = 1:40)
cas9_larv <- FindClusters(cas9_larv, resolution = 1.8)

### Run non-linear dimensional reduction
cas9_larv <- RunTSNE(cas9_larv, dims = 1:40)
DimPlot(cas9_larv, reduction = "tsne")

### Find marker genes
cas9_larv.markers <- FindAllMarkers(cas9_larv, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
#top50 <- cas9_larv.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
#write.table(top50,file="top_differential_genes_50pc.csv")

### tSNE 
tsne <- cas9_larv@reductions$tsne@cell.embeddings %>% as.data.frame() %>% rownames_to_column()
exp.tsne <- Idents(cas9_larv) %>% as.data.frame() %>% rownames_to_column() %>% set_colnames(c("cell", "cluster")) %>%
  left_join(tsne, by = c("cell" = "rowname")) %>% mutate_if(is.factor, as.character) %>% as_tibble()
