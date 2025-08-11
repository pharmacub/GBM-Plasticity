In single-cell RNA sequencing (scRNA-seq) experiments, data integration refers to the process of combining and harmonizing gene expression data from multiple samples or datasets into a single integrated dataset. This integration is necessary when dealing with scRNA-seq data from different experimental conditions, biological replicates, or sequencing batches, which may exhibit technical variations or batch effects.

The goal of data integration in scRNA-seq is to find common cell subgroups across both datasets. Determine markers indicative of cell types that remain consistent in both control and stimulated cells. Analyze the datasets to pinpoint how different cell types respond uniquely to stimulation.

First we start by loading the libraries and reading in the data.
**merge()** function merges multiple single-cell RNA-seq datasets into a single object.

```r

# Load Libraries

library(Rcpp)

library(harmony)

suppressPackageStartupMessages(library(Seurat))

suppressPackageStartupMessages(library(SeuratData))

suppressPackageStartupMessages(library(patchwork))

suppressPackageStartupMessages(library(dplyr))



# Load Data

# Files located in (".../Hamed_plasticity/samples_with_ccc/sample_1.rds")

#                  (".../Hamed_plasticity/samples_with_ccc/sample_2.rds")

#                  (".../Hamed_plasticity/samples_with_ccc/sample_3.rds")


cds1 = readRDS(".../Hamed_plasticity/samples_with_ccc/sample_1.rds")
cds1$old.ident <- Idents(cds1)

cds2 = readRDS(".../Hamed_plasticity/samples_with_ccc/sample_2.rds")

cds2$old.ident <- Idents(cds2)

cds3 = readRDS(".../Hamed_plasticity/samples_with_ccc/sample_3.rds")

cds3$old.ident <- Idents(cds3)

# Merge the Objects without Interation

merged = merge(x = cds1, y = list(cds2, cds3))

# Load the data from ccc after merging in RPCA integration method

rpca = readRDS(".../Hamed_plasticity/integration/final/rpca.rds")
```
```r

# Visualization

DimPlot(cds1 , reduction = "umap")

DimPlot(cds1, group.by = c("Phase", "seurat_clusters"), label = T)


```

```r

# Visualization

DimPlot(cds2 , reduction = "umap")

DimPlot(cds2, group.by = c("Phase", "seurat_clusters"))


```

```r

# Visualization

DimPlot(cds3 , reduction = "umap")

DimPlot(cds3, group.by = c("Phase", "seurat_clusters"), label = T)


```


#### **Perform analysis without integration:**

We can first analyze the dataset without integration. 
After loading and merging the data, the standard workflow of Seurat for single-cell RNA sequencing (scRNA-seq) data analysis performed followed by **normalizing** the data, and **highly variable genes are identified**. Principal Component Analysis **(PCA)** is used for dimensionality reduction, and cells are clustered using a Shared Nearest Neighbor **(SNN)** graph and **clustering** algorithms. visualization is performed with **UMAP**.

```r

# Run standard analysis workflow

# Normalization the data

merged <- NormalizeData(merged,verbose = F)

# Identify the most highly variable genes

merged <- FindVariableFeatures(merged,verbose = F)

# Scale the data

merged <- ScaleData(merged, verbose = F)

# linear dimensional reduction (running PCA)

merged <- RunPCA(merged, verbose = F)

```

```r

merged <- FindNeighbors(merged, dims = 1:30, reduction = "pca", verbose = F)

merged <- FindClusters(merged, resolution = 0.5, verbose = F)

merged <- RunUMAP(merged, dims = 1:30, reduction = "pca", verbose = F)

```

The resulting clusters are defined by three different samples, which creates challenges for downstream analysis.

```r

# Visualization 

DimPlot(merged, label = T)

DimPlot(merged, group.by = c("orig.ident", "seurat_clusters"),label=T)

DimPlot(merged, group.by = c("Phase", "seurat_clusters"))

```

#### **Perform integration:**

We now aim to integrate data from the three samples.

##### **Harmony integration method:**

Harmony's approach to **batch correction** is particularly effective because it aligns the data at the cluster level, maintaining the integrity of biological signals. Unlike some methods that apply global corrections, Harmony’s local adjustments within clusters ensure that the biological variations are not overshadowed by batch effects. This makes Harmony a robust and reliable method for integrating scRNA-seq data from multiple sources, enabling comprehensive and accurate downstream analyses.


**Benefits of Harmony:**

**Preservation of Biological Signal:** By focusing on aligning data within clusters, Harmony ensures that the biological signals are preserved.

**Efficiency:** Harmony is computationally efficient, capable of handling large datasets with multiple batches.

**Flexibility:** The method can be applied to various types of omics data and is adaptable to different analysis pipelines.


```r

# Performing integration with Harmony method 

merged_harmony <- RunHarmony(merged, group.by.vars = "orig.ident", theta = 3, verbose = F)

# re-join layers after integration

merged_harmony [["RNA"]] <- JoinLayers(merged_harmony [["RNA"]])


```


```r

merged_harmony <- FindNeighbors(merged_harmony, reduction = "harmony", dims = 1:30, verbose = F)
merged_harmony <- FindClusters(merged_harmony , resolution = 0.5, verbose = F)

merged_harmony <- RunUMAP(merged_harmony , dims = 1:30, reduction = "harmony", verbose = F)

```


```r

# Visualization

DimPlot(merged_harmony, reduction = "umap", group.by = c("orig.ident", "RNA_snn_res.0.5"),label = T)

DimPlot(merged_harmony, group.by = c("Phase", "RNA_snn_res.0.5"))


DimPlot(merged_harmony, reduction = "umap", split.by = c("orig.ident"))

```


```r

# Compare the integrated clusters with the original cluster labels

table_harmony <- table(merged_harmony$seurat_clusters, merged$old.ident, merged$orig.ident)

# Visualize a contingency table
print(table_harmony)
```

###### **Finding differentially expressed features (cluster biomarkers):**

```r

# find all markers distinguishing cluster 5 from cluster 10

cluster5vs10.markers <- FindMarkers(merged_harmony, ident.1 = 10, ident.2 = 5)

head(cluster5vs10.markers, n = 10)

```

```r

v1 = merged_harmony$RNA_snn_res.0.5
v2 = rpca$RNA_snn_res.0.7

table_rpca <- table(v1 , v2)

table(v1,merged_harmony$orig.ident)

tt = apply(table_rpca,2,function(x) x*100/sum(x))
library(pheatmap)

pheatmap(tt,cluster_rows = F,cluster_cols = F)

print(table_rpca)


DimPlot(merged_harmony, reduction = "umap")
DimPlot(rpca, reduction = "umap")
```

```r


dim(merged_harmony)
dim(cluster5vs10.markers)

up_in_cl5 = cluster5vs10.markers[which(cluster5vs10.markers$p_val_adj<0.05 & cluster5vs10.markers$avg_log2FC>3),]

write.table(row.names(up_in_cl5),file="up.txt",sep="\t",quote=F,row.names = F,col.names = F)

down_in_cl5 = cluster5vs10.markers[which(cluster5vs10.markers$p_val_adj<0.05 & cluster5vs10.markers$avg_log2FC < -3),]

write.table(row.names(down_in_cl5),file="down.txt",sep="\t",quote=F,row.names = F,col.names = F)
dim(down_in_cl5)

write.table(row.names(merged_harmony),file="all.txt",sep="\t",quote=F,row.names = F,col.names = F)

```

```r

saveRDS(merged_harmony,file = "final_integ.rds")
save.image(file = "final_integ.Rdata")
```
