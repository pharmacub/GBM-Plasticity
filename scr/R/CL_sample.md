---
title: "sample_1"
output: github_document
date: 
editor_options:
  chunk_output_type: console
---


```r

# Load Libraries


suppressPackageStartupMessages(library(Seurat))

suppressPackageStartupMessages(library(dplyr))

suppressPackageStartupMessages(library(patchwork))


# Load Data

# File located in (...)


cds = readRDS("~...")

```

#### **Standard Pre-processing:**

scRNA-seq data sets contain systematic and random noise (such as from poor-quality cells) that obscures the biological signal. Pre-processing of scRNA-seq data attempts to remove these confounding sources of variation. This involves quality control, normalization, data correction and feature selection.

##### **Quality Control:**

In single-cell RNA sequencing (scRNA-seq) data analysis, the percentage of reads that map to the mitochondrial genome is often used as a quality control metric. Mitochondrial contamination is indicative of low-quality or dying cells, as these cells tend to have higher mitochondrial RNA content due to mitochondrial dysfunction or increased cellular stress.

To calculate mitochondrial QC metrics, such as the percentage of counts originating from mitochondrial genes, you can use methods like the **PercentageFeatureSet()** function. In this case, you would use a set of features corresponding to mitochondrial genes. Typically, mitochondrial genes are prefixed with **"MT-"** or are explicitly annotated as mitochondrial genes in the reference genome.

By calculating the percentage of counts originating from mitochondrial genes, you can identify cells with high mitochondrial content, which may indicate low-quality or dying cells. These cells can then be filtered out from downstream analysis to ensure that the remaining cells are of higher quality and are more likely to represent viable cells with meaningful biological information.

```r

# QC mitochondrial

cds[["percent.mt"]] <- PercentageFeatureSet(cds, pattern = "^MT-")

```

Quality Control includes steps for identifying and filtering out low-quality cells and genes. Cells with low RNA content or high levels of ambient RNA can be removed to ensure that downstream analysis is performed on reliable data.

For this step we using **Violin plots** that are used when you want to observe the distribution of numeric data, and are especially useful when you want to make a comparison of distributions between multiple groups.

**nFeature_RNA:** The number of detected features (genes) in each cell.

**nCount_RNA:** The total number of molecules detected within a cell.


```r

# Visualize the QC matrices in Violin plot

VlnPlot(cds, features = c("nFeature_RNA","nCount_RNA","percent.mt") , ncol = 3, pt.size = 0)

```


```r
#visualize feature-feature relationships
plot1 <- FeatureScatter(cds, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cds, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```


**Based on the plots:**

1. Eliminate cells with more than **6000** features (doubles): This threshold likely aims to remove cells that have unusually high numbers of features, which could indicate potential technical artifacts such as doublets or multiplets.

2. Eliminate cells with less than **1500** features: Cells with very low numbers of features could be dead or empty, indicating low-quality data. Removing such cells helps ensure that the remaining cells are more likely to be viable and provide meaningful biological information.

3. Eliminate cells with more than **20%** mitochondrial gene expression: Cells with high percentages of mitochondrial gene expression are often associated with cell lysis after death or mitochondrial stress. Removing these cells helps filter out potentially low-quality or dying cells from the analysis.

By applying these criteria, we can focus the analysis on cells that are likely to be of higher quality and biologically relevant, while excluding those that may introduce noise or artifacts.

```r
cds <- subset(cds, subset = nFeature_RNA > 1500 & nFeature_RNA < 6000 & percent.mt < 20)
```

##### **Normalizing the data:**
The normalization process aims to make gene expression profiles comparable across cells, accounting for differences in sequencing depth and other technical biases.raw data are adjusted to account for factors that prevent direct comparison of expression measures.
Here we using **log normalize** method that Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor that is 10000 in here.


```r

# normalization the data

cds <- NormalizeData(cds, normalization.method = "LogNormalize", scale.factor = 10000)

```

##### **feature selection:**

Identifying genes with high cell-to-cell variation, often referred to as highly variable genes (HVGs), is a crucial step in single-cell RNA sequencing (scRNA-seq) analysis. These genes are informative because they capture biological heterogeneity across individual cells, helping to highlight meaningful biological signals in the dataset.

```r

cds <- FindVariableFeatures(cds, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(cds), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(cds)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

```

##### **Scaling the data:**

After normalization, each gene's expression values are typically log-transformed and scaled to ensure that each gene has a similar variance and mean across cells. This process is essential for downstream analyses such as dimensionality reduction and clustering, where differences in expression levels between genes could otherwise dominate the analysis.

After scaling the data, each gene's expression values are centered around **0** and scaled to have a standard deviation of **1** across cells. This ensures that genes with large expression values do not dominate downstream analyses and that the relative importance of genes is determined by their variation across cells rather than their absolute expression levels.

The **CellCycleScoring()** function calculates cell cycle scores for individual cells in scRNA-seq data based on the expression levels of genes associated with different phases of the cell cycle, such as **S** phase and **G2M** phase.

To avoid the Mitochondrial contamination and Cell-cycle effect using the **vars.to.regress** function needed.

```r

# scale the data

# remove the mitochondrial one

# cds <- ScaleData(cds, vars.to.regress = "percent.mt")

# to avoid the cell cycle effect

all.genes <- rownames(cds)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cds<- CellCycleScoring(cds, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

cds<- ScaleData(cds, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), features = all.genes)
```

#### **linear dimensional reduction:**

**Principal Component Analysis (PCA)** is a technique used to emphasize variation as well as similarity, and to bring out strong patterns in a dataset; it is one of the methods used for “dimensionality reduction”.

```r

#linear dimensional reduction (running PCA)

cds<- RunPCA(cds, features = VariableFeatures(cds))
```

The **elbow plot** visualizes the standard deviation of each PC. Where the elbow appears is usually the threshold for identifying the majority of the variation. However, this method can be a bit subjective about where the elbow is located.

```r
#to find how many PCA we wanted to use 
ElbowPlot(cds, ndims = 50)
```

##### **Cluster the cells:**

To find the nearest neighbors for a given dataset in the context of principal component analysis (PCA), we using **k-nearest neighbors (KNN)** method  based on the principal components (PCs) obtained from PCA.

In scRNA-seq analysis, **KNN** is commonly used to define cell-to-cell similarity or distance metrics based on gene expression profiles. Each cell's nearest neighbors are identified based on these metrics, and cells with similar expression profiles are grouped together into clusters. This approach helps in identifying cell populations and characterizing cell heterogeneity within a dataset.

The **KNN** algorithm itself is relatively simple. To classify a new data point, it calculates the distance between that point and all other points in the dataset. The K nearest neighbors of the new point are then identified based on these distances. For classification tasks, the most common class among these neighbors is assigned to the new point. For regression tasks, the output value is typically the average or weighted average of the values of the K nearest neighbors.

While KNN is a widely used and effective method, there are alternative algorithms for cell clustering and dimensionality reduction in scRNA-seq analysis. Some of these include hierarchical clustering, spectral clustering, and density-based clustering methods such as DBSCAN. Additionally, for dimensionality reduction, algorithms like principal component analysis **(PCA)**, independent component analysis **(ICA)**, and non-negative matrix factorization **(NMF)** are commonly used alternatives to **t-SNE** and **UMAP**. The choice of algorithm depends on the specific characteristics of the dataset, computational efficiency, and the desired output.

```r
cds <- FindNeighbors(cds, dims = 1:40)
cds <- FindClusters(cds, resolution = 0.44)
```

**UMAP or Uniform Manifold Approximation and Projection**, is a dimensionality reduction technique commonly used in single-cell RNA sequencing (scRNA-seq) data analysis. UMAP is used to visualize high-dimensional data in a lower-dimensional space while preserving the underlying structure and relationships between data points.

UMAP has gained popularity in scRNA-seq analysis due to its computational efficiency, scalability to large datasets, and ability to preserve both local and global structure. It aims to find a low-dimensional representation of the data that retains the intrinsic geometry of the high-dimensional space, making it well-suited for visualizing complex biological datasets such as scRNA-seq data.

In scRNA-seq analysis, UMAP is often used to visualize cellular heterogeneity, identify cell populations, and explore the relationships between cells. It complements other dimensionality reduction techniques like PCA, providing additional insights into the structure and organization of the data.

```r
cds <- RunUMAP(cds, dims = 1:40)
```

```r
DimPlot(cds, reduction = "umap")
```


#### **Finding differentially expressed features (cluster biomarkers):**

To figure out what types of **cells** are in a cluster and what makes each cluster unique, we look for genes that are highly expressed in each cluster in compression to the others. These genes, called **marker** genes, act like labels, helping us distinguish one cluster from another.

The **FindAllMarkers()** is identify genes that are differentially expressed between different clusters or groups of cells in an scRNA-seq dataset. By comparing gene expression levels across clusters, the function identifies genes that are significantly upregulated or downregulated in one cluster compared to others. These genes are often referred to as marker genes because they mark or characterize specific cell types, subpopulations, or biological states within the dataset.

**only.pos = TRUE:** This argument specifies that only positively differentially expressed genes should be returned as marker genes. These are genes that are significantly upregulated in a specific cluster compared to all other clusters.

**dplyr::filter(avg_log2FC > 1):** This filters the marker genes based on a threshold for average log2 fold change (avg_log2FC). Here, it selects only those marker genes where the average log2 fold change is greater than 1, indicating that they are significantly upregulated in the corresponding cluster compared to other clusters.

```r

# find markers for every cluster compared to all remaining cells, report only the positive ones


cds.markers <- FindAllMarkers(cds ,only.pos = TRUE, verbose = FALSE, min.pct = 0.1)
cds.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

```


The **heatmap** is a graphical representation of data where values in a matrix are represented as colors. It is particularly useful for visualizing large datasets with multiple variables.

In scRNA-seq analysis, heatmaps are frequently used to explore **gene expression** patterns across individual cells, identify clusters of cells with similar expression profiles, and discover genes that are differentially expressed between cell populations. They provide a visual summary of the underlying structure and heterogeneity within the dataset, aiding in the interpretation of complex biological phenomena.

Here we are plotting the top 10 markers for each cluster.

```r

cds.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(cds, features = top10$gene) 

```


```r

write.table(cds.markers, file = "~.../CL_markers.txt", quote = F, sep="\t" ) 

```


```r
sessionInfo()
```