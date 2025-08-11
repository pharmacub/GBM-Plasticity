---
output: html_document
editor_options: 
  chunk_output_type: console
---

## Metamodules

With this script we built a score from a signature of genes and looked at its expression in cells as if it were a single gene.

```r

# Load Libraries

suppressPackageStartupMessages(library(Seurat))

suppressPackageStartupMessages(library(SeuratData))

suppressPackageStartupMessages(library(patchwork))

suppressPackageStartupMessages(library(dplyr))

suppressPackageStartupMessages(library(ggplot2))

suppressPackageStartupMessages(library(RColorBrewer))

# Load Data

# Files located in (".../final_integ.rds")

cds = readRDS(".../final_integ.rds")

```



```r
# First load a file with the different signatures in columns
metamodules<-read.delim("metamodules.txt",h=T)

# Then  AddMoudleScore function in Seurat to include each of them in a Seurat object. In this script is called <cds>
library(stringi) # load this to remove empty spaces from a txt file with sigature genes

for (i in 1:ncol(metamodules)){
  mod1 <- as.character(stri_remove_empty(metamodules[,i]))
  cds <- AddModuleScore(cds,features = list(mod1), ctrl = length(mod1), name = colnames(metamodules)[i])
}

metamodule.names <- paste0(colnames(metamodules),"1")

# Here there is all the metamodules

pdf("Metamodules_FeaturePlot_1_1.pdf")
FeaturePlot(cds,features = metamodule.names) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
dev.off()
```


```r


VlnPlot(cds,features = metamodule.names,pt.size = 0)
DimPlot(cds,label = T)

```

# DcType and more resolution

To assign an Id to the different clusters. To do it, I increase resolution and then usd scType to assign automatically the scores to single clusters. If the same Id is given to different clusters, they are all merged into one single Id

```r
# cds <- FindNeighbors(cds, dims = 1:30)
cds <- FindClusters(cds, resolution = 2) # res=2 is high enough.

DimPlot(cds,label=T)
# I look at the second module with this new resolution
FeaturePlot(cds,features = metamodule.names[2], pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(cds,features = metamodule.names,pt.size = 0)

```

This is the result of increasing resolution. I do this to increase the granularity that will be used by scType to assigne groups of cells to cell types based on the different modules.

```r
# load libraries and functions

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), function (x) suppressPackageStartupMessages(library(x,character.only = T,quietly = T,)))

# Then I run the two modules from the internet
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Brain") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
```

I have loaded the 'Brain' database. Now I calculate the score. This is just to try the Brain cells signatures
 also needed to extract the scaled data to adjust to SeuratV5 object:

```r
# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(cds[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(cds[["RNA"]]$scale.data) else as.matrix(cds[["RNA"]]@scale.data)

# run ScType on the extracted scaled data
# it works on positive markers (gs_lists$gs_positive).

es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


```

and I prepared all the Ids in the object.

```r
# the clusters are stored in $seurat_clusters.

cL_resutls = do.call("rbind", lapply(unique(cds@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(cds@meta.data[cds@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(cds@meta.data$seurat_clusters==cl)), 10)
}))

# Then it gets the scores
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
# We can set this threshold to avoid having too many cells with unknown. Here they decide if it is less than 1/4 of the cells
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# Here the scores
print(sctype_scores[,1:3])
# Then it sets an Ident, called 'customclassif' to be filled with the Ids

cds@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  cds@meta.data$customclassif[cds@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Now a barplot

par(mar=c(14,6,3,3))
barplot(table(cds$customclassif),las=2)
as.data.frame(table(cds$customclassif))
DimPlot(cds, group.by = "customclassif")
```



```r
cds@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  cds@meta.data$sctype_classification[cds@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(cds, reduction = "umap", label = TRUE, group.by = 'sctype_classification')        

```

# Classification with metamodules

Now I include the metamodules in the ScTypeDB_full.xlsx excel file. I do it by naming the cell modules with 'metamodules' instead of a particular tissue, as 'Brain' in the example before.
And proceed as above:

```r
gs_list = gene_sets_prepare("ScTypeDB_full.xlsx", "metamodules")
# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, #gs2 = gs_list$gs_negative
                       )

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(cds@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(cds@meta.data[cds@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(cds@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

cds@meta.data$metamodules = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  cds@meta.data$metamodules[cds@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

par(mar=c(14,6,3,3))

barplot(table(cds$metamodules),las=2)
as.data.frame(table(cds$metamodules))
```

#ScType Classification Modules

I use the scType classification with the Metamodules 

```r
cds@meta.data$sctype_metamodules = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  cds@meta.data$sctype_metamodules[cds@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(cds, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_metamodules')   

DimPlot(cds,split.by = "RNA_snn_res.0.5",ncol = 3,group.by = "metamodules",label=T)
```


```r
library(loupeR)
create_loupe_from_seurat(cds,output_name="integrated_meta")
saveRDS(cds, file = "integrated_meta.rds")
save.image(file = "integrated.Rdata")
```

