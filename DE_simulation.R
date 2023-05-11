### Simulate single-cell RNA-Seq reads. See https://github.com/thecailab/SCRIP-study/blob/main/R/DE/DE.R

library(splatter)
library(SCRIP)
library (Seurat)
library (BPCells)

data(acinar.data)
raw_counts <- acinar.data
dim (raw_counts)
## 1000 genes x 80 cells


## Load a presaved version of 5000 sketched cells from the human prefrontal cortex
#obj <- readRDS ("seurat3.RDS")
#DefaultAssay(obj) <- "sketch"

## Sketch assay clustering on the 5000 sketched cells
#obj <- FindVariableFeatures(obj)
#obj <- ScaleData(obj)
#obj <- RunPCA(obj)
#obj <- FindNeighbors(obj, dims = 1:50)
#obj <- FindClusters(obj)
#obj <- RunUMAP(obj, dims = 1:50, return.model = T)

# UMAP representation
#DimPlot(obj, label = T, label.size = 5, reduction = "umap") + NoLegend()

## select cluster 0
#obj <- subset (obj, subset = seurat_clusters == 0)

#raw_counts <- obj[["sketch"]]$data
#raw_counts <- raw_counts [ ,1:100]
#raw_counts <- raw_counts[rev (order (apply (raw_counts, 1, sum))) , ]
#raw_counts <- raw_counts[1:12000, ]
#raw_counts <- data.matrix (raw_counts)
#dim (raw_counts)
# 12000   100


params <- splatEstimate (raw_counts)

a=getParams(params, c("mean.rate", "mean.shape"))
rate=a[[1]]
shape=a[[2]]

# install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library("TruncatedDistributions")

ngene <- dim (raw_counts)[1]  # total number of genes
nDE <- ngene/10   # number of DEG. Here we have 10% of total genes

set.seed(2020)
base_allcellmeans <- rtgamma (ngene, shape=shape, scale=1/rate, a=1, b=3)
DEgene <- sample(1:length(base_allcellmeans), nDE, replace=F)


fold <- 2  # linear fold changes
base_allcellmeansDE <- base_allcellmeans
base_allcellmeansDE[DEgene[1:(nDE/2)]] <- base_allcellmeansDE[DEgene[1:(nDE/2)]] *fold
base_allcellmeansDE[DEgene[(nDE/2+1):nDE]] <-base_allcellmeansDE[DEgene[(nDE/2+1):nDE]] *(1/fold)

n.cells <- 100   
drop.out <- 0
batchCells <- n.cells*3  # number of cells in total obtained from 3 mice. ie 100 cells x 3 replicates of mice x 1 condition = 60


# first batch (i.e control group)  
set.seed(2020)
sim <- SCRIPsimu(data=raw_counts, params=params, batchCells=  batchCells, method="single", mode = "GP-trendedBCV", libsize=NULL, bcv.shrink=1,
                 base_allcellmeans_SC= base_allcellmeans, Dropout_rate=drop.out) 

exps <- counts(sim)

# second group (i.e treated group) with some DEG introduced
set.seed(2020)
sim.dif <- SCRIPsimu(data=raw_counts, params=params, batchCells=  batchCells, method="single", mode = "GP-trendedBCV",libsize=NULL, bcv.shrink=1,
                 base_allcellmeans_SC= base_allcellmeansDE, Dropout_rate=drop.out) 

exps.dif <- counts(sim.dif)                     

# cbind the two groups (100 + 100 cells)
counts <- cbind(exps, exps.dif)
# for dropouts
counts[is.na (counts)] <- 0
colnames(counts) <- paste0("cell",1:ncol(counts))
rownames(counts) <- paste0("gene",1:nrow(counts))
# rownames(counts) [DEgene] <- paste (rownames(counts)[DEgene], "-DE", sep="")
rownames(counts) [DEgene[1:(nDE/2)]] <- paste (rownames(counts)[DEgene[1:(nDE/2)]], "-DE-increased", sep="")
rownames(counts) [DEgene[(nDE/2+1):nDE]] <- paste (rownames(counts)[DEgene[(nDE/2+1):nDE]], "-DE-decreased", sep="")



#######################
#### DEG at single cell

#### DESeq2 - single cell

library(DESeq2)

condition <- c (rep ("A", batchCells), rep ("B", batchCells) )
meta <- data.frame (sample= colnames (counts), condition= condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res <- res[order (res$padj), ]
res <- res[res$padj <= 0.05, ]

## precision (how many are real among the positives? == how good it is)
table (grepl ("DE", row.names (res)))[[2]] / dim (res)[1]
## sensitivity (how many are retrieved among the positives? == what we didn't loose)
table (grepl ("DE", row.names (res)))[[2]] / nDE

# verify the log fold change orientation
res.inc <- res[res$log2FoldChange > 0 & res$padj < 0.05, ]
length (grep ("increased", row.names (res.inc), value=TRUE)) / length (row.names (res.inc))

res.dec <- res[res$log2FoldChange < 0 & res$padj < 0.05, ]
length (grep ("decreased", row.names (res.dec), value=TRUE)) / length (row.names (res.dec))




#### Voom-limma - single cell

library (edgeR)
library (limma)

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
  
mm <- model.matrix(~0 + condition)
y <- voom(d0, mm, plot = F)
fit <- lmFit(y, mm)

contr <- makeContrasts(conditionB - conditionA , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res <- topTable(tmp, sort.by = "p", n = Inf) 
res <- res[res$adj.P.Val <= 0.05, ]

## precision (how many are real among the positives? == or how good it is)
table (grepl ("DE", row.names (res)))[[2]] / dim (res)[1]
## sensitivity (how many are retrieved among the initial positive set? == or what we didn't lose)
table (grepl ("DE", row.names (res)))[[2]] / nDE

# verify the log fold change orientation
res.inc <- res[res$logFC > 0 & res$adj.P.Val < 0.05, ]
length (grep ("increased", row.names (res.inc), value=TRUE)) / length (row.names (res.inc))

res.dec <- res[res$logFC < 0 & res$adj.P.Val < 0.05, ]
length (grep ("decreased", row.names (res.dec), value=TRUE)) / length (row.names (res.dec))





#### Using the Libra wrapper for Seurat object at the pseudobulk level. See https://github.com/neurorestore/Libra

library(dplyr)
library(ggplot2)
library(ggrepel)
library(Libra)

condition <- c (rep ("A", batchCells), rep ("B", batchCells) )

## replicate is the mouse (or tissue) where the single cells are collected from (to get the pseudobulk of this mouse)
## cell_type is the cell type to be analyzed
## label is the contrast (i.e treatment)

seurat <- CreateSeuratObject(as (counts, Class = "dgCMatrix") )
seurat@meta.data$cell_type <- "mycell_type"
seurat@meta.data$replicate <- c(rep ("mouse1", n.cells), rep("mouse2", n.cells), rep ("mouse3", n.cells), rep ("mouse1", n.cells), rep("mouse2", n.cells), rep ("mouse3", n.cells))
seurat@meta.data$label <- condition


## wilcoxon - single cell

res <- run_de(seurat, de_method = 'wilcox', de_family= "singlecell")
res <- data.frame (res)
res <- res[res$p_val_adj <= 0.05, ]
row.names (res) <- res$gene

## precision (how many are real among the positives? == or how good it is)
table (grepl ("DE", row.names (res)))[[2]] / dim (res)[1]
## sensitivity (how many are retrieved among the initial positive set? == or what we didn't lose)
table (grepl ("DE", row.names (res)))[[2]] / nDE


## verify the log fold change orientation
## here, it is A-B (and not B-A) !!

res.inc <- res[res$avg_logFC < 0 & res$p_val_adj < 0.05, ]
length (grep ("increased", row.names (res.inc), value=TRUE)) / length (row.names (res.inc))

res.dec <- res[res$avg_logFC > 0 & res$p_val_adj < 0.05, ]
length (grep ("decreased", row.names (res.dec), value=TRUE)) / length (row.names (res.dec))




## Likehood ratio test - single cell

res <- run_de(seurat, de_method = 'bimod', de_family= "singlecell")
res <- data.frame (res)
res <- res[res$p_val_adj <= 0.05, ]
row.names (res) <- res$gene

## precision (how many are real among the positives? == or how good it is)
table (grepl ("DE", row.names (res)))[[2]] / dim (res)[1]
## sensitivity (how many are retrieved among the initial positive set? == or what we didn't lose)
table (grepl ("DE", row.names (res)))[[2]] / nDE



## MAST - single cell

res <- run_de(seurat, de_method = 'MAST', de_family= "singlecell")
res <- data.frame (res)
res <- res[res$p_val_adj <= 0.05, ]
row.names (res) <- res$gene

## precision (how many are real among the positives? == or how good it is)
table (grepl ("DE", row.names (res)))[[2]] / dim (res)[1]
## sensitivity (how many are retrieved among the initial positive set? == or what we didn't lose)
table (grepl ("DE", row.names (res)))[[2]] / nDE



## DESeq2 pseudobulk

res <- run_de(seurat, de_method = 'DESeq2', de_type = 'Wald', de_family= "pseudobulk")
res <- data.frame (res)
res <- res[res$p_val_adj <= 0.05, ]
row.names (res) <- res$gene

## precision (how many are real among the positives? == or how good it is)
table (grepl ("DE", row.names (res)))[[2]] / dim (res)[1]
## sensitivity (how many are retrieved among the initial positive set? == or what we didn't lose)
table (grepl ("DE", row.names (res)))[[2]] / nDE


## Voom-limma pseudobulk

res <- run_de(seurat, de_method = 'limma', de_type = 'voom', de_family= "pseudobulk")
res <- data.frame (res)
res <- res[res$p_val_adj <= 0.05, ]
row.names (res) <- res$gene

## precision (how many are real among the positives? == or how good it is)
table (grepl ("DE", row.names (res)))[[2]] / dim (res)[1]
## sensitivity (how many are retrieved among the initial positive set? == or what we didn't lose)
table (grepl ("DE", row.names (res)))[[2]] / nDE

























