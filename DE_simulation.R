### Simulate single-cell RNA-Seq reads. See https://github.com/thecailab/SCRIP-study/blob/main/R/DE/DE.R
### cd /Volumes/king/ad_singlecell

library(splatter)
library(SCRIP)
library (Seurat)
library (BPCells)

#data(acinar.data)
#dim (acinar.data)
## 1000 genes x 80 cells

## Load a presaved version of 5000 sketched cells from the human prefrontal cortex
obj <- readRDS ("seurat3.RDS")
DefaultAssay(obj) <- "sketch"

## Sketch assay clustering on the 5000 sketched cells
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj)
obj <- RunUMAP(obj, dims = 1:50, return.model = T)

# UMAP representation
DimPlot(obj, label = T, label.size = 5, reduction = "umap") + NoLegend()

## select cluster 0
obj <- subset (obj, subset = seurat_clusters == 0)

raw_counts <- obj[["sketch"]]$data
raw_counts <- raw_counts [ ,1:100]
dim (raw_counts)
# 33538   100

params <- splatEstimate(data.matrix (raw_counts))

a=getParams(params, c("mean.rate", "mean.shape"))
rate=a[[1]]
shape=a[[2]]

# install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library("TruncatedDistributions")
set.seed(2020)

ngene <- dim (raw_counts)[1]  # total number of genes
nDE <- ngene/10   # number of DEG. Here we have 10% of total genes

base_allcellmeans <- rtgamma (ngene, shape=shape, scale=1/rate, a=1, b=3)
DEgene <- sample(1:length(base_allcellmeans), nDE, replace=F)


fold <- 2  # linear fold changes
base_allcellmeansDE <- base_allcellmeans
base_allcellmeansDE[DEgene[1:(nDE/2)]] <- base_allcellmeansDE[DEgene[1:(nDE/2)]] *fold
base_allcellmeansDE[DEgene[(nDE/2+1):nDE]] <-base_allcellmeansDE[DEgene[(nDE/2+1):nDE]] *(1/fold)

   
batchCells <- 60  # number of cells in total obtained from 3 mice. ie 20 cells x 3 replicates of mice x 1 condition = 60


# first batch (i.e control group)  
sim <- SCRIPsimu(data=data.matrix (raw_counts), params=params, batchCells=  batchCells, method="single", mode = "GP-trendedBCV", libsize=NULL, bcv.shrink=1,
                 base_allcellmeans_SC= base_allcellmeans, Dropout_rate=0) 

exps <- counts(sim)



# second group (i.e treated group) with some DEG introduced
sim.dif <- SCRIPsimu(data=data.matrix (raw_counts), params=params, batchCells=  batchCells, method="single", mode = "GP-trendedBCV",libsize=NULL, bcv.shrink=1,
                 base_allcellmeans_SC= base_allcellmeansDE, Dropout_rate=0) 

exps.dif <- counts(sim.dif)                     


# cbind the two groups (300 + 300 cells)
counts <- cbind(exps, exps.dif)

colnames(counts) <- paste0("cell",1:ncol(counts))
rownames(counts) <- paste0("gene",1:nrow(counts))
rownames(counts) [DEgene] <- paste (rownames(counts)[DEgene], "-DE", sep="")
#### DESeq2 - single cell


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

## sensitivity (how many are retrieved among the positives? == what we lost)
table (grepl ("DE", row.names (res)))[[2]] / nDE



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




#### Using the Libra wrapper for Seurat object at the pseudobulk level. See https://github.com/neurorestore/Libra

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(Libra)


## replicate is the mouse (or tissue) where the single cells are collected from (to get the pseudobulk of this mouse)
## cell_type is the cell type to be analyzed
## label is the contrast (i.e treatment)

seurat <- CreateSeuratObject(as (counts, Class = "dgCMatrix") )
seurat@meta.data$cell_type <- "mycell_type"
seurat@meta.data$replicate <- c(rep ("mouse1", 100), rep("mouse2", 100), rep ("mouse3", 100), rep ("mouse1", 100), rep("mouse2", 100), rep ("mouse3", 100))
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

























