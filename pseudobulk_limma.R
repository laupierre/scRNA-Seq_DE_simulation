library(splatter)
library(SCRIP)
library (Seurat)
library (BPCells)

data(acinar.data)
raw_counts <- acinar.data
dim (raw_counts)
## 1000 genes x 80 cells


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



##### Pseudobulk construction






