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

# number of cells
n.cells <- 100   
drop.out <- 0
batchCells <- n.cells*3  # number of cells in total obtained from 3 mice. ie 100 cells x 3 replicates of mice x 1 condition = 300 cells


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



##### Pseudobulk construction with Libra (option 1)

library (Libra)

# 100 cells from 3 mice in two conditions (A and B)
label <- c (rep ("A", batchCells), rep ("B", batchCells) )
cell_type <- "mycell_type"
replicate <- c(rep ("mouse1", n.cells), rep("mouse2", n.cells), rep ("mouse3", n.cells), rep ("mouse1", n.cells), rep("mouse2", n.cells), rep ("mouse3", n.cells))

meta <- data.frame (cbind (label, cell_type, replicate))
row.names (meta) <- colnames (counts)
head (meta)

pseudo.counts <- to_pseudobulk(counts, meta = meta)
pseudo.counts <- pseudo.counts$mycell_type
head (pseudo.counts)



#### Voom-limma

library (edgeR)
library (limma)

d0 <- DGEList(pseudo.counts)
d0 <- calcNormFactors(d0)
  
condition <- c (rep ("A", 3), rep ("B", 3))

mm <- model.matrix(~0 + condition)

y <- voom(d0, mm, plot = T)

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




##### Manual Pseudobulk construction  (option 2)

meta$mouse <- paste (meta$replicate, meta$label, sep=":")
mouse <- unique (meta$mouse)

mat <- list ()
for (i in (1:length (mouse))) {
mycells <- counts[ ,colnames (counts) %in% row.names (meta)[meta$mouse == mouse[i] ]]
a <- data.frame (rowSums (mycells))
row.names (a) <- row.names (mycells)
colnames (a)[1] <- mouse[i]
mat[[i]] <- a
}

pseudo.counts.2 <- do.call ("cbind", mat)
identical (pseudo.counts, pseudo.counts.2)
# TRUE




######
### Removing heteroscedasticity with Voom
### See https://github.com/YOU-k/voomByGroup and https://you-k.github.io
### git clone https://github.com/YOU-k/voomByGroup.git

# source ("~/voomByGroup/voomByGroup.R")


## voomQW with sample variability
y <- voomWithQualityWeights(d0, design = mm, plot = TRUE)
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


## voomQW with block variability
y <- voomWithQualityWeights(d0, design = mm, var.group=condition , plot = TRUE)
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


## voomQW with group variability
y <- voomByGroup(d0, design = mm, group=condition , plot = "all")
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








