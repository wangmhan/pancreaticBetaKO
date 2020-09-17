setwd("/Users/menwan/Documents/project/genewizPanc");
library(readr); library(Seurat); 
library(dplyr); library(scater)

read_data <- "panc_scReadCount_afterQC.txt"

expdata <- read.table(gzfile(read_data), header = TRUE, sep="\t", row.names=1)
dim(expdata); expdata <- as.matrix(expdata)
expdata <- expdata[,order(colnames(expdata))]

cell_type <- substr(colnames(expdata),1,5);
names(cell_type) <- colnames(expdata)

###########################################
#          input into Seurat  	          #
###########################################
panc <- CreateSeuratObject(raw.data = expdata, min.cells = 2, #filter genes
  project = "pancBetaCell")
VlnPlot(object = panc, features.plot = c("nGene", "nUMI"), group.by="orig.ident")

panc@meta.data$orig.ident <- substr(panc@meta.data$orig.ident,1,5)
panc <- NormalizeData(object = panc, 
  normalization.method = "LogNormalize", scale.factor = 10000)

#For example for UMI data normalized to a total of 10 000 molecules, 
#one would expect ~2,000 variable genes.
panc <- FindVariableGenes(object = panc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 0.25)
length(x = panc@var.genes)
panc <- ScaleData(object = panc, vars.to.regress = c("nUMI", "percent.mito"))

tmp1 <- factor(substr(panc@ident,1,5)); names(tmp1) <- names(panc@ident)

panc <- RunPCA(object = panc, pc.genes = panc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)

PCAPlot(object = panc, dim.1 = 1, dim.2 = 2, group.by="orig.ident")
PCHeatmap(object = panc, pc.use = 1:3, do.balanced = TRUE, label.columns = FALSE)

###########################################
#      canonical scRNA-seq pipeline       #
###########################################
panc <- JackStraw(object = panc, num.replicate = 1000, display.progress = FALSE)
JackStrawPlot(object = panc, PCs = 1:12)
PCElbowPlot(object = panc) #identified the significant PCs, determining which PCs 
panc <- FindClusters(object = panc, reduction.type = "pca", dims.use = 1:12, 
    resolution = 0.8, print.output = 0, save.SNN = TRUE) 
table(panc@ident)
panc <- RunTSNE(object = panc, dims.use = 1:12, do.fast = TRUE, max_iter = 2000)
TSNEPlot(object = panc, do.label = TRUE)

