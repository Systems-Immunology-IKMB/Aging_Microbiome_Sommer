 ## Summary
  # Samples
  ## Mice were exposed to sequential antibiotic treatment followed by colonization with either microbiome from same-aged mice (isochronic microbiome-iMB) or with the microbiome from 8 week old mice (young microbiome-yMB). At 134 weeks mice 6 mice were sacrificed and their intestinal mucosa prepared for scRNAseq technology.
  
  # Tissue separation
  ## The intestinal epithelium of each sample was separated from the lamina propria using a LP kit separation method. Thus, each sample is comprised of IECs and LP subsets.
library(Seurat)
library(monocle3)
library(dplyr)
library(Matrix)
library(gdata)
library(ggplot2)
library(reshape2)
library(devtools)
library(ggpubr)
library(SingleR)
library(RcisTarget)

library(clusterProfiler)
library(org.Mm.eg.db)
library(plyr)
library(RColorBrewer)
library(cowplot)
library(tidyr)
library(rstatix)
library(ggprism)
library(UpSetR)
library(ComplexHeatmap)
library(ggpmisc)
library(RcisTarget)
library(stringr)
library(rstatix)
set.seed(42)


# 1. Load samples from aging-MB iMB
dirname <- "~/input/"
counts_matrix_filename = paste0(dirname,"/aging-MB_Final_54_LP/raw_feature_bc_matrix/")
counts_matrix_filename2 = paste0(dirname,"/aging-MB_Final_36_LP/raw_feature_bc_matrix/")
counts_matrix_filename3 = paste0(dirname,"/aging-MB_Final_10_LP/raw_feature_bc_matrix/")

counts_matrix_filename4 = paste0(dirname,"/aging-MB_Final_54_IEC/raw_feature_bc_matrix/")
counts_matrix_filename5 = paste0(dirname,"/aging-MB_Final_36_IEC/raw_feature_bc_matrix/")
counts_matrix_filename6 = paste0(dirname,"/aging-MB_Final_10_IEC/raw_feature_bc_matrix/")

Sample1 <- Read10X(data.dir = counts_matrix_filename) 
Sample2 <- Read10X(data.dir = counts_matrix_filename2) 
Sample3 <- Read10X(data.dir = counts_matrix_filename3) 
Sample4 <- Read10X(data.dir = counts_matrix_filename4) 
Sample5 <- Read10X(data.dir = counts_matrix_filename5) 
Sample6 <- Read10X(data.dir = counts_matrix_filename6) 

Sample_aging_54_LP<-CreateSeuratObject(counts = Sample1, min.cells = 3, min.features = 200, project = "Sample_aging_54_LP")
Sample_aging_36_LP<-CreateSeuratObject(counts = Sample2, min.cells = 3, min.features = 200, project = "Sample_aging_36_LP")
Sample_aging_10_LP<-CreateSeuratObject(counts = Sample3, min.cells = 3, min.features = 200, project = "Sample_aging_10_LP")
Sample_aging_54_IEC<-CreateSeuratObject(counts = Sample4, min.cells = 3, min.features = 200, project = "Sample_aging_54_IEC")
Sample_aging_36_IEC<-CreateSeuratObject(counts = Sample5, min.cells = 3, min.features = 200, project = "Sample_aging_36_IEC")
Sample_aging_10_IEC<-CreateSeuratObject(counts = Sample6, min.cells = 3, min.features = 200, project = "Sample_aging_10_IEC")

#mitochodrial
Sample_aging_54_LP[["percent.mt"]] <- PercentageFeatureSet(Sample_aging_54_LP, pattern = "^mt-")
Sample_aging_36_LP[["percent.mt"]] <- PercentageFeatureSet(Sample_aging_36_LP, pattern = "^mt-")
Sample_aging_10_LP[["percent.mt"]] <- PercentageFeatureSet(Sample_aging_10_LP, pattern = "^mt-")
Sample_aging_54_IEC[["percent.mt"]] <- PercentageFeatureSet(Sample_aging_54_IEC, pattern = "^mt-")
Sample_aging_36_IEC[["percent.mt"]] <- PercentageFeatureSet(Sample_aging_36_IEC, pattern = "^mt-")
Sample_aging_10_IEC[["percent.mt"]] <- PercentageFeatureSet(Sample_aging_10_IEC, pattern = "^mt-")

# Load the the list of house keEping genes (tirosh paper)
Genes<-read.csv('~/tirosh_house_keeping.csv')
hkgenes <- as.vector(Genes$Mouse)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_aging_54_LP))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_aging_54_LP[hkgenes.found, ] )> 0)
Sample_aging_54_LP[['percentage.hk']] <- PercentageFeatureSet(Sample_aging_54_LP, features = n.expressed.hkgenes)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_aging_36_LP))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_aging_36_LP[hkgenes.found, ] )> 0)
Sample_aging_36_LP[['percentage.hk']] <- PercentageFeatureSet(Sample_aging_36_LP, features = n.expressed.hkgenes)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_aging_10_LP))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_aging_10_LP[hkgenes.found, ] )> 0)
Sample_aging_10_LP[['percentage.hk']] <- PercentageFeatureSet(Sample_aging_10_LP, features = n.expressed.hkgenes)

hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_aging_54_IEC))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_aging_54_IEC[hkgenes.found, ] )> 0)
Sample_aging_54_IEC[['percentage.hk']] <- PercentageFeatureSet(Sample_aging_54_IEC, features = n.expressed.hkgenes)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_aging_36_IEC))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_aging_36_IEC[hkgenes.found, ] )> 0)
Sample_aging_36_IEC[['percentage.hk']] <- PercentageFeatureSet(Sample_aging_36_IEC, features = n.expressed.hkgenes)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_aging_10_IEC))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_aging_10_IEC[hkgenes.found, ] )> 0)
Sample_aging_10_IEC[['percentage.hk']] <- PercentageFeatureSet(Sample_aging_10_IEC, features = n.expressed.hkgenes)

# Raw data visualization
pdf("Sample_aging_iMB_LP_raw.pdf", width = 10, height = 5)
VlnPlot(object = Sample_aging_54_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_aging_36_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_aging_10_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
dev.off()

pdf("Sample_aging_iMB_IEC_raw.pdf", width = 10, height = 5)
VlnPlot(object = Sample_aging_54_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_aging_36_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_aging_10_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
dev.off()


## Filter
Sample_aging_54_LP <- subset(Sample_aging_54_LP, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)
Sample_aging_36_LP <- subset(Sample_aging_36_LP, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)
Sample_aging_10_LP <- subset(Sample_aging_10_LP, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)

Sample_aging_54_IEC <- subset(Sample_aging_54_IEC, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)
Sample_aging_36_IEC<- subset(Sample_aging_36_IEC, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)
Sample_aging_10_IEC <- subset(Sample_aging_10_IEC, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)


# Filter data visualization
pdf("Sample_aging_iMB_LP_filter.pdf", width = 10, height = 5)
VlnPlot(object = Ageing_V3_clean,features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_aging_36_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
VlnPlot(object = Sample_aging_10_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
dev.off()

pdf("Sample_aging_iMB_IEC_filter.pdf", width = 20, height = 10)
VlnPlot(object = Sample_aging_54_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
VlnPlot(object = Sample_aging_36_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
VlnPlot(object = Sample_aging_10_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
dev.off()

Sample_aging_54_LP$Mouse='M54'
Sample_aging_36_LP$Mouse='M36'
Sample_aging_10_LP$Mouse='M10'
Sample_aging_54_IEC$Mouse='M54'
Sample_aging_36_IEC$Mouse='M36'
Sample_aging_10_IEC$Mouse='M10'

Sample_aging_54_LP$Tissue='LP'
Sample_aging_36_LP$Tissue='LP'
Sample_aging_10_LP$Tissue='LP'
Sample_aging_54_IEC$Tissue='IEC'
Sample_aging_36_IEC$Tissue='IEC'
Sample_aging_10_IEC$Tissue='IEC'

Sample_aging_54_LP$Treatment='iMB'
Sample_aging_36_LP$Treatment='iMB'
Sample_aging_10_LP$Treatment='iMB'
Sample_aging_54_IEC$Treatment='iMB'
Sample_aging_36_IEC$Treatment='iMB'
Sample_aging_10_IEC$Treatment='iMB'

save(Sample_aging_54_LP, Sample_aging_36_LP, Sample_aging_10_LP, Sample_aging_54_IEC, Sample_aging_36_IEC, Sample_aging_10_IEC,file ='iMB_IndSamples.rda' )


# 2. Load samples from aging-MB yMB
dirname <- "~/input/"
counts_matrix_filename = paste0(dirname,"/8w_Final_16_LP/raw_feature_bc_matrix/")
counts_matrix_filename2 = paste0(dirname,"/8w_Final_3_LP/raw_feature_bc_matrix/")
counts_matrix_filename3 = paste0(dirname,"/8w_Final_2_LP/raw_feature_bc_matrix/")

counts_matrix_filename4 = paste0(dirname,"/8w_Final_16_IEC/raw_feature_bc_matrix/")
counts_matrix_filename5 = paste0(dirname,"/8w_Final_3_IEC/raw_feature_bc_matrix/")
counts_matrix_filename6 = paste0(dirname,"/8w_Final_2_IEC/raw_feature_bc_matrix/")

Sample1 <- Read10X(data.dir = counts_matrix_filename) 
Sample2 <- Read10X(data.dir = counts_matrix_filename2) 
Sample3 <- Read10X(data.dir = counts_matrix_filename3) 
Sample4 <- Read10X(data.dir = counts_matrix_filename4) 
Sample5 <- Read10X(data.dir = counts_matrix_filename5) 
Sample6 <- Read10X(data.dir = counts_matrix_filename6) 

Sample_8w_16_LP<-CreateSeuratObject(counts = Sample1, min.cells = 3, min.features = 200, project = "Sample_8w_16_LP")
Sample_8w_3_LP<-CreateSeuratObject(counts = Sample2, min.cells = 3, min.features = 200, project = "Sample_8w_3_LP")
Sample_8w_2_LP<-CreateSeuratObject(counts = Sample3, min.cells = 3, min.features = 200, project = "Sample_8w_2_LP")
Sample_8w_16_IEC<-CreateSeuratObject(counts = Sample4, min.cells = 3, min.features = 200, project = "Sample_8w_16_IEC")
Sample_8w_3_IEC<-CreateSeuratObject(counts = Sample5, min.cells = 3, min.features = 200, project = "Sample_8w_3_IEC")
Sample_8w_2_IEC<-CreateSeuratObject(counts = Sample6, min.cells = 3, min.features = 200, project = "Sample_8w_2_IEC")

#mitochodrial
Sample_8w_16_LP[["percent.mt"]] <- PercentageFeatureSet(Sample_8w_16_LP, pattern = "^mt-")
Sample_8w_3_LP[["percent.mt"]] <- PercentageFeatureSet(Sample_8w_3_LP, pattern = "^mt-")
Sample_8w_2_LP[["percent.mt"]] <- PercentageFeatureSet(Sample_8w_2_LP, pattern = "^mt-")
Sample_8w_16_IEC[["percent.mt"]] <- PercentageFeatureSet(Sample_8w_16_IEC, pattern = "^mt-")
Sample_8w_3_IEC[["percent.mt"]] <- PercentageFeatureSet(Sample_8w_3_IEC, pattern = "^mt-")
Sample_8w_2_IEC[["percent.mt"]] <- PercentageFeatureSet(Sample_8w_2_IEC, pattern = "^mt-")

# Load the the list of house keEping genes (tirosh paper)
Genes<-read.csv('~/tirosh_house_keeping.csv')
hkgenes <- as.vector(Genes$Mouse)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_8w_16_LP))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_8w_16_LP[hkgenes.found, ] )> 0)
Sample_8w_16_LP[['percentage.hk']] <- PercentageFeatureSet(Sample_8w_16_LP, features = n.expressed.hkgenes)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_8w_3_LP))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_8w_3_LP[hkgenes.found, ] )> 0)
Sample_8w_3_LP[['percentage.hk']] <- PercentageFeatureSet(Sample_8w_3_LP, features = n.expressed.hkgenes)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_8w_2_LP))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_8w_2_LP[hkgenes.found, ] )> 0)
Sample_8w_2_LP[['percentage.hk']] <- PercentageFeatureSet(Sample_8w_2_LP, features = n.expressed.hkgenes)

hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_8w_16_IEC))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_8w_16_IEC[hkgenes.found, ] )> 0)
Sample_8w_16_IEC[['percentage.hk']] <- PercentageFeatureSet(Sample_8w_16_IEC, features = n.expressed.hkgenes)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_8w_3_IEC))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_8w_3_IEC[hkgenes.found, ] )> 0)
Sample_8w_3_IEC[['percentage.hk']] <- PercentageFeatureSet(Sample_8w_3_IEC, features = n.expressed.hkgenes)
hkgenes.found <- which(toupper(rownames(GetAssayData(object = Sample_8w_2_IEC))) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(GetAssayData(object = Sample_8w_2_IEC[hkgenes.found, ] )> 0)
Sample_8w_2_IEC[['percentage.hk']] <- PercentageFeatureSet(Sample_8w_2_IEC, features = n.expressed.hkgenes)

# Raw data visualization
pdf("Sample_aging_iMB_LP_raw.pdf", width = 10, height = 5)
VlnPlot(object = Sample_8w_16_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_8w_3_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_8w_2_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
dev.off()

pdf("Sample_aging_iMB_IEC_raw.pdf", width = 10, height = 5)
VlnPlot(object = Sample_8w_16_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_8w_3_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_8w_2_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
dev.off()


## Filter
Sample_8w_16_LP <- subset(Sample_8w_16_LP, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)
Sample_8w_3_LP <- subset(Sample_8w_3_LP, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)
Sample_8w_2_LP <- subset(Sample_8w_2_LP, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)

Sample_8w_16_IEC <- subset(Sample_8w_16_IEC, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)
Sample_8w_3_IEC<- subset(Sample_8w_3_IEC, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)
Sample_8w_2_IEC <- subset(Sample_8w_2_IEC, subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & percent.mt < 25 & percentage.hk > 10)


# Filter data
pdf("Sample_aging_ visualizationiMB_LP_filter.pdf", width = 10, height = 5)
VlnPlot(object = Ageing_V3_clean,features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0)
VlnPlot(object = Sample_8w_3_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
VlnPlot(object = Sample_8w_2_LP, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
dev.off()

pdf("Sample_aging_iMB_IEC_filter.pdf", width = 20, height = 10)
VlnPlot(object = Sample_8w_16_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
VlnPlot(object = Sample_8w_3_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
VlnPlot(object = Sample_8w_2_IEC, features = c( "nFeature_RNA", "nCount_RNA","percent.mt", 'percentage.hk'), ncol =  4, pt.size = 0.1)
dev.off()

Sample_8w_16_LP$Mouse='M16'
Sample_8w_3_LP$Mouse='M3'
Sample_8w_2_LP$Mouse='M2'
Sample_8w_16_IEC$Mouse='M16'
Sample_8w_3_IEC$Mouse='M3'
Sample_8w_2_IEC$Mouse='M2'

Sample_8w_16_LP$Tissue='LP'
Sample_8w_3_LP$Tissue='LP'
Sample_8w_2_LP$Tissue='LP'
Sample_8w_16_IEC$Tissue='IEC'
Sample_8w_3_IEC$Tissue='IEC'
Sample_8w_2_IEC$Tissue='IEC'

Sample_8w_16_LP$Treatment='yMB'
Sample_8w_3_LP$Treatment='yMB'
Sample_8w_2_LP$Treatment='yMB'
Sample_8w_16_IEC$Treatment='yMB'
Sample_8w_3_IEC$Treatment='yMB'
Sample_8w_2_IEC$Treatment='yMB'

save(Sample_8w_16_LP, Sample_8w_3_LP, Sample_8w_2_LP, Sample_8w_16_IEC, Sample_8w_3_IEC, Sample_8w_2_IEC,file ='yMB_IndSamples.rda' )


# Combine Treatments correcting for batch (Mice)
load('~/Documents/Documents - Joana’s MacBook Pro/IKMB/scRNAseq/Ageing/Manuscript/Analysis_July_2022/yMB_IndSamples.rda')
load('~/Documents/Documents - Joana’s MacBook Pro/IKMB/scRNAseq/Ageing/Manuscript/Analysis_July_2022/iMB_IndSamples.rda')
# load('~/Documents/IKMB/SC RNA Seq/Ageing/Manuscript/control_combine.rda')

# rename cells to avoid identical cell barcodes
Sample_aging_54_LP <- RenameCells(object = Sample_aging_54_LP, add.cell.id = "S54_LP")
Sample_aging_36_LP <- RenameCells(object = Sample_aging_36_LP, add.cell.id = "S36_LP")
Sample_aging_10_LP <- RenameCells(object = Sample_aging_10_LP, add.cell.id = "S10_LP")
Sample_aging_54_IEC <- RenameCells(object = Sample_aging_54_IEC, add.cell.id = "S54_IEC")
Sample_aging_36_IEC <- RenameCells(object = Sample_aging_36_IEC, add.cell.id = "S36_IEC")
Sample_aging_10_IEC <- RenameCells(object = Sample_aging_10_IEC, add.cell.id = "S10_IEC")

Sample_8w_16_LP <- RenameCells(object = Sample_8w_16_LP, add.cell.id = "S16_LP")
Sample_8w_3_LP <- RenameCells(object = Sample_8w_3_LP, add.cell.id = "S3_LP")
Sample_8w_2_LP <- RenameCells(object = Sample_8w_2_LP, add.cell.id = "S2_LP")
Sample_8w_16_IEC <- RenameCells(object = Sample_8w_16_IEC, add.cell.id = "S16_IEC")
Sample_8w_3_IEC <- RenameCells(object = Sample_8w_3_IEC, add.cell.id = "S3_IEC")
Sample_8w_2_IEC <- RenameCells(object = Sample_8w_2_IEC, add.cell.id = "S2_IEC")

# simple merge sample objects into one object
Treatment_pre<-merge(x = Sample_aging_54_LP, y = c(Sample_aging_54_IEC, Sample_aging_36_LP, Sample_aging_36_IEC, Sample_aging_10_LP, Sample_aging_10_IEC,
                                                   Sample_8w_16_LP,Sample_8w_16_IEC, Sample_8w_3_LP, Sample_8w_3_IEC,Sample_8w_2_LP,Sample_8w_2_IEC))

# integration of samples with batch correction for individual mice
Treatment_sep<-SplitObject(Treatment_pre, split.by = 'Mouse')
for(i in 1:length(Treatment_sep)) {
  Treatment_sep[[i]] <- NormalizeData(Treatment_sep[[i]], verbose = F)
  Treatment_sep[[i]] <- FindVariableFeatures(Treatment_sep[[i]], selection.method = 'vst', verbose = F)
}
ref_Treatment<-Treatment_sep[c('M2','M3','M16', 'M10','M36','M54')]
Sample_V3_anchors<-FindIntegrationAnchors(object.list = ref_Treatment, dims = 1:60)
Aging_V3<-IntegrateData(anchorset = Sample_V3_anchors, dims=1:60)

# save merged object
save(Aging_V3, file ='Aging-MB_V3_merged.rda' )

# analysis and separation into compartments
# clustering workflow analysis
all.genes<-rownames(Aging_V3)
Aging_V3 <- ScaleData(Aging_V3, features = all.genes)
Aging_V3 <- FindVariableFeatures(Aging_V3)
Aging_V3 <- RunPCA(Aging_V3,npcs = 80)
ElbowPlot(Aging_V3, ndims = 80)

Aging_V3 <- RunUMAP(object = Aging_V3, dims = 1:80)
Aging_V3 <- FindNeighbors(Aging_V3,  dims = 1:80)
Aging_V3 <- FindClusters(Aging_V3, resolution = 0.2)

# Visualization of merged object
pdf("UMAP_Aging-MB_V3_Sample_Treatment_Mouse.pdf", width = 8, height = 7)
DimPlot(Aging_V3, reduction = "umap",group.by = 'orig.ident', label=F, pt.size = 0.1)
DimPlot(Aging_V3, reduction = "umap",group.by = 'Treatment', label=F, pt.size = 0.1, order = c('yMB', 'iMB'),
        pt.size = 0.1, cols = c("#0072B2", "#D55E00"))
DimPlot(Aging_V3, reduction = "umap",group.by = 'Mouse', label=F, pt.size = 0.1)
dev.off()

pdf("UMAP_Aging-MB_V3_clusters.pdf", width = 13, height = 10)
DimPlot(Aging_V3, reduction = "umap",group.by = 'seurat_clusters', label=TRUE)
dev.off()

# Biomarkers for compartments (Epithelial, Stromal and Immune) in object clusters (res=0.2)
Epi_genes<-c('Epcam', 'Krt8', 'Krt18')
Stromal_genes<-c('Col1a1', 'Col1a2', 'Col6a1','Col6a2', 'Vwf', 'Plvap', 'Cdh5', 'S100b')
Immune_genes<-c('Cd52', 'Cd2', 'Cd3d','Cd3g', 'Cd3e', 'Cd79a', 'Cd79b', 'Cd14', 'Fcgr3','Cd68', 'Cd83', 'Csf1r', 'Fcer1g')
DotPlot(Aging_V3, features = c(Epi_genes,Stromal_genes,Immune_genes),dot.scale = 7, assay='RNA')+ RotatedAxis()

Idents(Aging_V3)<- Aging_V3$seurat_clusters
new.cluster.ids <- c("Epithelium", 'Epithelium','Epithelium','Epithelium','Immune', 'Epithelium',
                     'Immune', 'Stromal','Immune','Immune', 'Epithelium',
                     'Stromal', 'Immune','Immune','Immune','Stromal',
                     'Epithelium','Epithelium', 'Epithelium', 'Immune', 'Stromal',
                     'Epithelium', 'Immune','Stromal')
names(new.cluster.ids) <- levels(Aging_V3)
Aging_V3 <- RenameIdents(Aging_V3, new.cluster.ids)
Aging_V3$Compartment<-Idents(Aging_V3)

# Visualization of merged object compartments
pdf("UMAP_Aging-MB_V3_Compartment.pdf", width = 8, height = 7)
DimPlot(Aging_V3, reduction = "umap",group.by = 'Compartment', label=F,
        cols = c("bisque3","coral4",   "coral2"))
dev.off()

# SingleR cell-type calling
mouse.se<-MouseRNAseqData()
input <- as.matrix(GetAssayData(object = Aging_V3, slot = "data", assay = "RNA"))
singleR.list <- list()

# perform singleR classification
singleR.list$mouse <- SingleR(test = input, 
                              method="single",
                              fine.tune=FALSE,
                              ref = mouse.se, 
                              labels = mouse.se$label.main)

rm(input)
Ageing_V3_clean$SingleR <- singleR.list$mouse$labels

# Visualization of merged object SingleR cell-types
pdf("UMAP_Aging-MB_V3_SingleR.pdf", width = 8, height = 6)
DimPlot(object = Aging_V3, reduction = 'umap',  group.by ="SingleR", pt.size = 0.1, label = TRUE)
dev.off()

# Reference cell-type calling base on Tabula muris (2018) facs dataset. To note, this dataset only comprises large intestine samples so there is a mismatch with our object which is comprised of small intestinal epithelium cells
load("~/TM_Intestine_facs_Figure.rds")

reference<-SCTransform(TM_Int,  verbose = FALSE)
query <- SCTransform(Aging_V3,  verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals=FALSE,
  dims = 1:50
)

query <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(
    celltype.l1 = "free_annotation",
    celltype.l2 = "cell_ontology_class",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

# add cell-type metadata from reference
Aging_V3$celltype.l2 <- query$predicted.celltype.l2
Aging_V3$celltype.l1 <- query$predicted.celltype.l1

# visualization of cell-type annotaations base on Tabula Muris reference
pdf("UMAP_Aging-MB_V3_Celltype_Ref2_free_annotation.pdf", width = 12, height = 6)
DimPlot(Aging_V3, reduction = "umap",group.by = 'celltype.l2', label=TRUE)
dev.off()

pdf("UMAP_Aging-MB_V3_Celltype_Ref1_cell_ontology_class.pdf", width = 12, height = 6)
DimPlot(Aging_V3, reduction = "umap",group.by = 'celltype.l1', label=TRUE)
dev.off()

# Biomarkers for cell-types base on literature (Tabula Muris) and in-house markers
Idents(Aging_V3)<-Aging_V3$seurat_clusters
Stem<-c("Lgr5","Olfm4","Ascl2", "Slc12a1", "Gkn3")
Enterocytes<-c("Alpi","Fabp1","Apoa1", "Apoa4")
Goblet<-c("Muc2","Agr2","Clca3a1", "Tff3", "Clca3a2")
Paneth<-c("Lyz1","Defa17","Defa22", "Defa24", "Ang4")
EE<- c("Chga","Chgb","Tac1", "Tph1")
Tuft<-c("Dclk1","Trpm5","Gfi1b", "Il25")
TA<-c("Stmn1","Tubb5")
PE<-c('Sox4', 'Neurog3', 'Neurod2')
DotPlot(Aging_V3, features = c(Stem,Enterocytes,Goblet, Paneth, EE, Tuft, TA, PE),dot.scale = 7,assay='RNA')+ RotatedAxis()

EnteroMatProx<-c("Apoa4","Fabp1","Apoc2", "Rbp2", "Apoc3", 'Leap2')
EnteroImatProx<-c("Casp6")
EnteroMatDistal<- c("Tmigd1","Fabp6","Slc51a", "Slc51b", "Mep1a", 'Fam151a')
EnteroImatDistal<-c("Reg3g","Gsdmc4","Prss32", "Krt8")
DotPlot(Aging_V3, features = c(EnteroMatProx,EnteroImatProx,EnteroMatDistal, EnteroImatDistal),dot.scale = 7,assay='RNA')+ RotatedAxis()

Monocytes<-c("Itgam", 'Csf1r','Itgax', 'Ccr3', 'Cd14')
Monoinf<-c("Ly6c1",  'Spn',  'Siglecf')
nonclassic <-c("Rhoc", 'Fcgr3','Ms4a7','Cdkn1c', 'Aif1','Cotl1','Fcer1g')
Macrophages<-c("Trem2", 'Nupr1',  'C1qb')
Platelets<-c("Gp9", 'Pf4')
Myeloidpre<- c("Mydgf", 'Cd33')
DC<-c('Lef1','Ccr7','Foxp3','Il2ra', 'Il3ra', 'Cd303')
PlasmaDC<-c( 'Clec4b1', 'Gm14548', 'Gzmb')
DotPlot(Aging_V3, features = c(Monocytes,Monoinf,nonclassic, Macrophages),dot.scale = 7, assay = 'RNA')+ RotatedAxis()
DotPlot(Aging_V3, features = c(Platelets, Myeloidpre, DC, PlasmaDC),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

Tcells<-c("Cd3e", 'Cd4','Cd8a')
CD4<-c("Tcf7", 'Sell','Lef1', 'Ccr7')
Cytotoxic<- c("Gzmb", 'Ptprc','Prf1')
Treg<-c("Foxp3", 'Il2ra')
Tfh<-c("Cxcr5", 'Icos', 'Pdcd1')
ILC<-c("Il7r", 'Areg', 'Klrb1', 'Gata3', 'Itga1')
Naive<-c("Lef1")
NKcells<-c("Ncam1", 'Ncr1','Cd160','Fcgr3a','Ptprc')
DotPlot(Aging_V3, features = c(Tcells,CD4,Cytotoxic),dot.scale = 7, assay = 'RNA')+ RotatedAxis()
DotPlot(Aging_V3, features = c(Treg,Tfh,ILC, Naive,NKcells), assay = 'RNA',dot.scale = 7)+ RotatedAxis()

Bcells<-c('Cd19', 'Il2ra', 'Tnfrsf8')
PlasmaB<-c("Sdc1", 'Xbp1', 'Jchain')
Pre<-c('Mki67')
Epi<-c('Epcam')
MAST<-c('Kit', 'Enpp3')
DotPlot(Aging_V3, features = c(Bcells, PlasmaB, Pre, Epi, MAST),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

Stromal_C1<-c('C1qtnf3')
Stromal_DK<-c("Dkk2")
Stromal_Ac<-c('Ackr4')
Stromal_Cx<-c('Cxcl5')
Stromal_Ha<-c('Has1')
Stromal_Ed<-c('Ednrb')
Lymphatic<-c('Lyve1')
Pericyte<-c('Atcg2', 'Rgs4')
Dividing<-c('Top2a')
Macro<-c('Cd74')
Meso<-c('Slpi')
DotPlot(Aging_V3, features = c(Stromal_C1, Stromal_DK, Stromal_Ac, Stromal_Cx, Stromal_Ha,Stromal_Ed,
                                      Lymphatic, Pericyte, Dividing,Macro, Meso),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

Idents(Aging_V3)<- Aging_V3$seurat_clusters
new.cluster.ids <- c( '0','Enterocytes','Stem cells','Enterocytes','Macrophages','Stem cells',
                      'B cells','Stromal','CD4+ T cells','CD8+ T cells', 'Goblet cells',
                      'Pericyte', 'Tregs','Non-classic Mono','Plasma cells','Lymphatic node',
                      'Enteroendocrine', 'Erythroid cells','Paneth cells','CD8+ T cells','TA cells',
                      'Tuft cells','Macrophages','23')
names(new.cluster.ids) <- levels(Aging_V3)
Aging_V3 <- RenameIdents(Aging_V3, new.cluster.ids)
Aging_V3$CellTypes<-Idents(Aging_V3)

# From the different methods used to call for compartments and cell-types there were 2 clusters present (0 and 23) with a very low feature/count and high percentage of mithocondria.
# We deemed these clusters non-informative has there were ambigous in there gene expression and further remove them from the analysis
Idents(Aging_V3)<- Aging_V3$seurat_clusters
Ageing_V3_clean<-subset(Aging_V3, idents=c('0', '23'), invert=TRUE)
n.cells<-colnames(Aging_V3_clean)
load('Ageing_V3_merged.rda')
Aging_V3_clean<-subset(x = Aging_V3, cells = n.cells)

# we repeated the above analysis
# analysis and separation into compartments
# clustering workflow analysis
all.genes<-rownames(Aging_V3_clean)
Aging_V3_clean <- ScaleData(Aging_V3_clean, features = all.genes)
Aging_V3_clean <- FindVariableFeatures(Aging_V3_clean)
Aging_V3_clean <- RunPCA(Aging_V3_clean,npcs = 80)
ElbowPlot(Aging_V3_clean, ndims = 80)

Aging_V3_clean <- RunUMAP(object = Aging_V3_clean, dims = 1:80)
Aging_V3_clean <- FindNeighbors(Aging_V3_clean,  dims = 1:80)
Aging_V3_clean <- FindClusters(Aging_V3_clean, resolution = 0.2)

# Visualization of merged object
pdf("UMAP_Aging-MB_V3_clean_Sample_Treatment_Mouse.pdf", width = 8, height = 7)
DimPlot(Aging_V3_clean, reduction = "umap",group.by = 'orig.ident', label=F, pt.size = 0.1)
DimPlot(Aging_V3_clean, reduction = "umap",group.by = 'Treatment', label=F, pt.size = 0.1, order = c('yMB', 'iMB'),
        pt.size = 0.1, cols = c("#0072B2", "#D55E00"))
DimPlot(Aging_V3_clean, reduction = "umap",group.by = 'Mouse', label=F, pt.size = 0.1)
dev.off()

pdf("UMAP_Aging-MB_V3_clean_clusters.pdf", width = 13, height = 10)
DimPlot(Aging_V3_clean, reduction = "umap",group.by = 'seurat_clusters', label=TRUE)
dev.off()

# Biomarkers for compartments (Epithelial, Stromal and Immune) in object clusters (res=0.2)
Epi_genes<-c('Epcam', 'Krt8', 'Krt18')
Stromal_genes<-c('Col1a1', 'Col1a2', 'Col6a1','Col6a2', 'Vwf', 'Plvap', 'Cdh5', 'S100b')
Immune_genes<-c('Cd52', 'Cd2', 'Cd3d','Cd3g', 'Cd3e', 'Cd79a', 'Cd79b', 'Cd14', 'Fcgr3','Cd68', 'Cd83', 'Csf1r', 'Fcer1g')
DotPlot(Aging_V3_clean, features = c(Epi_genes,Stromal_genes,Immune_genes),dot.scale = 7, assay='RNA')+ RotatedAxis()

Idents(Aging_V3_clean)<- Aging_V3_clean$seurat_clusters
new.cluster.ids <- c("Epithelium", 'Epithelium','Epithelium','Immune','Epithelium', 'Immune',
                     'Stromal', 'Immune','Immune','Epithelium', 'Stromal',
                     'Immune', 'Immune','Epithelium','Immune','Epithelium',
                     'Stromal','Epithelium', 'Epithelium', 'Immune', 'Immune',
                     'Stromal', 'Epithelium','Immune')
names(new.cluster.ids) <- levels(Aging_V3_clean)
Aging_V3_clean <- RenameIdents(Aging_V3_clean, new.cluster.ids)
Aging_V3_clean$Compartment<-Idents(Aging_V3_clean)

# Visualization of merged object compartments
pdf("UMAP_Aging-MB_V3_clean_Compartment.pdf", width = 8, height = 7)
DimPlot(Aging_V3_clean, reduction = "umap",group.by = 'Compartment', label=F,
        cols = c("bisque3","coral4",   "coral2"))
dev.off()

# SingleR cell-type calling
mouse.se<-MouseRNAseqData()
input <- as.matrix(GetAssayData(object = Aging_V3_clean, slot = "data", assay = "RNA"))
singleR.list <- list()

# perform singleR classification
singleR.list$mouse <- SingleR(test = input, 
                              method="single",
                              fine.tune=FALSE,
                              ref = mouse.se, 
                              labels = mouse.se$label.main)

rm(input)
Ageing_V3_clean$SingleR <- singleR.list$mouse$labels

# Visualization of merged object SingleR cell-types
pdf("UMAP_Aging-MB_V3_cleam_SingleR.pdf", width = 8, height = 6)
DimPlot(object = Aging_V3_clean, reduction = 'umap',  group.by ="SingleR", pt.size = 0.1, label = TRUE)
dev.off()

# Reference cell-type calling base on Tabula muris (2018) facs dataset. To note, this dataset only comprises large intestine samples so there is a mismatch with our object which is comprised of small intestinal epithelium cells
load("~/TM_Intestine_facs_Figure.rds")

reference<-SCTransform(TM_Int,  verbose = FALSE)
query <- SCTransform(Aging_V3_clean,  verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals=FALSE,
  dims = 1:50
)

query <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(
    celltype.l1 = "free_annotation",
    celltype.l2 = "cell_ontology_class",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

# add cell-type metadata from reference
Aging_V3_clean$celltype.l2 <- query$predicted.celltype.l2
Aging_V3_clean$celltype.l1 <- query$predicted.celltype.l1

# visualization of cell-type annotaations base on Tabula Muris reference
pdf("UMAP_Aging-MB_V3_cleam_Celltype_Ref2_free_annotation.pdf", width = 12, height = 6)
DimPlot(Aging_V3_clean, reduction = "umap",group.by = 'celltype.l2', label=TRUE)
dev.off()

pdf("UMAP_Aging-MB_V3_clean_Celltype_Ref1_cell_ontology_class.pdf", width = 12, height = 6)
DimPlot(Aging_V3_clean, reduction = "umap",group.by = 'celltype.l1', label=TRUE)
dev.off()

# Biomarkers for cell-types base on literature (Tabula Muris) and in-house markers
Idents(Aging_V3_clean)<-Aging_V3_clean$seurat_clusters
Stem<-c("Lgr5","Olfm4","Ascl2", "Slc12a1", "Gkn3")
Enterocytes<-c("Alpi","Fabp1","Apoa1", "Apoa4")
Goblet<-c("Muc2","Agr2","Clca3a1", "Tff3", "Clca3a2")
Paneth<-c("Lyz1","Defa17","Defa22", "Defa24", "Ang4")
EE<- c("Chga","Chgb","Tac1", "Tph1")
Tuft<-c("Dclk1","Trpm5","Gfi1b", "Il25")
TA<-c("Stmn1","Tubb5")
PE<-c('Sox4', 'Neurog3', 'Neurod2')
DotPlot(Aging_V3_clean, features = c(Stem,Enterocytes,Goblet, Paneth, EE, Tuft, TA, PE),dot.scale = 7,assay='RNA')+ RotatedAxis()

EnteroMatProx<-c("Apoa4","Fabp1","Apoc2", "Rbp2", "Apoc3", 'Leap2')
EnteroImatProx<-c("Casp6")
EnteroMatDistal<- c("Tmigd1","Fabp6","Slc51a", "Slc51b", "Mep1a", 'Fam151a')
EnteroImatDistal<-c("Reg3g","Gsdmc4","Prss32", "Krt8")
DotPlot(Aging_V3_clean, features = c(EnteroMatProx,EnteroImatProx,EnteroMatDistal, EnteroImatDistal),dot.scale = 7,assay='RNA')+ RotatedAxis()

Monocytes<-c("Itgam", 'Csf1r','Itgax', 'Ccr3', 'Cd14')
Monoinf<-c("Ly6c1",  'Spn',  'Siglecf')
nonclassic <-c("Rhoc", 'Fcgr3','Ms4a7','Cdkn1c', 'Aif1','Cotl1','Fcer1g')
Macrophages<-c("Trem2", 'Nupr1',  'C1qb')
Platelets<-c("Gp9", 'Pf4')
Myeloidpre<- c("Mydgf", 'Cd33')
DC<-c('Lef1','Ccr7','Foxp3','Il2ra', 'Il3ra', 'Cd303')
PlasmaDC<-c( 'Clec4b1', 'Gm14548', 'Gzmb')
DotPlot(Aging_V3_clean, features = c(Monocytes,Monoinf,nonclassic, Macrophages),dot.scale = 7, assay = 'RNA')+ RotatedAxis()
DotPlot(Aging_V3_clean, features = c(Platelets, Myeloidpre, DC, PlasmaDC),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

Tcells<-c("Cd3e", 'Cd4','Cd8a')
CD4<-c("Tcf7", 'Sell','Lef1', 'Ccr7')
Cytotoxic<- c("Gzmb", 'Ptprc','Prf1')
Treg<-c("Foxp3", 'Il2ra')
Tfh<-c("Cxcr5", 'Icos', 'Pdcd1')
ILC<-c("Il7r", 'Areg', 'Klrb1', 'Gata3', 'Itga1')
Naive<-c("Lef1")
NKcells<-c("Ncam1", 'Ncr1','Cd160','Fcgr3a','Ptprc')
DotPlot(Aging_V3_clean, features = c(Tcells,CD4,Cytotoxic),dot.scale = 7, assay = 'RNA')+ RotatedAxis()
DotPlot(Aging_V3_clean, features = c(Treg,Tfh,ILC, Naive,NKcells), assay = 'RNA',dot.scale = 7)+ RotatedAxis()

Bcells<-c('Cd19', 'Il2ra', 'Tnfrsf8')
PlasmaB<-c("Sdc1", 'Xbp1', 'Jchain')
Pre<-c('Mki67')
Epi<-c('Epcam')
MAST<-c('Kit', 'Enpp3')
DotPlot(Aging_V3_clean, features = c(Bcells, PlasmaB, Pre, Epi, MAST),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

Stromal_C1<-c('C1qtnf3')
Stromal_DK<-c("Dkk2")
Stromal_Ac<-c('Ackr4')
Stromal_Cx<-c('Cxcl5')
Stromal_Ha<-c('Has1')
Stromal_Ed<-c('Ednrb')
Lymphatic<-c('Lyve1')
Pericyte<-c('Atcg2', 'Rgs4')
Dividing<-c('Top2a')
Macro<-c('Cd74')
Meso<-c('Slpi')
DotPlot(Aging_V3_clean, features = c(Stromal_C1, Stromal_DK, Stromal_Ac, Stromal_Cx, Stromal_Ha,Stromal_Ed,
                               Lymphatic, Pericyte, Dividing,Macro, Meso),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

Idents(Aging_V3_clean)<- Aging_V3_clean$seurat_clusters
new.cluster.ids <- c("Enterocytes", 'TA cells','Enterocytes','Macrophages','Stem cells', 'B cells',
                     'Stromal', 'CD4+ T cells','CD8+ T cells','Goblet cells', 'Pericytes',
                     'Monocytes', 'Plasma cells','Enterocytes','DC','Erythroid cells',
                     'Lymphatic node','Enteroendocrine', 'Paneth cells', 'DC', 'NK cells',
                     'Stromal', 'Tuft cells','Non-classic Mon')
names(new.cluster.ids) <- levels(Aging_V3_clean)
Aging_V3_clean <- RenameIdents(Aging_V3_clean, new.cluster.ids)
Aging_V3_clean$CellTypes<-Idents(Aging_V3_clean)

save(Aging_V3_clean, file = 'Aging-MB_V3_clean_analyzed.rda')
load('Aging-MB_V3_clean_analyzed.rda')

# Focus on the different compartments (Epithelium and Immune cells only in the manuscript)
split<-SplitObject(Ageing_V3_clean, split.by = 'Compartment')
Epi_Aging_V3<-split$Epithelium
Immu_Aging_V3<-split$Immune

## To facilitate the analysis we provide the high quality cell IDs after multiple rounds of removal of non-informative cells.
# Cell IDs provided for Epithelium
n.cells<-read.csv2('Epithelium_HighQuality_cells.csv')
Epi_Aging_V3[["CellName"]] <- colnames(Epi_Aging_V3)
Epi_Aging_V3_clean <- subset(Epi_Aging_V3, subset = CellName %in% n.cells$colnames.Epi_Ageing_V3_clean. )

# Cell IDs provided for Immune
n.cells<-read.csv2('Immune_HighQuality_cells.csv')
Immu_Aging_V3[["CellName"]] <- colnames(Immu_Aging_V3)
Immu_Aging_V3_clean <- subset(Immu_Aging_V3, subset = CellName %in% n.cells$colnames.Immu_Ageing_V3_clean. )

# analysis of Epithelium compartment separately
# clustering workflow analysis
all.genes<-rownames(Epi_Aging_V3_clean)
Epi_Aging_V3_clean <- ScaleData(Epi_Aging_V3_clean, features = all.genes)

Epi_Aging_V3_clean <- FindVariableFeatures(Epi_Aging_V3_clean)
Epi_Aging_V3_clean <- RunPCA(Epi_Aging_V3_clean, npcs = 80)
ElbowPlot(Epi_Aging_V3_clean, ndims = 80)

Epi_Aging_V3_clean <- RunUMAP(Epi_Aging_V3_clean, dims = 1:60)
Epi_Aging_V3_clean <- FindNeighbors(Epi_Aging_V3_clean,  dims = 1:60)
Epi_Aging_V3_clean <- FindClusters(Epi_Aging_V3_clean, resolution = 0.4)

# Visualization of epithelium object
pdf("UMAP_Epi_Aging-MB_V3_clean_cluster04_Sample_Mouse_Treatment_Tissue.pdf", width = 8, height = 7)
DimPlot(Epi_Aging_V3_clean, reduction = "umap", group.by = 'seurat_clusters', label = TRUE, pt.size = 0.5)
DimPlot(Epi_Aging_V3_clean, reduction = "umap", group.by = 'orig.ident',label=F, pt.size = 0.5)
DimPlot(Epi_Aging_V3_clean, reduction = "umap", group.by = 'Mouse',label=F, pt.size = 0.5)
DimPlot(Epi_Aging_V3_clean, reduction = "umap",group.by = 'Treatment', label=F, 
        order = c('iMB', 'yMB'),
        pt.size = 0.1, cols = c( "#D55E00","#0072B2"))
DimPlot(Epi_Aging_V3_clean, reduction = "umap",group.by = 'Tissue', label=F, 
        order = c('LP','IEC'))
dev.off()


DefaultAssay(Epi_Aging_V3_clean)<-'RNA'
Idents(Epi_Aging_V3_clean)<-Epi_Aging_V3_clean$seurat_clusters
Stem<-c("Lgr5","Olfm4","Ascl2", "Slc12a1", "Gkn3")
Enterocytes<-c("Alpi","Fabp1","Apoa1", "Apoa4")
Goblet<-c("Muc2","Agr2","Clca3a1", "Tff3", "Clca3a2")
Paneth<-c("Lyz1","Defa17","Defa22", "Defa24", "Ang4")
EE<- c("Chga","Chgb","Tac1", "Tph1")
Tuft<-c("Dclk1","Trpm5","Gfi1b", "Il25")
TA<-c("Stmn1","Tubb5")
PE<-c('Sox4', 'Neurog3', 'Neurod2')
DotPlot(Epi_Aging_V3_clean, features = c(Stem,Enterocytes,Goblet, Paneth, EE, Tuft, TA, PE),dot.scale = 7,assay='RNA')+ RotatedAxis()

EnteroMatProx<-c("Apoa4","Fabp1","Apoc2", "Rbp2", "Apoc3", 'Leap2')
EnteroImatProx<-c("Casp6")
EnteroMatDistal<- c("Tmigd1","Fabp6","Slc51a", "Slc51b", "Mep1a", 'Fam151a')
EnteroImatDistal<-c("Reg3g","Gsdmc4","Prss32", "Krt8")
DotPlot(Epi_Aging_V3_clean, features = c(EnteroMatProx,EnteroImatProx,EnteroMatDistal, EnteroImatDistal),dot.scale = 7,assay='RNA')+ RotatedAxis()


Idents(Epi_Aging_V3_clean)<-Epi_Aging_V3_clean$seurat_clusters
new.cluster.ids <- c("Stem cells (Lgr5+)", 'Enterocytes','Enterocytes','TA cells','Enterocytes', 'TA cells',
                     'Immature Enterocytes','Enterocytes Proximal','Goblet cells','Immature Enterocytes','Enteroendocrine', 
                     'Enterocytes Distal','Paneth cells','Goblet cells', 'Enterocytes', 'Tuft cells',
                     'Paneth cells','Enterocytes','Enterocytes Proxima')
names(new.cluster.ids) <- levels(Epi_Aging_V3_clean)
Epi_Aging_V3_clean <- RenameIdents(Epi_Aging_V3_clean, new.cluster.ids)
Epi_Aging_V3_clean$Celltype<-Idents(Epi_Aging_V3_clean)


# visualization of cell-types of epithelium compartment
singleR_colors <- c("Enterocytes"='#b35806',
                    "Immature Enterocytes"='#e3a064',
                    "Enterocytes Proximal"='#8b6914',
                    'Enterocytes Distal'='#bf812d',
                    'Paneth cells'='#27408b',
                    'Goblet cells'='#556b2f',
                    'Enteroendocrine'='#c2426b',
                    'Tuft cells'='#6A3D8A',
                    'TA cells'= '#bc6ead',
                    'Stem cells (Lgr5+)'='#ffe1ff')

pdf("UMAP_Epi_Aging-MB_V3_clean_Celltype_Figure.pdf", width = 9, height = 7)
DimPlot(Epi_Aging_V3_clean, reduction = "umap", group.by = 'Celltype',label=TRUE, pt.size = 0.5, 
        cols = singleR_colors)
dev.off()

save(Epi_Aging_V3_clean,  file ='Epithelium_Aging-MB_V3_clean_Figure.rda' )


# visualization of cell-types biomarkers of epithelium compartment
marker.genes <- rev(c('Apoa4','Casp6','Fabp6','Prss32','Stmn1','Tubb5','Lyz1','Agr2','Dclk1','Chga'))
pdf("Vlnplot_Epi_Ageing-MB_V3_clean_Celltype.pdf", width = 5, height = 3.5)
VlnPlot(object = Epi_Aging_V3_clean, 
        features = rev(marker.genes), 
        group.by ="Celltype",
        assay = "RNA", stack = TRUE, fill.by = "ident",
        cols= singleR_colors)+
  scale_y_discrete(limits=rev(c('Enterocytes', 'Enterocytes Proximal','Enterocytes Distal',
                                'Immature Enterocytes','Stem cells (Lgr5+)','TA cells',
                                'Paneth cells','Goblet cells','Tuft cells','Enteroendocrine')))+
  theme(axis.text.x.bottom = element_text(hjust=1),
        axis.title.x.top = element_text(angle=60)) + NoLegend()
dev.off()



# analysis of Immune compartment separately
# clustering workflow analysis
all.genes<-rownames(Immu_Aging_V3_clean)
Immu_Aging_V3_clean <- ScaleData(Immu_Aging_V3_clean, features = all.genes)
Immu_Aging_V3_clean <- FindVariableFeatures(Immu_Aging_V3_clean)
Immu_Aging_V3_clean <- RunPCA(Immu_Aging_V3_clean, npcs = 80)
ElbowPlot(Immu_Aging_V3_clean, ndims = 80)

Immu_Aging_V3_clean <- RunUMAP(Immu_Aging_V3_clean, dims = 1:50)
Immu_Aging_V3_clean <- FindNeighbors(Immu_Aging_V3_clean,  dims = 1:50)
Immu_Aging_V3_clean <- FindClusters(Immu_Aging_V3_clean, resolution = 0.4)

# Visualization of immune object
pdf("UMAP_Immu_Aging-MB_V3_clean_cluster04_Sample_Mouse_Treatment_Tissue.pdf", width = 8, height = 7)
DimPlot(Immu_Aging_V3_clean, reduction = "umap",  group.by = 'seurat_clusters', label = TRUE, pt.size = 0.5)
DimPlot(Immu_Aging_V3_clean, reduction = "umap", group.by = 'orig.ident',label=F, pt.size = 0.5)
DimPlot(Immu_Aging_V3_clean, reduction = "umap", group.by = 'Mouse',label=F, pt.size = 0.5)
DimPlot(Immu_Aging_V3_clean, reduction = "umap",group.by = 'Treatment', label=F, 
        order = c('iMB', 'yMB'),
        pt.size = 0.1, cols = c("#D55E00","#0072B2" ))
DimPlot(Immu_Aging_V3_clean, reduction = "umap",group.by = 'Tissue', label=F, 
        order = c('LP','IEC'))
dev.off()

# SingleR cell-type annotation
mouse.se<-MouseRNAseqData()
input <- as.matrix(GetAssayData(object = Immu_Aging_V3_clean, slot = "data", assay = "RNA"))
singleR.list <- list()
# perform singleR classification
singleR.list$mouse <- SingleR(test = input, 
                              method="single",
                              fine.tune=FALSE,
                              ref = mouse.se, 
                              labels = mouse.se$label.main)

rm(input)
Immu_Aging_V3_clean$SingleR <- singleR.list$mouse$labels
pdf("UMAP_Immu_Aging-MB__V3_clean_SingleR.pdf", width = 8, height = 6)
DimPlot(object = Immu_Aging_V3_clean, reduction = 'umap',  group.by ="SingleR", pt.size = 0.5, label = F)
dev.off()

DefaultAssay(Immu_Aging_V3_clean)<-'RNA'
Idents(Immu_Aging_V3_clean)<-Immu_Aging_V3_clean$seurat_clusters
Monocytes<-c("Itgam", 'Csf1r','Itgax', 'Ccr3', 'Cd14')
Monoinf<-c("Ly6c1",  'Spn',  'Siglecf')
nonclassic <-c("Rhoc", 'Fcgr3','Ms4a7','Cdkn1c', 'Aif1','Cotl1','Fcer1g')
Macrophages<-c("Trem2", 'Nupr1',  'C1qb')
Platelets<-c("Gp9", 'Pf4')
Myeloidpre<- c("Mydgf", 'Cd33')
DC<-c('Lef1','Ccr7','Foxp3','Il2ra')
PlasmaDC<-c("Il3ra", 'Clec4b1', 'Gm14548', 'Gzmb')
DotPlot(Immu_Aging_V3_clean, features = c(Monocytes,Monoinf,nonclassic, Macrophages),dot.scale = 7, assay = 'RNA')+ RotatedAxis()
DotPlot(Immu_Aging_V3_clean, features = c(Platelets, Myeloidpre, DC, PlasmaDC),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

Tcells<-c("Cd3e", 'Cd4','Cd8a')
CD4<-c("Tcf7", 'Sell','Lef1', 'Ccr7')
Cytotoxic<- c("Gzmb", 'Ptprc','Prf1')
Treg<-c("Foxp3", 'Il2ra')
Tfh<-c("Cxcr5", 'Icos', 'Pdcd1')
ILC<-c("Il7r", 'Areg', 'Klrb1', 'Gata3', 'Itga1')
Naive<-c("Lef1")
NKcells<-c("Ncam1", 'Ncr1','Cd160','Fcgr3a','Ptprc')
DotPlot(Immu_Aging_V3_clean, features = c(Tcells,CD4,Cytotoxic),dot.scale = 7, assay = 'RNA')+ RotatedAxis()
DotPlot(Immu_Aging_V3_clean, features = c(Treg,Tfh,ILC, Naive,NKcells), assay = 'RNA',dot.scale = 7)+ RotatedAxis()

Bcells<-c('Cd19', 'Il2ra', 'Tnfrsf8')
PlasmaB<-c("Sdc1", 'Xbp1', 'Jchain')
Pre<-c('Mki67')
Epi<-c('Epcam')
Granu<-c('Csf3r')
MAST<-c('Kit', 'Enpp3','Tpsab1', 'Cpa3','Srgn')
Ery<-c('Hba-a1', 'Hba-a2')
DotPlot(Immu_Aging_V3_clean, features = c(Bcells, PlasmaB, Pre, Epi, MAST, Ery, Granu),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

ILC1<-c('Il7r', 'Klrb1c', 'Ncr1', 'Itgae', 'Il1r1')
ILC2<-c( 'Kit', 'Thy1')
ILC3 <-c( 'Il23r')
DotPlot(Immu_Aging_V3_clean, features = c(ILC1, ILC2, ILC3),dot.scale = 7, assay = 'RNA')+ RotatedAxis()

Idents(Immu_Aging_V3_clean)<-Immu_Aging_V3_clean$seurat_clusters
new.cluster.ids <- c("B cells", 'Macrophages','CD8+ T cells','Plasma cells','Nonclassic', 'Macrophages',
                     'Erythroid cells', 'CD4+ T cells','MAIT cells','Tregs','ILC2', 
                     'ILC3','DC','NKT cells', 'CD4+ T cells', 'ILC1',
                     'pDC','CD8+ T cells','Prolif Lympho','Plasma cells','Prolif Lympho',
                     'Granulocytes','B cells' )
names(new.cluster.ids) <- levels(Immu_Aging_V3_clean)
Immu_Aging_V3_clean <- RenameIdents(Immu_Aging_V3_clean, new.cluster.ids)
Immu_Aging_V3_clean$Celltype<-Idents(Immu_Aging_V3_clean)

Idents(Immu_Aging_V3_clean)<-Immu_Aging_V3_clean$Celltype
Immu_Aging_V3_clean_sub<-subset(Immu_Aging_V3_clean, idents='Erythroid cells', invert=TRUE)

# Visualization of cell-types of immune compartment
singleR_colors <- c( "B cells" = "goldenrod1",
                     "Plasma cells" = "goldenrod4",
                     "CD4+ T cells" = "#EA811F",
                     "CD8+ T cells" = "maroon",
                     "ILC1" = "#6A3D8A",
                     "ILC2" = "#8a6a3d",
                     "ILC3" = "#3d8a6a",
                     "NKT cells" = "purple",
                     "Tregs"="#462609",
                     "MAIT cells" = "#1A3865",
                     "Macrophages" = "indianred",
                     "DC" = "orchid3",
                     "pDC"="orchid4",
                     "Nonclassic" = "#577676",
                     "Granulocytes" = "#1A3865",
                     'Prolif Lympho'='#8fdbaf')


pdf("UMAP_Immu_Ageing-MB__V3_clean_Celltype_Figure_sub.pdf", width = 9, height = 7)
DimPlot(Immu_Aging_V3_clean_sub, reduction = "umap", group.by = 'Celltype',label=TRUE, pt.size = 0.5, 
        cols = singleR_colors)
dev.off()


# visualization of cell-types biomarkers of immune compartment
marker.genes <- rev(c('Cd19','C3ar1','Cd8a','Cd4','Sdc1','Aif1','Eomes','Foxp3','Il7r','Thy1','Il23r',
                      'Cxcr3','Flt3','Ccr7','Mki67','Csf3r','Tfrc'))

pdf("Vlnplot_Immune_Ageing-MB_V3_clean_Celltype_sub.pdf", width = 5, height = 4.5)
VlnPlot(object = Immu_Aging_V3_clean_sub, 
        features = rev(marker.genes), 
        group.by ="Celltype",
        assay = "RNA", stack = TRUE, fill.by = "ident",
        cols= singleR_colors)+
  scale_y_discrete(limits=rev(c("B cells", 'Macrophages','CD8+ T cells','CD4+ T cells','Plasma cells','Nonclassic', 
                                'MAIT cells','Tregs','ILC1','ILC2', 
                                'ILC3','NKT cells','DC',
                                'pDC','Prolif Lympho',
                                'Granulocytes')))+
  theme(axis.text.x.bottom = element_text(hjust=1),
        axis.title.x.top = element_text(angle=60)) + NoLegend()
dev.off()

save(Immu_Aging_V3_clean_sub,  file ='Immune_Aging-MB_V3_clean_Figure.rda' )


# The following analysis was performed for both Epithelium and Immune compartments. To simplify we here provide the code for Epithelium object, but it could be substitute for Immune object.
# cell proportion differences between treatments (yMB vs. iMB)
## Calculation of cell numbers per mouse, per treatment for each identified cell-types
k=as.data.frame(table(Epi_Aging_V3_clean$Mouse, Epi_Aging_V3_clean$Treatment,Epi_Aging_V3_clean$Celltype))
write.csv2(k, 'Ref_doublet_Epi_Aging-MB_V3_clean_anno.csv')

# calculate cell-type proportions per mouse
Cell<-dcast(data = k,formula = Var1~Var3,fun.aggregate = sum,value.var = "Freq")
rownames(Cell)<-Cell$Var1
Cell$Var1=NULL
# normalize cell proportions (relative count)
Cell$total<-rowSums(Cell)
cell_props2<-as.data.frame(Cell/Cell$total)*100
cell_props2$total<-NULL
cell_props2$Treatment=c('iMB', 'yMB', 'yMB', 'yMB', 'iMB', 'iMB')
cell_props2$Mouse<-rownames(Cell)

cell2<-melt(cell_props2, id.vars=c('Mouse','Treatment'))
cell2$celltypes<-cell2$variable
cell2$variable<-NULL

df_p_val <- cell2 %>% 
  rstatix::group_by(celltypes) %>% 
  rstatix::t_test(value ~ Treatment) %>% 
  rstatix::add_xy_position(step.increase = 0.05)

#  visualization of cell proportion per treatment
pdf('Cellprop_SuppFigure1_detailed.pdf', width = 8, height = 5)
g <- ggplot(cell2, aes(x =Treatment, y =value))+
  geom_boxplot(aes(fill=Treatment), width = 0.5, )+
  theme_bw()+
  facet_wrap(vars(celltypes),scales = "free_y",   ncol = 4) +
  scale_fill_manual(values = c("yMB" = "#D55E00",
                               "iMB" = "#0072B2"))+
  scale_x_discrete(limits=c('iMB', 'yMB'))
g <- g + theme(panel.spacing = unit(0.25, "cm"))
g <- g + labs(x ="Treatment", y = "Cell proportions (%)")
g <- g + add_pvalue(df_p_val,label = "p",tip.length = 0, label.size = 2 )
g
dev.off()

# calculation of differentially expressed genes for each individual cell-type 
Epi_Ageing_V3_clean$Celltype_Treatment<-paste(Epi_Ageing_V3_clean$Celltype, Epi_Ageing_V3_clean$Treatment, sep='_')
Idents(Epi_Ageing_V3_clean)<-Epi_Ageing_V3_clean$Celltype_Treatment

Enterocytes_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Enterocytes_iMB', ident.2 = 'Enterocytes_yMB')
Enterocytes.Proximal_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Enterocytes Proximal_iMB', ident.2 = 'Enterocytes Proximal_yMB')
Enterocytes.Distal_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Enterocytes Distal_iMB', ident.2 = 'Enterocytes Distal_yMB')
Enterocytes.Immature_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Immature Enterocytes_iMB', ident.2 = 'Immature Enterocytes_yMB')
TA_cells_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'TA cells_iMB', ident.2 = 'TA cells_yMB')
Stem_cellsLgr5_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Stem cells (Lgr5+)_iMB', ident.2 = 'Stem cells (Lgr5+)_yMB')
Goblet_cells_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Goblet cells_iMB', ident.2 = 'Goblet cells_yMB')
Paneth_cells_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Paneth cells_iMB', ident.2 = 'Paneth cells_yMB')
Enteroendocrine_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Enteroendocrine_iMB', ident.2 = 'Enteroendocrine_yMB')
Tuft_Ageing <- FindMarkers(Epi_Ageing_V3_clean, ident.1 = 'Tuft cells_iMB', ident.2 = 'Tuft cells_yMB')


Enterocytes_Ageing$cluster<-'Enterocytes'
Enterocytes.Proximal_Ageing$cluster<-'EnteroP'
Enterocytes.Distal_Ageing$cluster<-'EnteroD'
Enterocytes.Immature_Ageing$cluster<-'ImmEntero'
TA_cells_Ageing$cluster<-'TA'
Stem_cellsLgr5_Ageing$cluster<-'StemLgr5'
Goblet_cells_Ageing$cluster<-'Goblet'
Paneth_cells_Ageing$cluster<-'Paneth'
Enteroendocrine_Ageing$cluster<-'EE'
Tuft_Ageing$cluster<-'Tuft'

Enterocytes_Ageing$genes<-rownames(Enterocytes_Ageing)
Enterocytes.Proximal_Ageing$genes<-rownames(Enterocytes.Proximal_Ageing)
Enterocytes.Distal_Ageing$genes<-rownames(Enterocytes.Distal_Ageing)
Enterocytes.Immature_Ageing$genes<-rownames(Enterocytes.Immature_Ageing)
TA_cells_Ageing$genes<-rownames(TA_cells_Ageing)
Stem_cellsLgr5_Ageing$genes<-rownames(Stem_cellsLgr5_Ageing)
Goblet_cells_Ageing$genes<-rownames(Goblet_cells_Ageing)
Paneth_cells_Ageing$genes<-rownames(Paneth_cells_Ageing)
Enteroendocrine_Ageing$genes<-rownames(Enteroendocrine_Ageing)
Tuft_Ageing$genes<-rownames(Tuft_Ageing)

Total<-do.call('rbind', list(Enterocytes_Ageing,Enterocytes.Proximal_Ageing, Enterocytes.Distal_Ageing,Enterocytes.Immature_Ageing,
                             TA_cells_Ageing,Stem_cellsLgr5_Ageing, Goblet_cells_Ageing, Paneth_cells_Ageing,
                             Enteroendocrine_Ageing,Tuft_Ageing ))
a<-subset(Total, Total$p_val_adj<0.05)
Up<-subset(a, a$avg_log2FC < 0)
Up$signal<-'Upregulated_yMB'
Down<-subset(a, a$avg_log2FC > 0)
Down$signal<-'Downregulated_yMB'
Total<-rbind(Up, Down)
write.csv2(Total, file='DEGs_Celltypes_Epithelium_yMBvsiMB_new.csv')

Total<-read.csv2('DEGs_Aging-MB_Celltypes_Epithelium_yMBvsiMB.csv')

# visualization of differentially expressed genes as barplot divided by down and upregulated genes
Total$group2 <- paste(Total$cluster, Total$signal, sep = "_")

Total2 <- Total %>% 
  group_by(cluster, signal) %>%
  dplyr::count(group2)

pdf("DEGs_Freq_Epithelium_Ageing_V3_clean_Celltypes.pdf", width = 5, height = 4)
ggplot(Total2, aes(x = cluster, y =n , fill = signal)) +
  #geom_point(aes(size = Count)) +
  geom_bar(colour='black',stat="identity", width=0.9,    
           position=position_dodge(0.9))+
  scale_fill_manual(values=c("cornflowerblue", "brown1")) +
  coord_flip()+
  theme_bw() + 
  scale_x_discrete(limits=rev(c("Enterocytes" ,"EnteroP","EnteroD" ,
                                "ImmEntero","StemLgr5",
                                "TA","Paneth", "Goblet" , "Tuft" ,
                                "EE"  )))+
  theme(text = element_text(size=14),axis.text.x = element_text(size=12,vjust = 0.5, hjust=1), axis.text.y = element_text(size=12))+
  labs(title="DEGs in yMB",
       x ="Cell-types", y = "Number of DE genes")
dev.off()

# visualization of common and unique differentially expressed genes by the means of a upset plot
#1201 unique genes, 2407 shared with at least 2
Up<-subset(Total, Total$signal %in% 'Upregulated_yMB')
Up2<-Up[,7:8]
Up3<-split( Up2 , f = Up$cluster, drop = TRUE )
y<-lapply(Up3, function(x) x[!(names(x) %in% c("cluster"))])

#3380 unique genes, 4465 shared with at least 2
Down<-subset(Total, Total$signal %in% 'Downregulated_yMB')
Down2<-Down[,7:8]
Down3<-split( Down2 , f = Down$cluster, drop = TRUE )
y<-lapply(Down3, function(x) x[!(names(x) %in% c("cluster"))])

x <- list(
  Entero = y$Enterocytes$genes, 
  ImmaEntero = y$ImmEntero$genes, 
  TA = y$TA$genes,
  Goblet = y$Goblet$genes,
  Paneth = y$Paneth$genes,
  EE = y$EE$genes
)
m = make_comb_mat(x)
ss = set_size(m)
cs = comb_size(m)
ht = UpSet(m, 
           comb_col = c("brown1"),
           set_order = order(ss),
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Gene Intersections" = anno_barplot(cs, 
                                                 ylim = c(0, max(cs)*1.1),
                                                 border = FALSE, 
                                                 gp = gpar(fill = "black"), 
                                                 height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           
           left_annotation = rowAnnotation(
             set_name = anno_text(set_name(m), 
                                  location = 0.6, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)
ht = draw(ht)

pdf('Upset_Count_UpCelltypes_Epithelium_yMBvsiMB_new.pdf', width = 7, height = 3)
ht
od = column_order(ht)
decorate_annotation("Gene Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
dev.off()

# gene ontology enrichment analysis of differentially expressed genes between treatment using comparison function generated by cluster profiler 
CompGO <- as.data.frame(compareCluster(ENTREZID~cluster+signal, data=Total,OrgDb='org.Mm.eg.db', fun="enrichGO", ont = "BP"))

ENTREZID<-as.data.frame(CompGO$geneID)
nmax <- max(stringr::str_count(ENTREZID$`CompGO$geneID`, "\\/")) + 1
ENTREZID2<-tidyr::separate(ENTREZID, `CompGO$geneID`, paste0("col", seq_len(nmax)), sep = "\\/", fill = "right")
rows<- nrow(ENTREZID2)
GENEID <- as.data.frame(matrix(mapIds(org.Mm.eg.db, as.character(unlist(ENTREZID2)), "SYMBOL","ENTREZID"), rows))
col<-colnames(GENEID)
GENEID2<-tidyr::unite(GENEID, col='GENEID', col, sep=',')
GENEID2$GENEID<-gsub("NULL,","",as.character(GENEID2$GENEID))
GENEID2$GENEID<-gsub(",NULL","",as.character(GENEID2$GENEID))
CompGO$GENEID<-GENEID2$GENEID


# save only top 25 terms for each cell-type
CompGO_sub<-CompGO %>% group_by(Cluster) %>% top_n(n= 25, wt = -p.adjust)
write.csv2(CompGO_sub, file='GOComparison_Aging-MB_Celltypes_Epithelium_yMBvsiMB.csv')

# visualization of curated GO term list base on hallmarks of aginf literature
CompGO<-read.csv2('GOComparison_Aging-MB_Celltypes_Epithelium_yMBvsiMB_interest.csv')
CompGO<-read.csv2('GOComparison_Aging-MB_Celltypes_Immune_yMBvsiMB_interest.csv')

CompGO %>% group_by(Cluster) %>% top_n(n= 10, wt = -p.adjust) -> terms
CompGOresults_sub<-subset(CompGO, CompGO$p.adjust <0.05)
CompGOresults_sub2<-subset(CompGOresults_sub, CompGOresults_sub$Count>=3)
CompGOresults_sub2$Description <- ifelse(nchar(CompGOresults_sub2$Description)>60,
                                         paste(substr(CompGOresults_sub2$Description, 1, 60),"[...]",sep=""),
                                         CompGOresults_sub2$Description)
CompGOresults_sub2$Description <- factor(CompGOresults_sub2$Description,levels=unique(CompGOresults_sub2$Description))
CompGOresults_sub3<-CompGOresults_sub2 %>% group_by(Cluster) %>% top_n(n= 5, wt = -p.adjust)
CompGOresults_sub4<-CompGOresults_sub3 %>% arrange(factor(cluster))
CompGOresults_sub5<-subset(CompGOresults_sub4, CompGOresults_sub4$cluster %in% c('EE', 'Enterocytes','ImmEntero','TA'))
CompGOresults_sub5<-subset(CompGOresults_sub4, CompGOresults_sub4$cluster %in% c('Macrophages', 'B cells','CD8+ T cells','MAIT cells'))

pdf("GOComparison_Celltypes_interest_Epithelium_yMBvsiMB_multiple_interest.pdf", width = 14, height =6)
ggplot(CompGOresults_sub5, aes(x = p.adjust, y = Description,size = Count, fill=signal)) +
  geom_point(shape=21) +
  scale_fill_manual(values=c("cornflowerblue", "brown1"))+
  scale_x_Pvalue()+
  ylab(NULL) +
  theme_bw() +
  facet_wrap(~ cluster, ncol = 2, scales = "free") +
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


# TFBS gene enrichment analysis for the differentially expressed genes between treatments (only for epithelium)
motifRankings_500bp <- importRankings("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
motifRankings_10kbp <- importRankings("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
data(motifAnnotations_mgi)

### subset for the cell types of interest
Idents(Epi_Aging_V3_clean) <- Epi_Aging_V3_clean$Celltype_Treatment
Interest<-subset(Epi_Aging_V3_clean, idents=c('Enterocytes_yMB', 'Enterocytes_iMB','TA cells_yMB', 'TA cells_iMB',
                                               'Immature Enterocytes_yMB','Immature Enterocytes_iMB','Paneth cells_yMB', 'Paneth cells_iMB','Enteroendocrine_yMB', 'Enteroendocrine_iMB'))

## Identify signature genes characteristic of the group comprised of cell-type and treatment
Interest_markers <- FindAllMarkers(object = Epi_Aging_V3_clean,
                                   only.pos = TRUE,
                                   min.pct = 0.2,
                                   assay = "RNA", 
                                   logfc.threshold = 0.15,
                                   min.diff.pct = 0.1,
                                   test.use = "wilcox"
)
Interest_markers2<-subset(Interest_markers, Interest_markers$p_val_adj < 0.5)

## enrichment analysis for TFBS
Rcistarget.list <- list()
for(i in unique(Interest_markers2$cluster)){
  print(i)
  markers <- Interest_markers2[Interest_markers2$cluster==i,]
  print(paste("cluster: ",i, ", marker genes: ",paste(markers$gene[1:10],collapse=", "),", ...",sep=""))
  
  motif_enrichment <- cisTarget(geneSets = markers$gene, 
                                motifRankings = motifRankings_500bp,
                                motifAnnot = motifAnnotations,
                                nesThreshold = 3,
                                motifAnnot_highConfCat = "directAnnotation"
  )
  
  if(nrow(motif_enrichment)>0){motif_enrichment$cluster <- as.character(i)}
  if(nrow(motif_enrichment)>0){Rcistarget.list[[paste(i)]] <- motif_enrichment}
}
sapply(Rcistarget.list, nrow)

enrichTF_list <- list()
for(j in 1: length(Rcistarget.list)){
  print(names(Rcistarget.list)[j])
  data <- Rcistarget.list[[j]]
  data <- data[order(data$NES, decreasing = T),]
  data$name <- data$TF_highConf
  data$name <- gsub('\\(.*?\\). ', '', data$name)
  data$description <- paste(ifelse(nchar(data$name)>50,
                                   paste(substr(data$name, 1, 50),"[...]",sep=""),
                                   data$name), 
                            "(",data$motif,")",sep="")
  tmp <- data[!duplicated(data$name),]
  data$plot <- ifelse(data$description %in% tmp$description, "plot highest NES per TF", "")
  data$cluster <- paste(names(Rcistarget.list)[j])
  data <- data[!data$TF_highConf == "",]
  enrichTF_list[[paste(names(Rcistarget.list)[j])]] <- data
}
# All databases
enrichTF <- do.call(rbind, enrichTF_list)
write.csv2(enrichTF, file='TFenrich_Epi_Ageing_V3_clean_Celltypes_Treatment_motifRankings_10k.csv')
write.csv2(enrichTF, file='TFenrich_Epi_Ageing_V3_clean_Celltypes_Treatment_motifRankings_500bp.csv')

Bp500<-read_csv2('TFenrich_Epi_Ageing_V3_clean_Celltypes_Treatment_motifRankings_500bp.csv')
Bp10k<-read_csv2('TFenrich_Epi_Ageing_V3_clean_Celltypes_Treatment_motifRankings_10k.csv')

# focused on the common TF sites identified using both windows
Common<-intersect(Bp500$TF_highConf, Bp10k$TF_highConf)
Common2<-subset(Bp500, Bp500$TF_highConf %in% Common)
TF_top <- Common2[order(Common2$cluster),]
TF_top$description <- factor(TF_top$description, levels = unique(TF_top$description))
TF_top2<-subset(TF_top, TF_top$plot %in% 'plot highest NES per TF')

# visualization of enriched TF
TF_top2$Description<-str_replace(TF_top2$description, " \\s*\\([^\\)]+\\)", "")
TF_top_all<-arrange(TF_top2, cluster)
TF_top_all2<-subset(TF_top_all, TF_top_all$cluster %in% c('Enterocytes_iMB','Enterocytes_yMB','Enterocytes Distal_iMB', 'Enterocytes Distal_yMB',
                                                          'Immature Enterocytes_iMB','Immature Enterocytes_yMB','TA cells_iMB', 'TA cells_yMB',
                                                          'Paneth cells_iMB', 'Paneth cells_yMB','Goblet cells_iMB', 'Goblet cells_yMB',
                                                          'Enteroendocrine_iMB', 'Enteroendocrine_yMB'))

pdf("TFenrich_Epi_Ageing_V3_clean_Interest_top10_motifRankings_Common.pdf", width = 4, height =4.5)
ggplot(TF_top_all2, aes(x = cluster, y = Description, color = NES)) +
  geom_point(aes(size = nEnrGenes)) + 
  scale_x_discrete(limits=c('Enterocytes_iMB','Enterocytes_yMB','Enterocytes Distal_iMB', 'Enterocytes Distal_yMB',
                            'Immature Enterocytes_iMB','Immature Enterocytes_yMB','TA cells_iMB', 'TA cells_yMB',
                            'Paneth cells_iMB', 'Paneth cells_yMB','Goblet cells_iMB', 'Goblet cells_yMB',
                            'Enteroendocrine_iMB', 'Enteroendocrine_yMB'))+
  ylab(NULL) +
  theme_bw() + 
  theme(text = element_text(size=12))+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

## cell-cell communication library (CellChat)
# initiatate the cell-cell communication calculation for epithelium and immune (only secreted) compatments seperately for yMB and iMB
Epi<-SplitObject(Epi_Aging_V3_clean, split.by = 'Treatment')
yMB<-Epi$yMB
iMB<-Epi$iMB

# example of individual object calculations
Ind<-yMB
meta<-Ind@meta.data
meta$Celltype = droplevels(meta$Celltype, exclude = setdiff(levels(meta$Celltype),unique(meta$Celltype)))
data.input <- Ind[["RNA"]]@data
# input is a Seurat object
## use the default cell identities of Seurat object
Ind_cellchat <- createCellChat(object = data.input, group.by = 'Celltype', assay = 'RNA', do.sparse = TRUE)

Ind_cellchat <- addMeta(Ind_cellchat, meta = meta)
Ind_cellchat <- setIdent(Ind_cellchat, ident.use = "Celltype") # set "labels" as default cell identity
levels(Ind_cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(Ind_cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the strCDture of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search =  list(c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact" , "Non-protein Signaling"))) # use all
# use all CellChatDB for cell-cell communication analysis

# set the used database in the object
Ind_cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
Ind_cellchat <- subsetData(Ind_cellchat) # This step is necessary even if using the whole database
Ind_cellchat <- identifyOverExpressedGenes(Ind_cellchat)
Ind_cellchat <- identifyOverExpressedInteractions(Ind_cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
Ind_cellchat <- computeCommunProb(Ind_cellchat, type = "triMean",population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
Ind_cellchat <- filterCommunication(Ind_cellchat, min.cells = 50)
Ind_cellchat <- computeCommunProbPathway(Ind_cellchat)
Ind_cellchat <- aggregateNet(Ind_cellchat)

groupSize <- as.numeric(table(Ind_cellchat@idents))

pdf("CellCell_inter_yMB_Epithelium_cellchat.pdf", width = 6, height = 6)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(Ind_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Ind_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("CellCell_inter_yMB_Epithelium_cellchat2.pdf", width = 6, height = 6)
mat <- Ind_cellchat@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# Compute the network centrality scores
Ind_cellchat <- netAnalysis_computeCentrality(Ind_cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
yMB_cellchat<-Ind_cellchat
save(yMB_cellchat, iMB_cellchat, file='Epithelium_Ageing_V3_clean_Cellchat.rda')

# Compare yMB vs. iMB
load('Epithelium_Ageing_V3_clean_Cellchat.rda')
yMB_cellchat<-updateCellChat(yMB_cellchat)
iMB_cellchat<-updateCellChat(iMB_cellchat)

object.list <- list( iMB = iMB_cellchat,yMB = yMB_cellchat)
Comp_cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
Comp_cellchat
# Compare total number of interactions and interaction strenght
gg1 <- compareInteractions(Comp_cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(Comp_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("CellCell_interactions_Comp_Epithelium_Ageing_V3_clean.pdf", width = 3, height = 2)
gg1 + gg2
dev.off()

# compare differential number of interactions: red increase in remission and blue dcrease in remission
par(mfrow = c(1,2), xpd=TRUE)
pdf("CellCell_interactions_Comp_Epithelium_Ageing_V3_clean2.pdf", width = 6, height = 6)
netVisual_diffInteraction(Comp_cellchat, weight.scale = T, measure = 'count')
netVisual_diffInteraction(Comp_cellchat, weight.scale = T, measure = "weight")
dev.off()


#cell population with sig changes
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+
    scale_y_continuous(limits = c(0,0.02)) +  scale_x_continuous(limits = c(0,0.02))
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf('CellCell_Celltypes_Comp_epithelium_Ageing_V3_clean.pdf', width = 6.5, height = 2.5)
patchwork::wrap_plots(plots = gg)
dev.off()


# Gene of interest general expression
## epithelium genes and mesenchymal markers
Strogenes<-c('Vim','Ctnnb1','Fn1','Aifm2','Tgfb1', 'Tgfbr1','Smad2','Smad3','Smad4')

# Calculation of gene module for literature genes related to EMT
Epi_Aging_V3_clean <- AddModuleScore(Epi_Aging_V3_clean, features = list(Strogenes), name="Stromal", assay = 'RNA')

pdf('Vln_Epithelium_stromalScore_Celtype_Treatment.pdf', width = 10, height = 5)
Strom<-VlnPlot(Epi_Aging_V3_clean, features = c( "Stromal1"), split.by = 'Treatment',
               group.by = 'Celltype',pt.size = 0, col=c("#0072B2", "#D55E00"))+
  scale_x_discrete(limits=c('Enterocytes', 'Enterocytes Proximal','Enterocytes Distal',
                            'Immature Enterocytes','Stem cells (Lgr5+)','TA cells',
                            'Paneth cells','Goblet cells','Tuft cells','Enteroendocrine'))
Strom
dev.off()

Test<-Strom[[1]][["data"]]
res<-Test  %>%
  group_by(ident) %>%
  wilcox_test(data =., Stromal1~split ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
res


## immune genes and inflammatory score markers

Chemo<-as.data.frame(c('Tnf', 'Ifng', 'Il1b','Il2','Il6','Cxcl15', 'Ccl20', 'Ccl9','Ccr1','Nfkb1','Myd88','Tlr6'))

Immu_Aging_V3_clean <- AddModuleScore(Immu_Aging_V3_clean, features = list(Chemo), name="Chemo", assay = 'RNA')

pdf('Vln_Immune_InflaScore_Celltype_Treatment.pdf', width = 6, height = 4)
a<-VlnPlot(Immu_Aging_V3_clean, features = c( "Chemo1"), split.by = 'Treatment',
           group.by = 'Celltype',pt.size = 0, col=c("#0072B2", "#D55E00"))+
  scale_x_discrete(limits=c('B cells', 'Macrophages','CD8+ T cells','CD4+ T cells',
                            'Plasma cells','Nonclassic','MAIT cells','Tregs',
                            'NKT cells'))+
  scale_y_continuous(limits =c(-0.25,1), breaks = c(0, 0.25, 0.5,0.75,1))
a
dev.off()

Test2<-a[[1]][["data"]]
Test3<-separate( Test2, ident, sep = '_', into = c('Celltype', 'Treatment') )
res2<-Test2  %>%
  group_by(ident) %>%
  wilcox_test(data =., Chemo1~split ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
res2
res

## telomere genes in epithelium and immune
telomere<-read.csv2('~/Documents/Documents - Joana’s MacBook Pro/IKMB/scRNAseq/Reference/Telomeres_Genes_Mmusculus.csv')
table(telomere$CellularFunction)
table(telomere$TM.classification)

Class<-subset(telomere, telomere$TM.classification %in% c('Telomerase','Telomerase & ALT', 'ALT'))
#Func<-subset(telomere, telomere$CellularFunction %in% c('Telomere biology'))

# Calculation of gene module for literature genes related to telomere (remove of Terc gene)
telo<-c("Terf1",   "Orc3",    "Klf4",    "Orc1",    "Lrwd1",   "Pot1a",   "Nr2c2",   "Nr2c1",   "Terf2",   "Terf2ip", "Tep1",    "Tinf2",   "Dnmt1",  
        "Atr",     "Ctnnb1",  "Aurkb",   "Tert",    "Pot1b",   "Tcof1"  )

DotPlot(Epi_Ageing_V3_clean, features = telo, assay = "RNA",group.by = 'Celltype_Treatment' )+
  scale_color_gradientn(colours  = rev(color)) + coord_flip() + scale_y_discrete(position = "right")+
  theme( axis.text.x = element_text(angle = 60, hjust = 0))


DefaultAssay(Epi_Ageing_V3_clean)<-'RNA'
Idents(Epi_Ageing_V3_clean)<-Epi_Ageing_V3_clean$Celltype_Treatment
Epi_Ageing_V3_clean <- AddModuleScore(Epi_Ageing_V3_clean, features =list(telo), name="Telomere", assay = 'RNA')

pdf('Vln_Epithelium_TelomereScore_Celtype_Treatment.pdf', width = 10, height = 5)
Tel<-VlnPlot(Epi_Ageing_V3_clean, features = c( "Telomere1"), split.by = 'Treatment',
               group.by = 'Celltype',pt.size = 0, col=c("#0072B2", "#D55E00"))+
  scale_x_discrete(limits=c('Enterocytes', 'Enterocytes Proximal','Enterocytes Distal',
                            'Immature Enterocytes','Stem cells (Lgr5+)','TA cells',
                            'Paneth cells','Goblet cells','Tuft cells','Enteroendocrine'))
Tel
dev.off()

Test<-Tel[[1]][["data"]]
res<-Test  %>%
  group_by(ident) %>%
  wilcox_test(data =., Telomere1~split ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
res


Immu_Ageing_V3_clean <- AddModuleScore(Immu_Ageing_V3_clean, features =list(telo), name="Telomere", assay = 'RNA')

pdf('Vln_Immune_TelomereScore_Celtype_Treatment.pdf', width = 10, height = 5)
Tel<-VlnPlot(Immu_Ageing_V3_clean, features = c( "Telomere1"), split.by = 'Treatment',
             group.by = 'Celltype',pt.size = 0, col=c("#0072B2", "#D55E00"))+
  scale_x_discrete(limits=c('B cells', 'Macrophages','CD8+ T cells','CD4+ T cells',
                            'Plasma cells','Nonclassic','MAIT cells','Tregs',
                            'NKT cells'))
Tel
dev.off()

Test<-Tel[[1]][["data"]]
res<-Test  %>%
  group_by(ident) %>%
  wilcox_test(data =., Telomere1~split ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
res


# Reference of Tabula Muris senis
## load Tabula muris senis object
# large intestine cells
Drop<-readRDS('c7c2124e-5d10-447e-89b9-7c5613f752c8.rds')
Facs<-readRDS('2e42118d-e35e-47b0-91cb-a11a98b2fa82.rds')

#spleen cells
Drop<-readRDS('7d090508-a291-4cc1-8b31-50b99bc5c670.rds')
Facs<-readRDS('244cc163-0376-4d58-aa96-7e61cb2fc9f4.rds')

Drop$data <- 'droplet'
Facs$data <- 'Facs'

Ref <- merge(x = Drop, y =Facs)

PEM_sep<-SplitObject(Ref, split.by = 'data')
for(i in 1:length(PEM_sep)) {
  PEM_sep[[i]] <- NormalizeData(PEM_sep[[i]], verbose = F)
  PEM_sep[[i]] <- FindVariableFeatures(PEM_sep[[i]], selection.method = 'vst', verbose = F)
}
ref_PEM<-PEM_sep[c('droplet','Facs')]
PEM_V3_anchors<-FindIntegrationAnchors(object.list = ref_PEM, dims = 1:30)
Ref_V3<-IntegrateData(anchorset = PEM_V3_anchors, dims=1:30)

DefaultAssay(Ref_V3)<-'RNA'
all.genes<-rownames(Ref_V3)
Ref_V3 <- ScaleData(Ref_V3, features = all.genes)

Ref_V3 <- FindVariableFeatures(Ref_V3)
Ref_V3 <- RunPCA(Ref_V3, npcs = 80)
ElbowPlot(Ref_V3, ndims = 80)

Ref_V3 <- FindNeighbors(Ref_V3,  dims = 1:40)
Ref_V3 <- FindClusters(Ref_V3, resolution = 0.1)
Ref_V3 <- RunUMAP(Ref_V3, dims = 1:40)


Idents(Ref_V3)<-Ref_V3$seurat_clusters
Stem<-c("Lgr5","Olfm4","Ascl2", "Slc12a1", "Gkn3")
Enterocytes<-c("Alpi","Fabp1","Apoa1", "Apoa4")
Goblet<-c("Muc2","Agr2","Clca3a1", "Tff3", "Clca3a2")
Paneth<-c("Lyz1","Defa17","Defa22", "Defa24", "Ang4")
EE<- c("Chga","Chgb","Tac1", "Tph1")
Tuft<-c("Dclk1","Trpm5","Gfi1b", "Il25")
TA<-c("Stmn1","Tubb5")
PE<-c('Sox4', 'Neurog3', 'Neurod2')
features = as.data.frame(c(Stem,Enterocytes,Goblet, Paneth, EE, Tuft, TA, PE))
markers_entrez <- bitr(features$`c(Stem, Enterocytes, Goblet, Paneth, EE, Tuft, TA, PE)`,
                       fromType = "SYMBOL",
                       toType="ENSEMBL",
                       OrgDb=org.Mm.eg.db)
features$SYMBOL<-features$`c(Stem, Enterocytes, Goblet, Paneth, EE, Tuft, TA, PE)`
features<-merge(features, markers_entrez,by.x='SYMBOL')
features$`c(Stem, Enterocytes, Goblet, Paneth, EE, Tuft, TA, PE)`=NULL

DotPlot(Ref_V3, features = features$ENSEMBL,dot.scale = 7,assay='RNA')+ RotatedAxis()

EnteroMatProx<-c("Apoa4","Fabp1","Apoc2", "Rbp2", "Apoc3", 'Leap2')
EnteroImatProx<-c("Casp6")
EnteroMatDistal<- c("Tmigd1","Fabp6","Slc51a", "Slc51b", "Mep1a", 'Fam151a')
EnteroImatDistal<-c("Reg3g","Gsdmc4","Prss32", "Krt8")
features = as.data.frame(c(EnteroMatProx,EnteroImatProx,EnteroMatDistal, EnteroImatDistal))
markers_entrez <- bitr(features$`c(EnteroMatProx, EnteroImatProx, EnteroMatDistal, EnteroImatDistal)`,
                       fromType = "SYMBOL",
                       toType="ENSEMBL",
                       OrgDb=org.Mm.eg.db)
features$SYMBOL<-features$`c(EnteroMatProx, EnteroImatProx, EnteroMatDistal, EnteroImatDistal)`
features<-merge(features, markers_entrez,by.x='SYMBOL')
features$`c(EnteroMatProx, EnteroImatProx, EnteroMatDistal, EnteroImatDistal)`=NULL
DotPlot(Ref_V3, features = features$ENSEMBL,dot.scale = 7,assay='RNA')+ RotatedAxis()

# annotation
Idents(Ref_V3)<-Ref_V3$seurat_clusters
new.cluster.ids <- c( 'Imma Enterocytes','Enterocytes','Imma Enterocytes','Imma Enterocytes',
                      'TA cells','TA cells','Goblet cells','Enterocytes',
                      "Enteroendocrine","Imma Enterocytes")
names(new.cluster.ids) <- levels(Ref_V3)
Ref_V3 <- RenameIdents(Ref_V3, new.cluster.ids)
Ref_V3$Celltype<-Idents(Ref_V3)

DimPlot(Ref_V3,  reduction = 'umap', group.by = 'Celltype')

# Finish annotation
Idents(Ref_V3)<-Ref_V3$free_annotation
new.cluster.ids <- c( 'Lgr5+ Stem cell','Enteroendocrine','Goblet cell','TA cell',
                      'Reg4/cKit+ deep crypt Distal_yMB','Enterocyte','Lgr5+ Stem cell','Enterocyte',
                      "Lgr5+ Stem cell","Goblet cell",'TA cell')
names(new.cluster.ids) <- levels(Ref_V3)
Ref_V3 <- RenameIdents(Ref_V3, new.cluster.ids)
Ref_V3$free_annotation2<-Idents(Ref_V3)

DimPlot(Ref_V3,  reduction = 'umap', group.by = 'free_annotation')

save(Ref_V3, file ='Ref_V3_Spleen_object_analyzed_droplet_facs.rda')
save(Ref_V3, file ='Ref_V3_LargeInt_object_analyzed_droplet_facs.rda')

#use only 3m and 30m
Idents(Ref_V3)<-Ref_V3$age
Age_Ref<-subset(Ref_V3, idents=c('3m','30m'))


# Calculation of gene module for literature genes related to EMT in large intestina
Strogenes<-as.data.frame(c('Vim','Ctnnb1','Fn1','Aifm2','Tgfb1', 'Tgfbr1','Smad2','Smad3','Smad4'))

markers_entrez <- bitr(Strogenes$`c("Vim", "Ctnnb1", "Fn1", "Aifm2", "Tgfb1", "Tgfbr1", "Smad2", "Smad3", "Smad4")`,
                       fromType = "SYMBOL",
                       toType="ENSEMBL",
                       OrgDb=org.Mm.eg.db)
Strogenes$SYMBOL<-Strogenes$`c("Vim", "Ctnnb1", "Fn1", "Aifm2", "Tgfb1", "Tgfbr1", "Smad2", "Smad3", "Smad4")`
Strogenes<-merge(Strogenes, markers_entrez,by.x='SYMBOL')
Strogenes$`c("Vim", "Ctnnb1", "Fn1", "Aifm2", "Tgfb1", "Tgfbr1", "Smad2", "Smad3", "Smad4")`=NULL

Ref_V3 <- AddModuleScore(Age_Ref, features = list(Strogenes$ENSEMBL), name="Stromal", assay = 'RNA')
pdf('Vln_Ref_V3_freeanoo_stromalScore_Celtyype_Treatment.pdf', width = 4, height = 3.5)
VlnPlot(Age_Ref, features = c( "Stromal1"), split.by = 'age',assay = 'RNA',
        group.by = 'free_annotation2',pt.size = 0, col=c("#0072B2", "lightblue"))+
  scale_x_discrete(limits=rev(c("Lgr5+ Stem cell" ,'TA cell',"Goblet cell" ,"Enterocyte")))
dev.off()

Test<-Strom[[1]][["data"]]
Test2<-subset(Test, !Test$ident %in% c('Reg4/cKit+ deep crypt Distal_yMB', 'Enteroendocrine'))
res<-Test2  %>%
  group_by(ident) %>%
  wilcox_test(data =., Stromal1~split ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
res

Strom<-VlnPlot(Ref_V3, features = c( "Stromal1"), group.by = 'age', pt.size = 0, assay = 'RNA')
Test<-Strom[[1]][["data"]]
res<-Test  %>%
  wilcox_test(data =., Stromal1~ident ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") 
res
stat.test <- stat.test %>% add_xy_position(x = "group")

mean2<-aggregate(Stromal1 ~ ident, Test, mean)

# Calculation of gene module for literature genes related to EMT in spleen
Chemo<-as.data.frame(c('Tnf', 'Ifng', 'Il1b','Il2','Il6','Cxcl15', 'Ccl20', 'Ccl9','Ccr1','Nfkb1','Myd88','Tlr6'))

markers_entrez <- bitr(Chemo$`c("Tnf", "Ifng", "Il1b", "Il2", "Il6", "Cxcl15", "Ccl20", "Ccl9", "Ccr1", "Nfkb1", "Myd88", "Tlr6")`,
                       fromType = "SYMBOL",
                       toType="ENSEMBL",
                       OrgDb=org.Mm.eg.db)
Chemo$SYMBOL<-Chemo$`c("Tnf", "Ifng", "Il1b", "Il2", "Il6", "Cxcl15", "Ccl20", "Ccl9", "Ccr1", "Nfkb1", "Myd88", "Tlr6")`
Chemo<-merge(Chemo, markers_entrez,by.x='SYMBOL')
Chemo$`c("Tnf", "Ifng", "Il1b", "Il2", "Il6", "Cxcl15", "Ccl20", "Ccl9", "Ccr1", "Nfkb1", "Myd88", "Tlr6")`=NULL

Age_Ref <- AddModuleScore(Age_Ref, features = list(Ccgenes$ENSEMBL), name="Chemo", assay = 'RNA')
pdf('Vln_Ref_V3_Spleen_cell_types_ChemoScore_Celtype_Treatment.pdf', width = 4, height = 4)
a<-VlnPlot(Age_Ref, features = c( "Chemo1"), split.by = 'age',assay = 'RNA',
           group.by = 'cell_type',pt.size = 0, col=c("#0072B2", "lightblue"))+
  scale_x_discrete(limits=c("B cell" ,'macrophage','T cell','plasma cell',
                            "mature NK T cell" ))+
  scale_y_continuous(limits =c(-0.25,1), breaks = c(0, 0.25, 0.5,0.75,1))
a
dev.off()

Test2<-a[[1]][["data"]]
Test3<-subset( Test2, !Test2$ident %in% c('CD4-positive, alpha-beta T cell',
                                          'CD8-positive, alpha-beta T cell'))

Test3<-separate( Test2, ident, sep = '_', into = c('cell_type', 'age') )
res2<-Test3  %>%
  group_by(ident) %>%
  wilcox_test(data =., Chemo1~split ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
res2
res


