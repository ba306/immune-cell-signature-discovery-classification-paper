# Preparing reference dataset: Hao PBMC
# downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378

library(Seurat)
library(dplyr)
library(qs)

dir_name="./data/Hao/"

# load data
sobj=Seurat::Read10X(data.dir =dir_name)

meta=read.csv(paste0(dir_name,"GSE164378_sc.meta.data_3P.csv"),row.names = 1)

# create sobj
hao=CreateSeuratObject(counts = sobj,meta.data = meta)

hao[["percent.mt"]] <- PercentageFeatureSet(hao, pattern = "^MT-")
VlnPlot(tiss_immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(hao$nFeature_RNA)#min 501, max 5990 
summary(hao$percent.mt) #0-15


# Remove doublet, erythrocytes and platelets 

remove=
  hao$celltype.l2 !="Doublet"&
  hao$celltype.l2 !="Eryth" &
  hao$celltype.l2 !="Platelet"

hao=hao[,remove]

hao=NormalizeData(hao)


# Harmonize cell types

hao$harmonized_celltype=hao$celltype.l2
hao$harmonized_celltype=gsub("B intermediate","B",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("B memory","B",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("B naive","B",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("Plasmablast","B",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("cDC1","DC",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("cDC2","DC",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD14 Mono","Mono",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD16 Mono","Mono",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("NK_CD56bright","NK",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("NK Proliferating","NK",hao$harmonized_celltype)

hao$harmonized_celltype=gsub("CD4 TCM","T CD4",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD4 TEM","T CD4",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD4 Naive","T CD4",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD4 CTL","T CD4",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD4 Proliferating","T CD4",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("Treg","T CD4",hao$harmonized_celltype)

hao$harmonized_celltype=gsub("gdT","Other cells",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("dnT","Other cells",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("MAIT","Other cells",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("ILC","Other cells",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("HSPC","Other cells",hao$harmonized_celltype)

hao$harmonized_celltype=gsub("CD8 TCM","T CD8",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD8 TEM","T CD8",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD8 Naive","T CD8",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD8 Proliferating","T CD8",hao$harmonized_celltype)

hao$harmonized_celltype=gsub("ASDC","DC",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("pDC","DC",hao$harmonized_celltype)


dim(hao)
#cells vs genes
#158783 33538 

# average number of genes per cell
mean(hao@meta.data$nFeature_RNA)
#2206.678

# save qs for referencing purposes
qsave(hao,paste0(dir_name,"pbmc_hao_ref_up.qs"))
