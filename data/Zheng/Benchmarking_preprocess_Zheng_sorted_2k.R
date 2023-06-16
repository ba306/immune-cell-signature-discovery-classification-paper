# Preparing benchmarking dataset: Zheng FACS sorted 10 X PBMC
# Downloaded from 10X and saved data from each cell type --> put into different folders for each cell type

library(Seurat)
library(varhandle)
library(foreach)
library(qs)

set.seed(42)

dataset_dir="./data/Zheng/"

# list all files containing matrices and gene names for each cell
dir_klist=list.dirs(dataset_dir)[-1]

cell_names=gsub( "./data/Zheng//","",dir_klist)
cell_Types=c()

# compile data from each cell type
mtx_all=foreach(i=seq_along(cell_names),.combine = cbind)%do%{
  print(i)
  sobj_celltypes=Read10X(data.dir = dir_klist[i])
  cell_Types=c(cell_Types,rep(cell_names[i],ncol(sobj_celltypes)))
  sobj_celltypes
}
cell_Types_df=data.frame(cell_Types)

dim(mtx_all)
#32738 85423

rownames(cell_Types_df)=make.unique(colnames(mtx_all))
colnames(mtx_all)=rownames(cell_Types_df)

# create sobj
sobj=CreateSeuratObject(counts = mtx_all,meta.data = cell_Types_df, min.cells = 3, min.features = 200)
dim(sobj)
#18171 85394
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sobj, features = c("nFeature_RNA"), pt.size = 0,y.max = 2000)
VlnPlot(sobj, features = c("percent.mt"), pt.size = 0,y.max = 20)

sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & 
                 percent.mt < 5)

sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

sobj$cell_Types=gsub("Thelper","TCD4_helper",sobj$cell_Types,fixed = T)
sobj$cell_Types=gsub("Treg","TCD4_reg",sobj$cell_Types,fixed = T)
sobj$cell_Types=gsub("cd14_monoytes","Monocytes_CD14",sobj$cell_Types,fixed = T)
sobj$cell_Types=gsub("cytotoxic_tcd8","TCD8_cytotoxic",sobj$cell_Types,fixed = T)
sobj$cell_Types=gsub("naive_cytotoxic","TCD8_naive",sobj$cell_Types,fixed = T)
sobj$cell_Types=gsub("nk","NK",sobj$cell_Types,fixed = T)

dim(sobj)
#18171 84738
table(sobj$cell_Types)

# 2k samples from each cell for benchmarking

Idents(sobj)=sobj$cell_Types
sobj <- subset(sobj, downsample = 2000)

dim(sobj)
#[1] 18171 18000

# Harmonize cell types

sobj$harmonized_celltype=sobj$cell_Types
sobj$harmonized_celltype=gsub("Monocytes_CD14","Mono",sobj$harmonized_celltype,fixed = T)
sobj$harmonized_celltype=gsub("TCD4_memory","T CD4",sobj$harmonized_celltype,fixed = T)
sobj$harmonized_celltype=gsub("TCD4_naive","T CD4",sobj$harmonized_celltype,fixed = T)
sobj$harmonized_celltype=gsub("TCD4_reg","T CD4",sobj$harmonized_celltype,fixed = T)
sobj$harmonized_celltype=gsub("TCD4_helper","T CD4",sobj$harmonized_celltype,fixed = T)

sobj$harmonized_celltype=gsub("TCD8_cytotoxic","T CD8",sobj$harmonized_celltype,fixed = T)
sobj$harmonized_celltype=gsub("TCD8_naive","T CD8",sobj$harmonized_celltype,fixed = T)

# average number of genes per cell
mean(sobj@meta.data$nFeature_RNA)
#557.2158

qsave(sobj,file = paste0(dataset_dir,"pbmc_zheng_sorted_2k.qs"))

