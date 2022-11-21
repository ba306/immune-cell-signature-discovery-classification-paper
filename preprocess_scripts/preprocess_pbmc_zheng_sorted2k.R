library(Seurat)
library(varhandle)
#zheng sorted pbmc downloaded from 10X
set.seed(42)

dataset_dir="N:/datasets/PBMC/zheng_sorted/"
dir_klist=list.dirs(dataset_dir)[-1]

cell_names=gsub( "N:/datasets/PBMC/zheng_sorted//","",dir_klist)
cell_Types=c()

library(foreach)

mtx_all=foreach(i=seq_along(cell_names),.combine = cbind)%do%{
  print(i)
  s=Read10X(data.dir = dir_klist[i])
  cell_Types=c(cell_Types,rep(cell_names[i],ncol(s)))
  s
}
cell_Types_df=data.frame(cell_Types)

dim(mtx_all)
#32738 63947

rownames(cell_Types_df)=make.unique(colnames(mtx_all))
colnames(mtx_all)=rownames(cell_Types_df)

sobj=CreateSeuratObject(counts = mtx_all,meta.data = cell_Types_df, min.cells = 3, min.features = 200)
dim(sobj)

sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sobj, features = c("nFeature_RNA"), pt.size = 0,y.max = 2000)
VlnPlot(sobj, features = c("percent.mt"), pt.size = 0,y.max = 20)

sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 5)

sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

table(sobj$cell_Types)
sobj$cell_Types=gsub("cd14_monoytes","Monocytes_CD14",sobj$cell_Types,fixed = T)
sobj$cell_Types=gsub("cytotoxic_tcd8","TCD8_cytotoxic",sobj$cell_Types,fixed = T)
sobj$cell_Types=gsub("naive_cytotoxic","TCD8_naive",sobj$cell_Types,fixed = T)
sobj$cell_Types=gsub("nk","NK",sobj$cell_Types,fixed = T)

saveRDS(sobj,file = paste0(dataset_dir,"pbmc_zheng_sorted_sobj.RDS"))

#2k cell per cell type 

sobj=readRDS(file = paste0(dataset_dir,"pbmc_zheng_sorted_sobj.RDS"))
Idents(sobj)=sobj$cell_Types
sobj_2k <- subset(sobj, downsample = 2000)
table(sobj_2k$cell_Types)
saveRDS(sobj_2k,file = paste0(dataset_dir,"pbmc_zheng_sorted_sobj2k.RDS"))
