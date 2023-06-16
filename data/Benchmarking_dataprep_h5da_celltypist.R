library(anndata)
library(Seurat)
library(qs)

# Preparing benchmarking datasets for celltypist format

# Convert seurat to anndata
# only raw counts 

########### Hao PBMC reference ########### 
reference_seurat <- qread("./data/Hao/pbmc_hao_ref_up.qs")
reference_seurat=reference_seurat[,reference_seurat$harmonized_celltype!="Other cells"]

mat=reference_seurat@assays$RNA@data
mat=t(as.matrix(mat))
ann=anndata::AnnData(mat)
ann$obs$harmonized_celltype=reference_seurat$harmonized_celltype
write_h5ad(anndata = ann,filename = "./data/Hao/pbmc_hao_ref_up.h5ad")

###########  Zheng sorted FACS 2k per cell type- benchmarking  ########### 
reference_zheng <- qread("./data/Zheng/pbmc_zheng_sorted_2k.qs")

mat=reference_zheng@assays$RNA@data
mat=t(as.matrix(mat))
ann=anndata::AnnData(mat)
write_h5ad(anndata = ann,filename = "./data/Zheng/pbmc_zheng_sorted_2k.h5ad")

###########  Kotliaraov PBMC - benchmarking  ########### 
kotliarov <- qread("./data/Kotliarov/kotliarov_pbmc.qs")

mat=kotliarov@assays$RNA@data
mat=t(as.matrix(mat))
ann=anndata::AnnData(mat)
write_h5ad(anndata = ann,filename = "./data/Kotliarov/kotliarov_pbmc.h5ad")

########### Kartha PBMC interferon ########### 
kartha <- qs::qread("./data/Kartha/Kartha_sobj.qs")

mat=kartha@assays$RNA@data
mat=t(as.matrix(mat))
ann=anndata::AnnData(mat)
write_h5ad(anndata = ann,filename = "./data/Kartha/Kartha.h5ad")
