
library(dplyr)
library(varhandle)
library(foreach)
library("cluster")
library(Matrix)
library(matrixStats)
library(Seurat)
library(fst)
set.seed(42)

# Script for integration of 4 cancer scRNAseq discovery datasets
dirname="N:./integrate_4/v3/w_seurat_exp/"
save_dir="N:./immunsig_scripts_manuscript/immune_signature_discovery_workflow/"

##################################################################
############### Load and prepare datasets ########################
##################################################################

# Load datasets
li_seurat= readRDS(file = "N:./datasets/Li_et_al_Cell_2018_Melanoma_20190312.RDS")
melanoma_2_seurat= readRDS(file = "N:./datasets/Melanoma_2018_GSE115978_20190514.RDS")
puram= readRDS(file = "N:./datasets/puram/HNSCC_2_filtered_sobj.RDS")
hepato=readRDS(file = "N:./datasets/hepato/hepato_counts_filtered_sobj_2.RDS")

# save the datasets except hepato in a list
data_list= list(li_seurat,melanoma_2_seurat,puram)

# Calculate log2(x+1) for datasets except hepato 
# find 10k variable features individually
data_list= lapply(data_list, function(x){
  x@assays$RNA@data=Matrix(log2(x@assays$RNA@data+1),sparse = T)
  x <- FindVariableFeatures(x, selection.method = "vst", 
                            nfeatures = 10000, verbose = T)
  x
})

# normalize and find 10k variable features in hepato dataset
hepato <- NormalizeData(hepato, verbose = T)
hepato <- FindVariableFeatures(hepato, selection.method = "vst", 
                               nfeatures = 10000, verbose = T)

# add hepato to the data list 
data_list[[4]]= hepato
names(data_list)=c("li","melanoma_2","puram","hepato")
#rm(li_seurat,melanoma_2_seurat,puram,hepato)

#find common intgeration features (genes)
features <- SelectIntegrationFeatures(object.list = data_list,nfeatures = 10000)
save(features,file = paste0(dirname,"features.RData"))

# In order to use unified annotations, we adapted the annotations accordingly.

puram@meta.data[["non.cancer.cell.type"]]=unfactor(puram@meta.data[["non.cancer.cell.type"]])
puram@meta.data[["non.cancer.cell.type"]] <- replace(puram@meta.data[["non.cancer.cell.type"]], 
                                                     puram@meta.data[["non.cancer.cell.type"]] == "B cell", 
                                                     "B.cell")
puram@meta.data[["non.cancer.cell.type"]] <- replace(puram@meta.data[["non.cancer.cell.type"]], 
                                                     puram@meta.data[["non.cancer.cell.type"]] == "T cell", 
                                                     "T.cell")
puram@meta.data[["non.cancer.cell.type"]] <- replace(puram@meta.data[["non.cancer.cell.type"]], 
                                                     puram@meta.data[["non.cancer.cell.type"]] == "Dendritic", 
                                                     "DC")
puram@meta.data[["non.cancer.cell.type"]] <- replace(puram@meta.data[["non.cancer.cell.type"]], 
                                                     puram@meta.data[["non.cancer.cell.type"]] == "0", 
                                                     "NA")
puram@meta.data[["non.cancer.cell.type"]] <- replace(puram@meta.data[["non.cancer.cell.type"]], 
                                                     puram@meta.data[["non.cancer.cell.type"]] == "Endothelial", 
                                                     "Endo.")

cell_types=  c(li_seurat$cell_type ,melanoma_2_seurat$cell_type,puram$non.cancer.cell.type,
               hepato$celltype_sub) 

cell_types= replace(cell_types,cell_types=="Dysfunctional_CD8",
                    "T.CD8.dysfunctional")
cell_types= replace(cell_types,cell_types=="Treg",
                    "T.reg")
cell_types= replace(cell_types,cell_types=="Naive",
                    "T.naive")
cell_types= replace(cell_types,cell_types=="Cytotoxic_CD8_effector",
                    "T.CD8.cytotoxic_effector")
cell_types= replace(cell_types,cell_types=="Transitional_CD8_effector",
                    "T.CD8.transitional.effector")
cell_types= replace(cell_types,cell_types=="Memory_T_cell",
                    "T.memory")
cell_types= replace(cell_types,cell_types=="Dysfunctional_CD4",
                    "T.CD4.dysfunctional")
cell_types= replace(cell_types,cell_types=="non-classic-monocyte",
                    "monocyte.non.classic")
cell_types= replace(cell_types,cell_types=="immature-DC",
                    "DC.immature")
cell_types= replace(cell_types,cell_types=="mature-DC",
                    "DC.mature")
cell_types=replace(cell_types,cell_types=="MÏ???","Macrophage")

# save cell type annotation metadata 
# save(cell_types,file=paste0(dirname,"cell_types.RData"))
# load(file=paste0(dirname,"cell_types.RData"))

cell_types_list= foreach(i=seq_along(data_list))%do% {
  cell_types[colnames(data_list[[i]])]
}
save(cell_types_list,file=paste0(dirname,"cell_types_list.RData"))


# scale datasets individually
# run pca on 10k genes individually
data_list <- lapply(data_list, function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})

##################################################################
######################## Integration #############################
##################################################################

# Find anchors in 4 datasets
# reciprocal pca ('rpca') method was used 
# as reference 2nd dataset in the datalist was used
anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:30, 
                                  anchor.features = 10000,reference =2, reduction = "rpca")

# Integrate datasets using the anchors
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(integrated) <- "integrated"

# Scale integrated dataset
integrated <- ScaleData(integrated, verbose = T)

# Check integration on UMAP plot
integrated <- RunPCA(integrated, npcs = 30, verbose = T)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)

integrated@meta.data[["active.ident"]]= unfactor(integrated@active.ident)

pdf(paste0(save_dir,"Fig2a_umap_integration_celltypes_seurat.pdf"),width = 25,height = 18)
DimPlot(integrated, reduction = "umap", label = T, pt.size = .1,
        label.size = 10,group.by = "cell_types",repel = T)
dev.off()

# save scaled integrated expression matrix as fst
# cell types and gene names have been saved as cell_types.RData and features.RData,respectively
seurat_exp=  as.data.frame(integrated@assays[["integrated"]]@scale.data)
write.fst(seurat_exp,path = paste0(dirname,"seurat_exp.fst"))

#umap according to datasets
ident=c(rep("GSE123139",ncol(li_seurat)),
        rep("GSE115978",ncol(melanoma_2_seurat)),
        rep("GSE103322",ncol(puram)),
        rep("GSE140228",ncol(hepato)))
integrated@meta.data[["proj"]]=ident
pdf(paste0(save_dir,"Fig2b_umap_integration_datasets_seurat.pdf"),width = 20,height = 15)

DimPlot(integrated, reduction = "umap", label = F, pt.size = .1,group.by = "proj",repel = T)

dev.off()

