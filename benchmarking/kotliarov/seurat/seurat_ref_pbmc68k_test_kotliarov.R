library(Seurat)
library(dplyr)
source("N:./immunsig_scripts_manuscript/benchmarking/kotliarov/kotliarov_harmonize_prediction_scores_func.R")

set.seed(42)

# Script for cell type annotation/prediction using Seurat package
# reference pbmc68k
dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/kotliarov/seurat/"


# Load query data
# Load input scRNAseq data 
tiss_immune <- readRDS("N:/datasets/PBMC/kolitrov/kotliarov_pbmc_lognorm_sobj2.RDS")

# Load reference data
reference_seurat=readRDS("N:/datasets/PBMC/pbmc68k/PBMC_68k.rds")
reference_seurat=subset(reference_seurat,celltype!="CD34+")
reference_seurat=subset(reference_seurat,celltype!="CD4+/CD25 T Reg")
reference_seurat=subset(reference_seurat,celltype!="CD4+ T Helper2")

table(reference_seurat$celltype)

reference_seurat$celltype=gsub("CD19+ B","B",reference_seurat$celltype,fixed = T)
reference_seurat$celltype=gsub("CD14+ Monocyte","Mono",reference_seurat$celltype,fixed = T)
reference_seurat$celltype=gsub("CD4+/CD45RA+/CD25- Naive T","TCD4 memory_naive",reference_seurat$celltype,fixed = T)
reference_seurat$celltype=gsub("CD4+/CD45RO+ Memory","TCD4 memory_naive",reference_seurat$celltype,fixed = T)
reference_seurat$celltype=gsub("CD56+ NK","NK",reference_seurat$celltype,fixed = T)
reference_seurat$celltype=gsub("CD8+ Cytotoxic T","TCD8",reference_seurat$celltype,fixed = T)
reference_seurat$celltype=gsub("CD8+/CD45RA+ Naive Cytotoxic","TCD8",reference_seurat$celltype,fixed = T)
reference_seurat$celltype=gsub("Dendritic","DC",reference_seurat$celltype,fixed = T)

#standard pipeline
# normalize reference
# find 2k variable genes
# scale using 2k variable genes
reference_seurat <- NormalizeData(object = reference_seurat)
reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = 'vst', nfeatures = 2000)
reference_seurat <- ScaleData(reference_seurat)

##prediction###
start.time <- Sys.time()

sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = tiss_immune,
                                   dims = 1:30)
predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$celltype,
                            dims = 1:30)
query_seurat <- AddMetaData(object = tiss_immune, metadata = predictions)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#save predictions
saveRDS(predictions,file = paste0(dir_analysis,"seurat_pbmc68k_kotliarov.RDS"))

predictions=readRDS(file = paste0(dir_analysis,"seurat_pbmc68k_kotliarov.RDS"))

##accuracy
general_harmonnize_pred_stat(new_clusters=predictions$predicted.id,tiss_immune)
