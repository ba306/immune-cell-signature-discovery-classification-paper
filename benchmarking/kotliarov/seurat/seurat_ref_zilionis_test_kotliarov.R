library(Seurat)
library(dplyr)
source("N:./immunsig_scripts_manuscript/benchmarking/kotliarov/kotliarov_harmonize_prediction_scores_func.R")

set.seed(42)

# Script for cell type annotation/prediction using Seurat package
# reference zilionis tumor
dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/kotliarov/seurat/"

# Load query data
# Load input scRNAseq data 
tiss_immune <- readRDS("N:/datasets/PBMC/kolitrov/kotliarov_pbmc_lognorm_sobj2.RDS")

# Load reference data
# Lung scRNAseq data from Zilionis et al. 2019
reference_seurat <- readRDS("N:/datasets/myeloid/myeloid_raw_tumor_filtered_sobj.RDS")
table(reference_seurat$cell.types)
reference_seurat=reference_seurat[,reference_seurat$cell.types!="Eosinophils"]
reference_seurat=reference_seurat[,reference_seurat$cell.types!="Neutrophils"]
reference_seurat=reference_seurat[,reference_seurat$cell.types!="Mast cells activated"]
reference_seurat=reference_seurat[,reference_seurat$cell.types!="Mast cells resting"]
reference_seurat=reference_seurat[,reference_seurat$cell.types!="Macrophages M0"]
reference_seurat=reference_seurat[,reference_seurat$cell.types!="Macrophages M1"]
reference_seurat=reference_seurat[,reference_seurat$cell.types!="Macrophages M2"]

reference_seurat=reference_seurat[,reference_seurat$cell.types !="T cells follicular helper" &
                                    reference_seurat$cell.types !="T cells regulatory (Tregs)" ]
reference_seurat$cell.types=gsub("B cells memory","B",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("B cells naive","B",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("Plasma cells","B",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("Dendritic cells activated","DC",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("Dendritic cells resting","DC",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("NK cells activated","NK",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("NK cells resting","NK",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("T cells CD4 naive","T CD4 naive_memory",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("T cells CD4 memory resting","T CD4 naive_memory",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("T cells CD4 memory activated","T CD4 naive_memory",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("T cells CD8","T CD8",reference_seurat$cell.types)
reference_seurat$cell.types=gsub("Monocytes","Mono",reference_seurat$cell.types)

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
predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$cell.types,
                            dims = 1:30)
query_seurat <- AddMetaData(object = tiss_immune, metadata = predictions)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#save predictions
save(predictions,file = paste0(dir_analysis,"seurat_zil_kotliarov.RDS"))

predictions=readRDS(file = paste0(dir_analysis,"seurat_zil_kotliarov.RDS"))
#accuracy harmonize
general_harmonnize_pred_stat(new_clusters=predictions$predicted.id,tiss_immune)
# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1       79.05 77.75 93.67       90.11    89.16    78.39