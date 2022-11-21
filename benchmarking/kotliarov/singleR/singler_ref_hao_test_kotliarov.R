############################################################
# 	Author: Bogac Aybey                                    #
#             Merck KGaA - DKFZ                            #
############################################################
library(Seurat)
library(SeuratDisk)
library(dplyr)
set.seed(42)

source("N:./immunsig_scripts_manuscript/benchmarking/kotliarov/kotliarov_harmonize_prediction_scores_func.R")
# Script for cell type annotation/prediction using singleR package

# reference Hao tumor
dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/kotliarov/singleR/"


# Load reference data
reference_seurat <- readRDS("Z:/pbmc_hao_ref.RDS")

# Load query data
# Load input scRNAseq data 
tiss_immune <- readRDS("N:/datasets/PBMC/kolitrov/kotliarov_pbmc_lognorm_sobj2.RDS")


# normalized query assay loaded
test_assay <- as.SingleCellExperiment(tiss_immune)
# create single cell experiment for reference data
ref_sce <- as.SingleCellExperiment(reference_seurat)

rm(tiss_immune)
rm(reference_seurat)
# prediction
library(SingleR)
start.time <- Sys.time()

pred.hesc <- SingleR(test = test_assay, ref = ref_sce,   assay.type.test = "logcounts",
                     assay.type.ref = "logcounts",
                     labels = ref_sce$harmonized_celltype,fine.tune = F)
# save prediction
saveRDS(pred.hesc,file = paste0(dir_analysis,"singleR_prediction_ref_hao_kotliarov.RDS"))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken # 4.8 h

pred.hesc=readRDS(file = paste0(dir_analysis,"singleR_prediction_ref_hao_kotliarov.RDS"))


#accuracy harmonize

general_harmonnize_pred_stat(new_clusters=pred.hesc$pruned.labels,test_assay)


# Sensitivity  PPV   NPV Specificity Accuracy f1_score
# 1       88.18 89.2 96.24       94.96    94.09    88.69
