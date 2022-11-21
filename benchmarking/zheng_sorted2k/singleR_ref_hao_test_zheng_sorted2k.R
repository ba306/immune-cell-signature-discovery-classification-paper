library(Seurat)
library(dplyr)
set.seed(42)

source("N:./immunsig_scripts_manuscript/benchmarking/celltype_annot_stats.R")
# Script for cell type annotation/prediction using singleR package

# reference zilionis tumor
dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/zheng_sorted2k/"


# Load reference data
reference_seurat <- readRDS("Z:/pbmc_hao_ref.RDS")

# Load query data
# Load input scRNAseq data 
tiss_immune <- readRDS("N:/datasets/PBMC/zheng_sorted/pbmc_zheng_sorted_sobj2k.RDS")

# normalized query assay loaded
test_assay <- as.SingleCellExperiment(tiss_immune)
# create single cell experiment for reference data
ref_sce <- as.SingleCellExperiment(reference_seurat)

# rm(tiss_immune)
rm(reference_seurat)
# prediction
library(SingleR)
start.time <- Sys.time()

pred.hesc <- SingleR(test = test_assay, ref = ref_sce,   assay.type.test = "logcounts",
                     assay.type.ref = "logcounts",
                     labels = ref_sce$harmonized_celltype,fine.tune = F)
# save prediction
saveRDS(pred.hesc,file = paste0(dir_analysis,"singleR_prediction_ref_hao_zhengsorted2k.RDS"))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken #2.3 h

pred.hesc=readRDS(file = paste0(dir_analysis,"singleR_prediction_ref_hao_zhengsorted2k.RDS"))

tiss_immune$K1=tiss_immune$cell_Types
tiss_immune$K1=gsub("Monocytes_CD14","Mono",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD4_memory","T CD4 memory_naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD4_naive","T CD4 memory_naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD8_cytotoxic","T CD8",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD8_naive","T CD8",tiss_immune$K1,fixed = T)


#accuracy harmonize

df_filter= data.frame(original= tiss_immune$K1,predicted=pred.hesc$pruned.labels)
df_filter$original=unfactor(df_filter$original)
df_filter$predicted=unfactor(df_filter$predicted)
table(df_filter)
Cell_type_stats(df_filter)
# Sensitivity   PPV NPV Specificity Accuracy f1_score
# 1       81.71 90.48 94.69       94.61    91.61    85.87
