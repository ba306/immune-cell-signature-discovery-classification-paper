#Seurat 180 HVGs

library(Seurat)
library(dplyr)
library(foreach)
source("N:./immunsig_scripts_manuscript/benchmarking/celltype_annot_stats.R")
set.seed(42)

# Script for cell type annotation/prediction using Seurat package
dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/zheng_sorted2k/"


# Load query data
# Load input scRNAseq data 
tiss_immune <- readRDS("N:/datasets/PBMC/zheng_sorted/pbmc_zheng_sorted_sobj2k.RDS")

# Load reference data
reference_seurat <- readRDS("Z:/pbmc_hao_ref.RDS")

##prediction###
  reference_seurat <- FindVariableFeatures(object = reference_seurat, 
                                           selection.method = 'vst', nfeatures = 180)
  reference_seurat <- ScaleData(reference_seurat,features = VariableFeatures(reference_seurat))

  start.time <- Sys.time()
  
  
  sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = tiss_immune,
                                     dims = 1:30)
  predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$harmonized_celltype,
                              dims = 1:30)
  #query_seurat <- AddMetaData(object = tiss_immune, metadata = predictions)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken %>% print()

saveRDS(predictions,file = paste0(dir_analysis,"seurat_hao_zhengsorted180HVGs.RDS"))

predictions=readRDS(file = paste0(dir_analysis,"seurat_hao_zhengsorted180HVGs.RDS"))

tiss_immune$K1=tiss_immune$cell_Types
tiss_immune$K1=gsub("Monocytes_CD14","Mono",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD4_memory","T CD4 memory_naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD4_naive","T CD4 memory_naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD8_cytotoxic","T CD8",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD8_naive","T CD8",tiss_immune$K1,fixed = T)

  df_filter= data.frame(original= tiss_immune$K1,predicted=predictions$predicted.id)
  df_filter$original=unfactor(df_filter$original)
  df_filter$predicted=unfactor(df_filter$predicted)
table(df_filter)
  Cell_type_stats(df_filter)
  # Sensitivity   PPV NPV Specificity Accuracy f1_score
  # 1       77.38 81.44  93       91.54    88.09    79.36