library(Seurat)
library(dplyr)
library(foreach)
library(ggplot2)
library(reshape2)
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
pred_list=foreach(i=seq(100,2000,100)) %do%{
  print(i)
  reference_seurat <- FindVariableFeatures(object = reference_seurat, 
                                           selection.method = 'vst', nfeatures = i)
  reference_seurat <- ScaleData(reference_seurat,features = VariableFeatures(reference_seurat))
  #tiss_immune <- ScaleData(tiss_immune,features = VariableFeatures(reference_seurat))
  
  start.time <- Sys.time()
  
  
  sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = tiss_immune,
                                     dims = 1:30)
  predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$harmonized_celltype,
                              dims = 1:30)
  #query_seurat <- AddMetaData(object = tiss_immune, metadata = predictions)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken %>% print()
  # saveRDS(predictions,file = paste0(dir_analysis,as.character(i),"_seurat_hao_zhengsorted.RDS"))
  
  predictions
}
saveRDS(pred_list,file = paste0(dir_analysis,"prediction_list_seurat_hao_zhengsorted.RDS"))

pred_list=readRDS(file = paste0(dir_analysis,"prediction_list_seurat_hao_zhengsorted.RDS"))
names(pred_list)=seq(100,2000,100)

tiss_immune$K1=tiss_immune$cell_Types
tiss_immune$K1=gsub("Monocytes_CD14","Mono",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD4_memory","T CD4 memory_naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD4_naive","T CD4 memory_naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD8_cytotoxic","T CD8",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD8_naive","T CD8",tiss_immune$K1,fixed = T)


stats_df=foreach(i=seq_along(pred_list),.combine=rbind) %do%{
  
  df_filter= data.frame(original= tiss_immune$K1,predicted=pred_list[[i]]$predicted.id)
  df_filter$original=unfactor(df_filter$original)
  df_filter$predicted=unfactor(df_filter$predicted)
  
  stats=Cell_type_stats(df_filter)
  data.frame(noHVGs=as.character(names(pred_list)[i]),stats)
  
}


stats_df_melt=melt(stats_df)
stats_df_melt$noHVGs=as.numeric(unfactor(stats_df_melt$noHVGs))
stats_df_melt$variable=gsub("f1_score","F1 score",stats_df_melt$variable)


# svg(paste0(dir_analysis,"Supp_seurat_diffHVGs_zhengsorted_ref_hao.svg"),width = 6,height = 4)
plot=ggplot(stats_df_melt, aes(noHVGs, value,colour=variable)) + geom_point() + geom_line()+ 
  scale_x_continuous(breaks = seq(0, 2000, by = 400))+
  xlab("Number of HVGs") + ylab("Metrics score")+ labs(colour = "Metrics")+theme_bw(base_size=22)+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))+ theme(legend.position="top",
                                                      legend.margin=margin(0,0,0,0),
                                                      legend.box.margin=margin(-10,-10,-10,-10))
# plot
# dev.off()

saveRDS(plot,paste0(dir_analysis,
                    "Supp_seurat_diffHVGs_zhengsorted_ref_hao.RDS"))




