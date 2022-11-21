source("N:./immunsig_scripts_manuscript/benchmarking/celltype_annot_stats.R")
library(dplyr)
library(varhandle)

general_harmonnize_pred_stat=function(new_clusters,sobj){
  original=sobj$K1
  original[sobj$K3=="DC"]="DC"
  
  original=gsub("pDC","DC",original,fixed = T)
  original=gsub("Monocyte CD14 / DC","Mono",original,fixed = T)
  original=gsub("Monocyte CD16","Mono",original,fixed = T)
  original=gsub("T CD4 memory","TCD4 memory_naive",original,fixed = T)
  original=gsub("T CD4 naive","TCD4 memory_naive",original,fixed = T)
  original=gsub("T CD8 memory","TCD8",original,fixed = T)
  original=gsub("T CD8 naive","TCD8",original,fixed = T)
  
  new_clusters=gsub("T CD8","TCD8",new_clusters,fixed = T)
  new_clusters=gsub("T CD4 memory_naive","TCD4 memory_naive",new_clusters,fixed = T)
  new_clusters=gsub("T CD4 naive_memory","TCD4 memory_naive",new_clusters,fixed = T)
  
  df_filter= data.frame(original=original,predicted=new_clusters)
  df_filter$original=unfactor(df_filter$original)
  df_filter$predicted=unfactor(df_filter$predicted)
  print(sort(unique(df_filter$original)))
        print(sort(unique(df_filter$predicted)))
              
    data.frame(Cell_type_stats(df_filter))
}
