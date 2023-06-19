library(ProfilerAPI2)
library("qs")
library(SeuratObject)
library(dplyr)
library(reshape2)
library(foreach)

data_dir="./data/TME/"
#### Dataset pull-download ####
# Connection to XOP
# pull datasets from Profiler

api = ProfilerAPI2::profiler_api(profile = "myprofile")   ## needed for Rmarkdown connection to ProfilerAPI2
con <- api$conn

Views_single_cell <- api$studies$list() %>% 
  filter(grepl("Publication - Single Cell", name)) %>% 
  pull("id")

profiler_folder_id <- api$data_lake$object_id(Views_single_cell, 'TISCH2')

# Datasets to take
# Curated meta data
curated_metalist <- readRDS("./data/TME/curated_metalist.RDS")

dataset_c= melt(curated_metalist) %>% select(dataset) %>% unique() %>% .$dataset
names(curated_metalist)=dataset_c

# Search for the datasets in the TISCH2 folder
Datasets <-foreach(i=dataset_c,.combine=rbind)%do% {
print(i)
  api$data_lake$list(profiler_folder_id, recursive = TRUE) %>%
  dplyr::filter(grepl(".qs", name)) %>%
  dplyr::filter(grepl(i, name)) %>%
  as.data.frame()
}

# Download the datasets
for(i in unique(Datasets$name)){
  
  api$data_lake$download_file(Datasets$id[grep(i,Datasets$name)], path = data_dir, overwrite = TRUE,
                              show_progress = T)
}


# Load the discovery datasets and add meta data
unq_files=list.files(data_dir, pattern = "processed.qs",full.names = T)


data_list=foreach(d=seq_along(unq_files))%do%{
  
  # read seurat files
  seurat <- qread(unq_files[d])
  
  #read metadata
  meta_tisch=seurat@meta.data
  
  #read metadata from original papers
  
  match_dataset=gsub(paste0(data_dir,"/Seurat_"),"",unq_files[d]) %>% gsub("_processed.qs","",.)
  meta_df=curated_metalist[[names(curated_metalist)[grepl(match_dataset,names(curated_metalist))]]]
  meta_df$Cell=rownames(meta_df)
  
  if(unique(meta_df$dataset)=="KIRC_GSE139555"){
    
    intersected_cells=intersect( meta_tisch$Cell,meta_df$Cell)
    meta_tisch_selected=meta_tisch[ meta_tisch$Cell %in% intersected_cells,]
    #
    joined_metadata= left_join(meta_tisch_selected,meta_df,by="Cell")
    rownames(joined_metadata) <- rownames(meta_tisch_selected)
    
    print(paste0(length(intersected_cells)," / ", length(selected_cells)))
    
    seurat_selected=seurat[, rownames(joined_metadata)]
    
    AddMetaData(seurat_selected,
                joined_metadata[,c("dataset","major_celltypes","minor_celltypes")])
  }else{
  
  # Filter only selected columns and extract the rownames  
  selected_cells=sub(".*@", "", rownames(meta_df))
  meta_df=meta_df[,c("Cell","dataset","major_celltypes","minor_celltypes")]
  
  print(unq_files[d])
  print(unique(meta_df$dataset))
  
  # In some rownames @ was added to the rownames in the tisch 
  meta_tisch$Cell=sub(".*@", "", meta_tisch[,"Cell"])
  
  if(unique(meta_df$dataset)=="SKCM_GSE123139"){
    meta_tisch$Cell= sapply(meta_tisch$Cell, function(x) strsplit(x, "_")[[1]][1]) %>%
      unname()
  }
  
  #
  
  intersected_cells=intersect( meta_tisch$Cell,selected_cells)
  meta_tisch_selected=meta_tisch[ meta_tisch$Cell %in% intersected_cells,]
  #
  joined_metadata= left_join(meta_tisch_selected,meta_df,by="Cell")
  rownames(joined_metadata) <- rownames(meta_tisch_selected)

  print(paste0(length(intersected_cells)," / ", length(selected_cells)))
  
seurat_selected=seurat[, rownames(joined_metadata)]

AddMetaData(seurat_selected,
                            joined_metadata[,c("dataset","major_celltypes","minor_celltypes")])
}
}
names(data_list)=gsub(paste0(data_dir,"/Seurat_"),"",unq_files) %>% gsub("_processed.qs","",.)
  
qs::qsave(data_list,paste0(data_dir,"discovery_datalist.qs"))

# Calculate and print out average number of genes per cell
lapply(data_list, function(x){
  cat(paste0("Dataset:" , unique(x@meta.data$dataset)), "\n")
  cat(paste0("No. patients: " , length(unique(x@meta.data$Patient))), "\n")
  cat(paste0("No. cells: " , ncol(x)), "\n")
  cat(paste0("Average no of genes per cell: " ,round(mean(x@meta.data$nFeature_RNA)),1), "\n")
  # cat("mtPercentage summary: ","\n")
  # cat(summary(x@meta.data$percent_mt), "\n")
  # cat("nFeature summary: ","\n")
  # cat(summary(x@meta.data$nFeature_RNA), "\n")
  
}
)

# Dataset:BRCA_GSE176078 
# No. patients: 26 
# No. cells: 43140 
# Average no of genes per cell: 12661 
# Dataset:CRC_GSE166555 
# No. patients: 12 
# No. cells: 13369 
# Average no of genes per cell: 12681 
# Dataset:KIRC_GSE139555 
# No. patients: 3 
# No. cells: 18120 
# Average no of genes per cell: 11861 
# Dataset:LIHC_GSE140228_10X 
# No. patients: 5 
# No. cells: 16724 
# Average no of genes per cell: 13751 
# Dataset:LIHC_GSE140228_Smartseq2 
# No. patients: 6 
# No. cells: 2351 
# Average no of genes per cell: 39531 
# Dataset:NSCLC_GSE131907 
# No. patients: 11 
# No. cells: 25915 
# Average no of genes per cell: 14051 
# Dataset:SKCM_GSE123139 
# No. patients: 8 
# No. cells: 4817 
# Average no of genes per cell: 8861 




