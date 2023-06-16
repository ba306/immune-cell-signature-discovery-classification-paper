# Integrate each sample in IFN stimulation PBMC dataset from Kartha
# for further input into Seurat and scType 

library(dplyr)
library(ggplot2)
library(qs)
library(Seurat)
library(tibble)

save_dir="./data/Kartha/"
pbmc_cell_type=qread(paste0(save_dir,"Kartha_sobj.qs"))


# Take only IFNg and control
take=pbmc_cell_type$StimType=="Control" | pbmc_cell_type$StimType=="IFN"
pbmc_cell_type=pbmc_cell_type[,take]

table(pbmc_cell_type$Condition,pbmc_cell_type$Donor)

# Remove ControlGolgiPlug_6h and IFNGolgiPlug_6h
# remove donor 4 since this only has IFNg treatment

remove=pbmc_cell_type$Condition=="ControlGolgiPlug_6h" | pbmc_cell_type$Condition=="IFNGolgiPlug_6h" | pbmc_cell_type$Donor =="Donor4"
pbmc_cell_type=pbmc_cell_type[,!remove]

# Follow Seurat RPCA integration workflow
# Separate by the conditions/samples
pbmc_list <- SplitObject(pbmc_cell_type, split.by = "Condition")

for (i in 1:length(pbmc_list)) {
  print(i)
  pbmc_list[[i]] <- NormalizeData(pbmc_list[[i]], verbose = F)
  pbmc_list[[i]] <- FindVariableFeatures(pbmc_list[[i]], selection.method = "vst", nfeatures = 2000,
                                             verbose = F)
}

features <- SelectIntegrationFeatures(object.list = pbmc_list)

pbmc_list <- lapply(X = pbmc_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = pbmc_list, reduction = "rpca", 
                                  dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(pancreas.integrated)="integrated"


qsave(pancreas.integrated,file =paste0(save_dir, "Kartha_integrated.qs"))
