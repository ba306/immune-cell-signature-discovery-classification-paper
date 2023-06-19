# Integration script for discovery datasets
# Seurat RPCA integration based on 3k HVGs 

library(Seurat)
library(qs)

set.seed(42)

# Data directory for discovery datasets
data_dir="./data/TME/"

# Analysis
dir="./discovery/"

# Functions for discovery workflow
source("./functions/discovery_functions.R")

# define parameters
no_HVGs=3000


############### Load and prepare datasets ########################

# Load the datasets 

data_list <- qs::qread(paste0(data_dir,"discovery_datalist.qs"))

#check the cell type compositions from each dataset
# Violin plot for cell type compositions from different datasets

#remove NKT cells --> only annotated at BRCA dataset
data_list$BRCA_GSE176078=data_list$BRCA_GSE176078[,data_list$BRCA_GSE176078$major_celltypes!="NKT"]

# Stacked bar plot for the contribution of each dataset to each major cell type
png(paste0(dir,"Vln_celltyecomposition.png"),width = 800,height = 400)
violin_celltype_dataperc(data_list,by = "dataset")
dev.off()

######################## Integration #############################

# find variable features 
data_list <- lapply(data_list, function(x) {
  FindVariableFeatures(x, selection.method = "vst", 
                       nfeatures = no_HVGs, verbose = T)
  
})

#find common intgeration features (genes)
features <- SelectIntegrationFeatures(object.list = data_list,nfeatures = no_HVGs)

# scale datasets individually
# run pca on  selected HVGs individually

data_list <- lapply(data_list, function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T,npcs = 20)
})

# reciprocal pca ('rpca') method was used  for integration
anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:20, 
                                  anchor.features = features,
                                  reduction = "rpca")

# Integrate datasets using the anchors
integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(integrated) <- "integrated"

# Scale integrated dataset
integrated <- ScaleData(integrated, verbose = T)

# # Check integration on UMAP plot
integrated <- RunPCA(integrated, npcs = 20, verbose = T)
ElbowPlot(integrated)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:10)


plot_dataset=DimPlot(integrated,pt.size = 0.1,group.by = "dataset",reduction = "umap",
                     label=F,repel = T,cols=chr )+
  guides(shape = guide_legend(override.aes = list(size = 2),
                              nrow=1),
         fill = guide_legend(override.aes = list(size = 2),
                             nrow=1))+ 
  theme(legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),legend.position = "top",
        plot.title = element_blank())

plot_majorcelltypes=DimPlot(integrated,pt.size = 0.1,group.by = "major_celltypes",reduction = "umap",
                            label=T,repel = T,cols=chr,label.box = T,label.color = "white",label.size = 3) +
  guides(shape = guide_legend(override.aes = list(size = 2),
                              nrow=1),
         fill = guide_legend(override.aes = list(size = 2),
                             nrow=1))+ 
  theme(legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),legend.position = "top",
        plot.title = element_blank())

png(paste0(dir,"umap_integration_plots_all.png"),
    width = 25, height = 12, units = "cm", res = 300)

grid.arrange(plot_majorcelltypes, plot_dataset,nrow=1)
dev.off()

#save the integrated sobj for downstream analysis
qsave(integrated, paste0(dir,"integrated_discovery.qs"))

