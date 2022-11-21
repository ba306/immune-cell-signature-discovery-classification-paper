
library(fst)
library(dplyr)
library(varhandle)
library(foreach)
library(Seurat)

dirname= "N:/"

# Script for pre-processing and saving 2 melanoma scRNAseq datasets as Seurat objects
# GEO accessions: GSE115978 and GSE123139
# Datasets previously processed and saved by Dr. Sheng Zhao
  # computesum() applied previously

# load fst files
melanoma_2= read.fst("Melanoma_2018_GSE115978_20190514.fst")

li= read.fst("Li_et_al_Cell_2018_Melanoma_20190312.fst")

# get metadata especially cell types 
metadata1 <- dplyr::select(li, cell, cell_type)
metadata2 <- dplyr::select(melanoma_2, cells, samples, cell_type, treatment.group, Cohort, no.of.genes, no.of.reads)


metadata1$cell_type <- replace(metadata1$cell_type, metadata1$cell_type == "B", "B.cell")
metadata1$cell_type <- replace(metadata1$cell_type, metadata1$cell_type == "macrophage", "Macrophage")
metadata2$cell_type <- replace(metadata2$cell_type, metadata2$cell_type == "?", "NA")

rownames(metadata1)=metadata1$cell
rownames(metadata2)=metadata2$cells


li <- t(as.matrix(dplyr::select(li, -cell, -cell_type)))
colnames(li)=metadata1$cell

melanoma_2 <- t(as.matrix(dplyr::select(melanoma_2, -(1:7))))
colnames(melanoma_2)=metadata2$cells


# Create seurat objects

li_seurat <- CreateSeuratObject(li, assay = "RNA", meta.data = metadata1, 
                                min.cells = 3, min.features = 200)
#save seurat

li_seurat[["percent.mt"]] <- PercentageFeatureSet(li_seurat, pattern = "^MT-")

# Filter cells and genes based on number of unique features (genes) and mitochondrial gene percentage 
pdf(paste0(dirname,"QC_seurat_filter_plots_li.pdf"))
VlnPlot(li_seurat,features = c("nFeature_RNA"), ncol = 2)
VlnPlot(li_seurat,features = c("nFeature_RNA"), ncol = 2,y.max = 4000)
VlnPlot(li_seurat,features = c("percent.mt"), ncol = 2,y.max = 10)
dev.off()

li_seurat <- subset(li_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 )
saveRDS(li_seurat,"N:/Li_et_al_Cell_2018_Melanoma_20190312.RDS")


# Create seurat objects

melanoma_2_seurat <- CreateSeuratObject(melanoma_2, assay = "RNA", meta.data = metadata2, 
                                        min.cells = 3, min.features = 200)

# Filter cells and genes based on number of unique features (genes) and mitochondrial gene percentage 

melanoma_2_seurat[["percent.mt"]] <- PercentageFeatureSet(melanoma_2_seurat, pattern = "^MT-")

pdf(paste0(dirname,"QC_seurat_filter_plots_melanoma_2.pdf"))
VlnPlot(melanoma_2_seurat,features = c("nFeature_RNA"), ncol = 2,y.max = 12000)
VlnPlot(melanoma_2_seurat,features = c("percent.mt"), ncol = 2,y.max = 10)
dev.off()

melanoma_2_seurat <- subset(melanoma_2_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 )
saveRDS(melanoma_2_seurat,"N:/Melanoma_2018_GSE115978_20190514.RDS")



