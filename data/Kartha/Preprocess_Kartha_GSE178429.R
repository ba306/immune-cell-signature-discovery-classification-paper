# Preprocessing PBMC treatment data with different stimuli
# IFN-g, IFN GolgiPlug, LPS, LPSGolgiPlug, PMA and PMAGolgiPlug
# We will for now preprocess all and then only focus on IFNg and control

library(tidyverse)
library(stringr)
library(data.table)
library(Matrix)
library(qs)
library(Seurat)

file_dir="./data/Kartha/"

# Get the meta data
meta=read.table(paste0(file_dir,"GSE178429_PBMCs_stim_scRNAseq_cellMeta.txt"),header = T)
meta$meta_number=rownames(meta) %>% as.integer()
meta$time=gsub(x = meta$Condition,pattern =   "(.*)_(.*)","\\2")
meta$Stimulation=gsub(x = meta$Condition,pattern =   "(.*)_(.*)","\\1")
rownames(meta)=meta$cellBarcode

# Get the gene names
genes=read.table(paste0(file_dir,"GSE178429_PBMCs_stim_scRNAseq_geneNames.txt"))
genes$gene_number=rownames(genes) %>% as.integer()

# Get the count matric
counts=read.table(paste0(file_dir,"GSE178429_PBMCs_stim_scRNAseq_counts.txt"),header = T,skip = 1)
# first column gene number, second column cell number and third column counts
colnames(counts)=c("gene_number","meta_number","count")

counts_genes_meta=left_join(genes,counts) %>% 
  left_join(meta) %>%
  select(V1,count,cellBarcode)

# Merge Count matrix with annotations
counts_genes_meta_dt=as.data.table(counts_genes_meta)
counts_combined=data.table::dcast(counts_genes_meta_dt, V1  ~ cellBarcode, value.var = "count")%>% 
  column_to_rownames(var="V1")%>% rename_with(~str_remove(., 'count.'))

#NA values are 0 counts
counts_combined <- mutate_all(counts_combined, ~replace(., is.na(.), 0))
# counts_combined_sparse=Matrix(as.matrix(counts_combined), sparse=TRUE)

dim(counts_combined)
#[1] 19222 23754

pbmc.big <- CreateSeuratObject(counts = counts_combined,meta.data =meta ,min.cells = 3, 
                               min.features = 200) 
pbmc.big$project="IFN"

dim(pbmc.big)
#[1] 16949 23754

pbmc.big[["percent.mt"]] <- PercentageFeatureSet(pbmc.big, pattern = "^MT-")

VlnPlot(pbmc.big, features = "nFeature_RNA",pt.size = 0,group.by = "orig.ident",split.by = "Condition")
VlnPlot(pbmc.big, features = "percent.mt",pt.size = 0,group.by = "orig.ident",split.by = "Condition")

# log Normalize
pbmc.big=NormalizeData(pbmc.big)

# save as qs
qsave(pbmc.big,paste0(file_dir,"Kartha_sobj.qs"))

