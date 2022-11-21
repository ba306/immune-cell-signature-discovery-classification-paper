library(Seurat)
library(ggplot2)
options(ggrepel.max.overlaps = Inf)
# pbmc68k preprocess


dataset_dir="N:/datasets/PBMC/pbmc68k/"

##################################################################
########### Load  scRNAseq data and our gene sets ############

# Load input scRNAseq data 
#https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/pbmc68k_data.rds

tiss_immune=readRDS(paste0(dataset_dir,"pbmc68k_data.rds"))
mat=tiss_immune$all_data$`17820`$hg19$mat
mat = as(mat, "dgCMatrix") 
mat=t(mat)
dimnames(mat) = list(tiss_immune$all_data$`17820`$hg19$gene_symbols,
                     tiss_immune$all_data$`17820`$hg19$barcodes)
mat[1:2,1:2]

meta=read.table(paste0(dataset_dir,"68k_pbmc_barcodes_annotation.tsv"),sep = "\t",
                header=T,row.names = 3)
meta=meta[colnames(mat),]

pbmc <- CreateSeuratObject(counts = mat, project = "Pbmc-68k", 
                           min.cells = 3,
                           min.features = 200,meta.data = meta)
dim(pbmc)
#[1] 17788 68551
# already prefiltered by the curators
# saveRDS(pbmc,paste0(dataset_dir,"PBMC_68k.rds"))

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
summary(pbmc$nFeature_RNA)#min 200, max 2739 
summary(pbmc$percent.mt) #0-26.4
pbmc <- subset(pbmc, subset = nFeature_RNA < 2000 & percent.mt < 10)

pbmc=NormalizeData(pbmc)
saveRDS(pbmc,paste0(dataset_dir,"PBMC_68k_preprocessed.rds"))
