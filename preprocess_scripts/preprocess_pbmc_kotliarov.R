library(Seurat)
library(dplyr)
options(ggrepel.max.overlaps = Inf)
set.seed(42)
dataset_dir="N:/datasets/PBMC/kolitrov/"


dir_analysis="N:/PBMC_class/"
dir.create(paste0(dir_analysis,"kolitrov/"))
dirname=paste0(dir_analysis,"kolitrov/")

#################################################################
########### Load  scRNAseq data and our gene sets ############

# Load input scRNAseq data 
##from the paper https://figshare.com/articles/dataset/CITE-seq_protein-mRNA_single_cell_data_from_high_and_low_vaccine_responders_to_reproduce_Figs_4-6_and_associated_Extended_Data_Figs_/11349761

pbmc=readRDS( paste0(dataset_dir,"H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"))
exp_ma= pbmc@raw.data
dim(exp_ma)
meta= pbmc@meta.data
dim(meta)
exp_ma=exp_ma[,rownames(meta)]

# pbmc_so=CreateSeuratObject(counts = exp_ma,meta.data = meta)
# pbmc_so@assays$RNA@data=pbmc@data[,rownames(meta)]
# saveRDS(pbmc_so,paste0(dataset_dir,"kolit_pbmc_paper_norm.RDS"))

tiss_immune=CreateSeuratObject(counts = exp_ma,
                               meta.data = meta)
tiss_immune[["percent.mt"]] <- PercentageFeatureSet(tiss_immune, pattern = "^MT-")

# VlnPlot(tiss_immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(tiss_immune$nCount_RNA)
summary(tiss_immune$nFeature_RNA)#min 120, max 2342 
summary(tiss_immune$percent.mt)

tiss_immune <- subset(tiss_immune, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)


# tiss_immune@assays$RNA@data=pbmc@data[,rownames(meta)]
tiss_immune=NormalizeData(tiss_immune)

tiss_immune=tiss_immune[,-grep("C7",tiss_immune$K1)]
tiss_immune=tiss_immune[,-grep("C7.0.0",tiss_immune$K3)]


tiss_immune$K1=gsub("C0","T CD4 naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C1","T CD4 memory",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C2","Monocyte CD14 / DC",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C3","B",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C4","T CD8 memory",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C5","NK",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C6","T CD8 naive",tiss_immune$K1,fixed = T)
# tiss_immune$K1=gsub("C7","T unconventional",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C8","Monocyte CD16",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C9","pDC",tiss_immune$K1,fixed = T)

tiss_immune$K3=gsub("C0.0.0","T CD4 naive",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C0.1.0","T CD4 naive double negative T",tiss_immune$K3,fixed = T)

tiss_immune$K3=gsub("C1.0.0","T CD4 central/transitional memory",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C1.1.0","T CD4 TEMRA/effector memory",tiss_immune$K3,fixed = T)

tiss_immune$K3=gsub("C2.0.0","Monocyte CD14",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C2.0.1","Monocyte CD14 IgA",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C2.0.2","HSC",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C2.1.0","DC",tiss_immune$K3,fixed = T)

tiss_immune$K3=gsub("C3.0.0","B transitional",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C3.0.1","B unswitched",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C3.1.0","B switched",tiss_immune$K3,fixed = T)


tiss_immune$K3=gsub("C4.0.0","T CD8 central/transitional memory",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C4.0.1","T CD8 TEMRA/effector memory",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C4.0.2","T CD8 NKT-like",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C4.0.3","T CD8 CD103+",tiss_immune$K3,fixed = T)


tiss_immune$K3=gsub("C5.0.0","NK CD16hi",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C5.0.1","NK CD16lo/CD56hi",tiss_immune$K3,fixed = T)

tiss_immune$K3=gsub("C6.0.0","T CD8 naive",tiss_immune$K3,fixed = T)
# tiss_immune$K3=gsub("C7.0.0","T unconv CD161+/CD3+",tiss_immune$K3,fixed = T)
# tiss_immune$K3=gsub("C7.0.1","T unconv CD161hi/CD3+/CD8+",tiss_immune$K3,fixed = T)

tiss_immune$K3=gsub("C8.0.0","Monocyte CD16",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C9.0.0","pDC",tiss_immune$K3,fixed = T)

saveRDS(tiss_immune,paste0(dataset_dir,"kotliarov_pbmc_lognorm_sobj2.RDS"))





