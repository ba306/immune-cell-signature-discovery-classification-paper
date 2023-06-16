# Preparing benchmarking dataset: Kotliarov PBMC

library(Seurat)
library(dplyr)
library(qs)

dataset_dir="./data/Kotliarov/"

#################################################################
########### Load  scRNAseq data and our gene sets ############

# Load input scRNAseq data 
##from the paper https://figshare.com/articles/dataset/CITE-seq_protein-mRNA_single_cell_data_from_high_and_low_vaccine_responders_to_reproduce_Figs_4-6_and_associated_Extended_Data_Figs_/11349761

pbmc=readRDS( paste0(dataset_dir,"H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"))

# Expression values
exp_ma= pbmc@raw.data
dim(exp_ma)
#Meta data
meta= pbmc@meta.data
dim(meta)
#same order
exp_ma=exp_ma[,rownames(meta)]

# create sobj
tiss_immune=CreateSeuratObject(counts = exp_ma,
                               meta.data = meta)
tiss_immune[["percent.mt"]] <- PercentageFeatureSet(tiss_immune, pattern = "^MT-")

VlnPlot(tiss_immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

tiss_immune <- subset(tiss_immune, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# Cell type harmonization
# Annotate clusters from the paper

tiss_immune$K1=gsub("C0","T CD4 naive/DNT",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C1","T CD4 memory",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C2","Monocyte CD14 / DC",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C3","B",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C4","T CD8 memory",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C5","NK",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C6","T CD8 naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C7","T unconventional",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C8","Monocyte CD16",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("C9","pDC",tiss_immune$K1,fixed = T)

tiss_immune$K3=gsub("C0.0.0","T CD4 naive",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C0.1.0","T dn",tiss_immune$K3,fixed = T)

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
tiss_immune$K3=gsub("C7.0.0","T unconv CD161+/CD3+ dn T",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C7.0.1","T unconv CD161hi/CD3+/CD8+",tiss_immune$K3,fixed = T)

tiss_immune$K3=gsub("C8.0.0","Monocyte CD16",tiss_immune$K3,fixed = T)
tiss_immune$K3=gsub("C9.0.0","pDC",tiss_immune$K3,fixed = T)

# Harmonnized cell types
original=tiss_immune$K1
table(original)
original[tiss_immune$K3=="DC"]="DC"
original[tiss_immune$K3=="T dn"]="T unconv"
original[tiss_immune$K3=="T unconv CD161+/CD3+ dn T"]="T unconv"
original[tiss_immune$K3=="T unconv CD161hi/CD3+/CD8+"]="T unconv"
original[tiss_immune$K3=="T CD8 NKT-like"]="T unconv"

original=gsub("T unconventional","T unconv",original,fixed = T)

original[tiss_immune$K3=="HSC"]="HSC"

original=gsub("pDC","DC",original,fixed = T)
original=gsub("Monocyte CD14 / DC","Mono",original,fixed = T)
original=gsub("Monocyte CD16","Mono",original,fixed = T)
original=gsub("T CD4 memory","T CD4",original,fixed = T)
original=gsub("T CD4 naive/DNT","T CD4",original,fixed = T)
original=gsub("T CD8 memory","T CD8",original,fixed = T)
original=gsub("T CD8 naive","T CD8",original,fixed = T)

# check if grouped well

tiss_immune$harmonized_celltype=original
table(tiss_immune$K3,original)

# 2 cells seem to be annotated wrong... remove those cells
colnames(tiss_immune)[original=="T unconv" & tiss_immune$K1=="Monocyte CD14 / DC"]
tiss_immune@meta.data["CACATTTGTGGAAAGA_H1B2ln5",]
tiss_immune=tiss_immune[,-grep("CACATTTGTGGAAAGA_H1B2ln5",colnames(tiss_immune))]
#Monocyte CD14 / DC C2.0 T unconv CD161+/CD3+ dn T 

colnames(tiss_immune)[original=="HSC" & tiss_immune$K1=="T CD4 naive/DNT"]
tiss_immune@meta.data["GTATTCTAGAAGGACA_H1B2ln1",]
tiss_immune=tiss_immune[,-grep("GTATTCTAGAAGGACA_H1B2ln1",colnames(tiss_immune))]

tiss_immune=NormalizeData(tiss_immune)
dim(tiss_immune)
#32738 genes and 52849 cells

# average number of genes per cell
mean(tiss_immune@meta.data$nFeature_RNA)
#748.0894

#save
qsave(tiss_immune,paste0(dataset_dir,"kotliarov_pbmc.qs"))

