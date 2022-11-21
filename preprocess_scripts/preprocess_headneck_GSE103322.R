
library(dplyr)
library(Seurat)
library(foreach)
library(Matrix)

# Script for pre-processing head and neck scRNAseq dataset
# Downloaded from with GEO Accession ID: GSE103322


#directory 
dirname <- "N:./datasets/puram/"

#load data
data= read.table(paste0(dirname,"GSE103322_HNSCC_all_data.txt"),
                 stringsAsFactors = F,sep = "\t",
                 header = T,row.names = 1)

dim(data)
#[1] 23691  5902

row.names(data)[1:5]
#[1] "processed by Maxima enzyme"     "Lymph node"                    
#[3] "classified  as cancer cell"     "classified as non-cancer cells"
#[5] "non-cancer cell type"     

# first five rows are the annotations
# save annotations/metadata
df= data[1:5,]
t_df=as.data.frame(t(df))

# data without metadata saved 
data=data[-c(1:5),]


# save data as sparse matrix
data_s= Matrix(as.matrix(data),sparse=T)

#take tumor not lymph node !
t_df=t_df[t_df$`Lymph node`==0,]
#and also take immune cells 
t_df=t_df[t_df$`non-cancer cell type`!=0,]
data_s=data_s[,rownames(t_df)]


#seurat object 
pbmc <- CreateSeuratObject(counts = data_s, project = "HNSCC", 
                           min.cells = 3, min.features = 200,meta.data = t_df[,"non-cancer cell type",drop=F])
dim(pbmc)
#[1] 20011  2817

proj= rep("HNSCC",length(levels(pbmc)))
names(proj) <- levels(pbmc)
pbmc= RenameIdents(pbmc,proj)

# Filter cells and genes based on number of unique features (genes) and mitochondrial gene percentage 

#mt content percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pdf(paste0(dirname,"QC_seurat_filter_plots_puram.pdf"))
VlnPlot(pbmc,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(pbmc,features = c("nFeature_RNA"), ncol = 2,y.max = 10000)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
dim(pbmc)
#[1] 20011  2790

#save seurat object 
saveRDS(pbmc,paste0(dirname,"HNSCC_2_filtered_sobj.RDS"))
