
library(Seurat)
library(Matrix)
library(dplyr)

set.seed(42)

# Script for pre-processing hepato scRNAseq dataset
# Downloaded from GEO with GEO Accession ID: GSE140228
# Smartseq2 downloaded

dirname= "N:./datasets/hepato/"

# Load raw read counts
drop_count= read.table(file=paste0(dirname,"GSE140228_read_counts_Smartseq2.csv"),sep=",",header = T,
                       stringsAsFactors = F)
# convert to sparse matrix
drop_count=Matrix(as.matrix(drop_count[,-1]),sparse = T)
# save as sparse matrix for memory and efficiency
writeMM(Matrix(as.matrix(drop_count[,-1]),sparse = T),
        file = paste0(dirname,"GSE140228_read_counts_Smartseq2.mtx"))

# load gene names
drop_features= read.table(file=paste0(dirname,"GSE140228_gene_info_Smartseq2.tsv"),sep="\t",header = T,
                          stringsAsFactors = F)$SYMBOL
# load metadata
barcodes= read.table(file=paste0(dirname,"GSE140228_cell_info_Smartseq2.tsv"),sep="\t",header = T,
                     stringsAsFactors = F,check.names = F)
rownames(drop_count)=drop_features
colnames(drop_count)=barcodes$Barcode

#getting cell types 
df= barcodes
unique(df$celltype_sub)

#choose only tumor! 
df_f= df %>%
  subset(.$Tissue %in% c("Tumor"))

#cell type annotations
# modify cell type annotations according to the cell types 
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "Lymphoid-B", 
                             "B.cell")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "Lymphoid-B-Plasma", 
                             "plasma")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD4-C3-ANXA1", 
                             "T.CD4.memory")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD4-C4-IL7R", 
                             "T.CD4.memory")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD4-C5-TCF7", 
                             "T.CD4.naive")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD4-C6-CXCL13", 
                             "T.CD4.exhausted")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD4-C7-FOXP3", 
                             "T.CD4.reg")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD8-C2-MKI67", 
                             "T.CD8.proliferative")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD8-C4-CX3CR1", 
                             "T.CD8.effector")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD8-C5-SELL", 
                             "T.CD8.memory")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD8-C6-GZMK", 
                             "T.CD8.effector.memory")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD8-C7-KLRD1", 
                             "T.CD8.tissue.resident.memory")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD8-C8-PDCD1", 
                             "T.CD8.exhausted")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "CD8-C9-SLC4A10", 
                             "T.CD8.MAIT")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "Mono-C1-CD14", 
                             "monocyte")
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "Mono-C2-FCGR3A", 
                             "monocyte.non.classic")
#getting the cell types such as CD8 not other subgrouoings such as CD8-C6-GZMk...
df_f$celltype_sub=unlist(lapply(strsplit(x = df_f$celltype_sub,split = "-"),'[[',1))

# characters for MÏ??? change according to encoding, be careful !
df_f$celltype_sub <- replace(df_f$celltype_sub, 
                             df_f$celltype_sub == "MÏ???", 
                             "Macrophage")
rownames(df_f)= df_f$Barcode

# take samples with cell type info
drop_count=drop_count[,df_f$Barcode]

pbmc <- CreateSeuratObject(counts = drop_count, project = "hepatoma", min.cells = 3,
                           min.features = 200,meta.data = df_f)
dim(pbmc)
#[1] 36263  2423

proj= rep("hepatoma",length(levels(pbmc)))
names(proj) <- levels(pbmc)
pbmc= RenameIdents(pbmc,proj)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter cells and genes based on number of unique features (genes) and mitochondrial gene percentage 
pdf(paste0(dirname,"QC_seurat_filter_plots.pdf"))
VlnPlot(pbmc,features = c("nFeature_RNA"), ncol = 1)
VlnPlot(pbmc,features = c("nFeature_RNA"), ncol = 1,y.max = 10000)

VlnPlot(pbmc,features = c("percent.mt"), ncol = 1)
VlnPlot(pbmc,features = c("percent.mt"), ncol = 1,y.max = 10)

dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 8)
dim(pbmc)
#[1] 36263  2086

saveRDS(pbmc,paste0(dirname,"hepato_counts_filtered_sobj_2.RDS"))

