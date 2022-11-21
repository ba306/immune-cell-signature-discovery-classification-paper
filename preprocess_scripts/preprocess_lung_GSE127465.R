
library(Matrix)
library(foreach)
library(Seurat)
library(dplyr)


## Script for pre-processing lung scRNAseq dataset
#directory 
dirname= "N:./datasets/myeloid/"


# Downloaded from GEO with GEO Accession ID:  GSE127465
# only human samples were used!

# Combine all expression files 

files= list.files(paste0(dirname,'GSE127465_RAW/'),pattern = ".tsv",full.names = T)

# Save tsv files into sparse matrix (e.g. GSM3635278_human_p1t1_raw_counts.mtx)
# and also save barcode names for each patient and sample (e.g. "p1t1_bcFSNQ") into a vector barcodes_names 
barcodes_names=foreach(i=seq_along(files),.combine = c)%do%{
print(i)
count_matrix=Matrix(as.matrix(read.table(file=files[i],sep="\t",header = F,
                              stringsAsFactors = F,row.names = 1,skip = 1)),sparse=T)
file_name=gsub("N:./myeloid/GSE127465_RAW/(.*).tsv","\\1",files[i])
writeMM(obj = count_matrix,
        file = paste0(dirname,"matrix/",file_name,".mtx"))
patient_sample= sapply(strsplit(files[i],split = "human_*_"),'[[',2)
patient_sample=sapply(strsplit(patient_sample,split = "_raw"),'[[',1)
print(paste0(patient_sample))
paste0(patient_sample,"_",row.names(count_matrix))
}

#save patient_sample barcodes 
df= data.frame(barcodes_names)
write.table(df,file=paste0(dirname,"final_test/barcodes.tsv"),sep = "\t",row.names = F)

# Read all expression files and combine them into sparse matrix format matrix.mtx
files_m= list.files(paste0(dirname,'matrix/'),pattern = ".mtx",full.names = T)

mm= foreach(i=seq_along(files_m),.combine = rbind) %do% {
  print(i)
  a=readMM(files_m[i])
}

# Genes x cell
writeMM(t(mm),file = paste0(dirname,"/final_test/matrix.mtx"))

# Preprocess dataset and save as seurat object
#read counts
# one can proceed with using t(mm) as counts if memory is not a problem
counts= readMM(paste0(dirname,"/final_test/matrix.mtx"))
dim(counts)
#[1] 41861 173954

# Gene names GSE127465_gene_names_human_41861.tsv have been saved as features.tsv 
features= read.table(file=paste0(dirname,"/final_test/features.tsv"),sep="\t",header = F,
                     stringsAsFactors = F)
# Read barcodes 
#barcodes= read.table(file=paste0(dirname,"/final_test/barcodes.tsv"),sep="\t",header = T,
 #                    stringsAsFactors = F,check.names = F)
rownames(counts)=features$V1
colnames(counts)=barcodes$barcodes_names

#Meta data human_cell_metadata_54773x25.tsv have been saved as annot.tsv
df= read.table(file=paste0(dirname,"final_test/annot.tsv"),sep="\t",header = T,
               stringsAsFactors = F)
lib_bar= paste0(df$Library,"_",df$Barcode)

# take samples with annotation 
counts_final= counts[,lib_bar]
dim(counts_final)
#[1] 41861 54773

#cell type annotations
meta.data= data.frame(cell.types=df$Most.likely.LM22.cell.type)
rownames(meta.data)=lib_bar

#seurat object 
pbmc <- CreateSeuratObject(counts = counts_final, project = "myeloid_lung", min.cells = 3,
                           min.features = 200,meta.data = meta.data)
dim(pbmc)
#[1] 36555 53215

proj= rep("myeloid_lung",length(levels(pbmc)))
names(proj) <- levels(pbmc)
pbmc= RenameIdents(pbmc,proj)

# Filter cells and genes based on number of unique features (genes) and mitochondrial gene percentage 

#mt content percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pdf(paste0(dirname,"QC_seurat_filter_plots.pdf"))
#VlnPlot(pbmc,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(pbmc,features = c("nFeature_RNA"), ncol = 1,y.max = 5000)
VlnPlot(pbmc,features = c("percent.mt"), ncol = 1,y.max = 20)

dev.off()

pdf(paste0(dirname,"QC_seurat_filter_plots_2.pdf"))
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(pbmc)
#[1] 36555 40875

#take only tumor samples
df_tumor= df%>%
  subset(Tissue=='tumor')
tumor= paste0(df_tumor$Library,"_",df_tumor$Barcode)

inter= intersect(colnames(pbmc),tumor)
pbmc_tumor=pbmc[,inter]
dim(pbmc_tumor)
#[1] 36555 27473

saveRDS(pbmc_tumor,paste0(dirname,"myeloid_raw_tumor_filtered_sobj.RDS"))


