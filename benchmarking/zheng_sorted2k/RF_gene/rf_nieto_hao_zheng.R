# RF Nieto genes
# Ref Hao
# Benchmarking Zheng

library(Seurat)
library(dplyr)
library(readxl)
library(foreach)
library(randomForest)
source("N:./immunsig_scripts_manuscript/benchmarking/celltype_annot_stats.R")

set.seed(42)

dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/zheng_sorted2k/RF_gene/"


Nieto= read_excel(path = paste0("N:./immunsig_scripts_manuscript/gene_lists_papers/",
                                "/nieto_immune_2021.xlsx"),sheet = 1)
selected=c("B cells"       ,  "Plasma B cells"   ,        "T cells naive"   ,
           "CD4 effector memory"    ,"CD4 transitional memory" ,
           "CD4 naïve-memory"    ,       "CD8 effector memory"  ,"NK"     ,
           "Monocytes"      ,          "cDC"        ,              "pDC"      ,                "mDC" )
genes=Nieto[Nieto$`Cell type` %in% selected,] %>%.$Markers %>%
  strsplit(", ") %>%unlist() %>% unique()
length(genes) #85

# Load input scRNAseq data
tiss_immune <- readRDS("N:/datasets/PBMC/zheng_sorted/pbmc_zheng_sorted_sobj2k.RDS")

genes_take_cell_types=intersect(rownames(tiss_immune),genes)

paste(length(genes_take_cell_types),"/", length(genes))
# [1] "79 / 85"

tiss_immune=tiss_immune[genes_take_cell_types,]
#
# # Load reference data
reference_seurat <- readRDS("Z:/pbmc_hao_ref.RDS")

genes_take_cell_types=intersect(rownames(reference_seurat),genes)

paste(length(genes_take_cell_types),"/", length(genes))
#[1] "83 / 85"
reference_seurat=reference_seurat[genes_take_cell_types,]

common_genes=intersect(
  rownames(tiss_immune),
  rownames(reference_seurat)
)

length(common_genes) #79


mean_sig_func=function(data,assay,common_genes){
  data=data[common_genes,]
  data <- ScaleData(data, features = rownames(data))
  t(data@assays[[assay]]@scale.data)

}
#
hao_mean=mean_sig_func(reference_seurat,"RNA",common_genes)

hao_mean=cbind(hao_mean,data.frame(cell_type=reference_seurat$harmonized_celltype))
saveRDS(hao_mean,paste0(dir_analysis,"hao_nietogenes.RDS"))
hao_mean=readRDS(paste0(dir_analysis,"hao_nietogenes.RDS"))

### RF signatures ##
mean_all=hao_mean
table(mean_all$cell_type)
# Prepare data for random forest
input= mean_all
input=data.frame(input)

input$cell_type=factor(input$cell_type)

#down_sample=sample(1:nrow(input), 50000, replace=TRUE)
#table(input[down_sample,"cell_type"])
#input=input[down_sample,]
# Random Forest 
# Seperate cells into train and test set 67:33
ss <- sample(1:2,size=nrow(input),replace=TRUE,prob=c(0.67,0.33))

train= input[ss==1,]
test= input[ss==2,]

#train.class=input[ss==1,which(colnames(input) %in% "cell")]
#input[,"cell"]

#train=cbind(train,train.class)

# Train the data 
start.time <- Sys.time()

rf <- randomForest(
  cell_type ~ .,
  data=train,importance = TRUE
)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken  %>% print()

saveRDS(rf,paste0(dir_analysis,"rf_nietogenes.RDS"))

rf=readRDS(paste0(dir_analysis,"rf_nietogenes.RDS"))
table(rf$predicted)

train=input[names(rf$predicted),]
test=input[setdiff(rownames(input),names(rf$predicted)),]

table(train$cell_type)


df_filter= data.frame(original=rf$y,predicted=rf$predicted)
df_filter$original=unfactor(df_filter$original)
df_filter$predicted=unfactor(df_filter$predicted)
Cell_type_stats(df_filter)

# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1    93.99 94.03 98.65       98.27    97.59    94.01

# Predict cell types using trained RF model on test set
pred = predict(rf, test)


df_filter= data.frame(original=test$cell_type,predicted=pred)
df_filter$original=unfactor(df_filter$original)
df_filter$predicted=unfactor(df_filter$predicted)
Cell_type_stats(df_filter)
# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1      94.02 94.03 98.64       98.27    97.58    94.02

### pred 
common_genes=colnames(hao_mean)[-ncol(hao_mean)]
table(tiss_immune$cell_Types)

tiss_immune$K1=tiss_immune$cell_Types
tiss_immune$K1=gsub("Monocytes_CD14","Mono",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD4_memory","T CD4 memory_naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD4_naive","T CD4 memory_naive",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD8_cytotoxic","T CD8",tiss_immune$K1,fixed = T)
tiss_immune$K1=gsub("TCD8_naive","T CD8",tiss_immune$K1,fixed = T)
table(tiss_immune$K1)

ko_mean=mean_sig_func(tiss_immune,"RNA",common_genes)
ko_mean=cbind(ko_mean,data.frame(cell_type=tiss_immune$K1))
colnames(ko_mean)=gsub("-",".",colnames(ko_mean),fixed = T)

table(ko_mean$cell_type)

start.time <- Sys.time()

pred = predict(rf, ko_mean[,-ncol(ko_mean)])
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# saveRDS(pred,paste0(dir_anaylsise,"pred_ko_mean.RDS"))

#accuracy harmonize
df_filter= data.frame(original= ko_mean$cell_type,predicted=pred)
df_filter$original=unfactor(df_filter$original)
df_filter$predicted=unfactor(df_filter$predicted)

table(df_filter$original)
table(df_filter$predicted)
table(pred)
table(df_filter$original,df_filter$predicted)

Cell_type_stats(df_filter)
# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1      79.23 81.28 92.71       92.18    88.73    80.24
