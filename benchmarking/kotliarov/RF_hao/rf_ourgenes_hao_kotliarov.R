# RF our genes
# Ref Hao
# Benchmarking Kotliarov

library(Seurat)
library(dplyr)
library(readxl)
library(foreach)
library(randomForest)
library(varhandle)
source("N:./immunsig_scripts_manuscript/benchmarking/celltype_annot_stats.R")
source("N:./immunsig_scripts_manuscript/benchmarking/kotliarov/kotliarov_harmonize_prediction_scores_func.R")

set.seed(42)

dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/kotliarov/RF_hao/"

#Load our gene signatures
db_genes= read.table("N:./immunsig_scripts_manuscript/immune_signature_discovery_workflow/dbscan_cluster_genes_param_0_GO_annot.tsv",sep = "\t",
                     stringsAsFactors = F)

db_genes_filter= db_genes %>%
  subset(Go %in%unique(db_genes$Go)[c(4,5,6,7,8,11,16,19)] )
genes= rownames(db_genes_filter)

length(genes)

# Load input scRNAseq data
tiss_immune <- readRDS("N:/datasets/PBMC/kolitrov/kotliarov_pbmc_lognorm_sobj2.RDS")

genes_take_cell_types=intersect(rownames(tiss_immune),genes)

paste(length(genes_take_cell_types),"/", length(unique(genes)))
#[1] "201 / 209"


tiss_immune=tiss_immune[genes_take_cell_types,]
#
# # Load reference data
reference_seurat <- readRDS("Z:/pbmc_hao_ref.RDS")

genes_take_cell_types=intersect(rownames(reference_seurat),genes)

paste(length(genes_take_cell_types),"/", length(unique(genes)))

reference_seurat=reference_seurat[genes_take_cell_types,]
#[1] "200 / 209"
common_genes=intersect(
  rownames(tiss_immune),
  rownames(reference_seurat)
)

length(common_genes)
#199

mean_sig_func=function(data,assay,common_genes){
  data=data[common_genes,]
  data <- ScaleData(data, features = rownames(data))
  t(data@assays[[assay]]@scale.data)

}
#
hao_mean=mean_sig_func(reference_seurat,"RNA",common_genes)

hao_mean=cbind(hao_mean,data.frame(cell_type=reference_seurat$harmonized_celltype))
saveRDS(hao_mean,paste0(dir_analysis,"hao_ourgenes.RDS"))
hao_mean=readRDS(paste0(dir_analysis,"hao_ourgenes.RDS"))

### RF signatures ##
mean_all=hao_mean
table(mean_all$cell_type)
# Prepare data for random forest
input= mean_all
input=data.frame(input)

input$cell_type=factor(input$cell_type)

# Random Forest 
# Seperate cells into train and test set 67:33
ss <- sample(1:2,size=nrow(input),replace=TRUE,prob=c(0.67,0.33))

train= input[ss==1,]
test= input[ss==2,]


# Train the data 
start.time <- Sys.time()

rf <- randomForest(
  cell_type ~ .,
  data=train,importance = TRUE
)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken  %>% print()

saveRDS(rf,paste0(dir_analysis,"rf_ourgenes.RDS"))


rf=readRDS(paste0(dir_analysis,"rf_ourgenes.RDS"))
table(rf$predicted)

train=input[names(rf$predicted),]
test=input[setdiff(rownames(input),names(rf$predicted)),]

table(train$cell_type)


df_filter= data.frame(original=rf$y,predicted=rf$predicted)
df_filter$original=unfactor(df_filter$original)
df_filter$predicted=unfactor(df_filter$predicted)
Cell_type_stats(df_filter)

# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1       98.17 98.17 99.59       99.51     99.3    98.17

# Predict cell types using trained RF model on test set
pred = predict(rf, test)


df_filter= data.frame(original=test$cell_type,predicted=pred)
df_filter$original=unfactor(df_filter$original)
df_filter$predicted=unfactor(df_filter$predicted)
Cell_type_stats(df_filter)


### pred 
common_genes=colnames(hao_mean)[-ncol(hao_mean)]

ko_mean=mean_sig_func(tiss_immune,"RNA",common_genes)
ko_mean=cbind(ko_mean,data.frame(cell_type=tiss_immune$K1))
colnames(ko_mean)=gsub("-",".",colnames(ko_mean),fixed = T)

table(ko_mean$cell_type)

start.time <- Sys.time()

pred = predict(rf, ko_mean[,-ncol(ko_mean)])
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#accuracy harmonize

general_harmonnize_pred_stat(new_clusters=pred,tiss_immune)
# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1          90 90.17 95.77       96.62     94.9    90.08


