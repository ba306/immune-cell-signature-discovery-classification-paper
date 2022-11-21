
library(Seurat)
library(dplyr)
set.seed(42)

source("N:./immunsig_scripts_manuscript/benchmarking/kotliarov/kotliarov_harmonize_prediction_scores_func.R")

# Script for cell type annotation/prediction using singleR package

# reference pbmc68k tumor
dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/kotliarov/singleR/"



# Load query data
tiss_immune <- readRDS("N:/datasets/PBMC/kolitrov/kotliarov_pbmc_lognorm_sobj2.RDS")


# normalized query assay loaded
test_assay <- as.SingleCellExperiment(tiss_immune)

rm(tiss_immune)

# Load reference data
hpca.se=readRDS("N:/datasets/PBMC/pbmc68k/PBMC_68k.rds")
hpca.se=subset(hpca.se,celltype!="CD34+")
hpca.se=subset(hpca.se,celltype!="CD4+/CD25 T Reg")
hpca.se=subset(hpca.se,celltype!="CD4+ T Helper2")

table(hpca.se$celltype)

hpca.se$celltype=gsub("CD19+ B","B",hpca.se$celltype,fixed = T)
hpca.se$celltype=gsub("CD14+ Monocyte","Mono",hpca.se$celltype,fixed = T)
hpca.se$celltype=gsub("CD4+/CD45RA+/CD25- Naive T","TCD4 memory_naive",hpca.se$celltype,fixed = T)
hpca.se$celltype=gsub("CD4+/CD45RO+ Memory","TCD4 memory_naive",hpca.se$celltype,fixed = T)
hpca.se$celltype=gsub("CD56+ NK","NK",hpca.se$celltype,fixed = T)
hpca.se$celltype=gsub("CD8+ Cytotoxic T","TCD8",hpca.se$celltype,fixed = T)
hpca.se$celltype=gsub("CD8+/CD45RA+ Naive Cytotoxic","TCD8",hpca.se$celltype,fixed = T)
hpca.se$celltype=gsub("Dendritic","DC",hpca.se$celltype,fixed = T)
table(hpca.se$celltype)
# normalize reference
hpca.se <- NormalizeData(hpca.se, normalization.method = "LogNormalize", scale.factor = 10000)
# create single cell experiment for reference data
ref_sce <- as.SingleCellExperiment(hpca.se)

rm(hpca.se)


# prediction
library(SingleR)
start.time <- Sys.time()

pred.hesc <- SingleR(test = test_assay, ref = ref_sce,   assay.type.test = "logcounts",
                     assay.type.ref = "logcounts",
                     labels = ref_sce$celltype,fine.tune = F)
# save prediction
saveRDS(pred.hesc,file = paste0(dir_analysis,"singleR_prediction_ref_pbmc68k_kotliarov.RDS"))

end.time <- Sys.time() #41 min
time.taken <- end.time - start.time
time.taken

pred.hesc=readRDS(file = paste0(dir_analysis,"singleR_prediction_ref_pbmc68k_kotliarov.RDS"))
#accuracy harmonize
general_harmonnize_pred_stat(new_clusters=pred.hesc$pruned.labels,test_assay)
# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1       71.1 80.58 91.37       91.63     88.3    75.54

