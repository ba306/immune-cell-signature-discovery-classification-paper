
library(Seurat)
library(dplyr)
set.seed(42)

source("N:./immunsig_scripts_manuscript/benchmarking/kotliarov/kotliarov_harmonize_prediction_scores_func.R")

# reference zilionis tumor
dir_analysis="N:./immunsig_scripts_manuscript/benchmarking/kotliarov/singleR/"


# Load query data
tiss_immune <- readRDS("N:/datasets/PBMC/kolitrov/kotliarov_pbmc_lognorm_sobj2.RDS")

# normalized query assay loaded
test_assay <- as.SingleCellExperiment(tiss_immune)

rm(tiss_immune)

# Load reference data
# Lung scRNAseq data from Zilionis et al. 2019
hpca.se <- readRDS("N:/datasets/myeloid/myeloid_raw_tumor_filtered_sobj.RDS")
table(hpca.se$cell.types)
hpca.se=hpca.se[,hpca.se$cell.types!="Eosinophils"]
hpca.se=hpca.se[,hpca.se$cell.types!="Neutrophils"]
hpca.se=hpca.se[,hpca.se$cell.types!="Mast cells activated"]
hpca.se=hpca.se[,hpca.se$cell.types!="Mast cells resting"]
hpca.se=hpca.se[,hpca.se$cell.types!="Macrophages M0"]
hpca.se=hpca.se[,hpca.se$cell.types!="Macrophages M1"]
hpca.se=hpca.se[,hpca.se$cell.types!="Macrophages M2"]

hpca.se=hpca.se[,hpca.se$cell.types !="T cells follicular helper" &
                  hpca.se$cell.types !="T cells regulatory (Tregs)" ]
hpca.se$cell.types=gsub("B cells memory","B",hpca.se$cell.types)
hpca.se$cell.types=gsub("B cells naive","B",hpca.se$cell.types)
hpca.se$cell.types=gsub("Plasma cells","B",hpca.se$cell.types)
hpca.se$cell.types=gsub("Dendritic cells activated","DC",hpca.se$cell.types)
hpca.se$cell.types=gsub("Dendritic cells resting","DC",hpca.se$cell.types)
hpca.se$cell.types=gsub("NK cells activated","NK",hpca.se$cell.types)
hpca.se$cell.types=gsub("NK cells resting","NK",hpca.se$cell.types)
hpca.se$cell.types=gsub("T cells CD4 naive","T CD4 naive_memory",hpca.se$cell.types)
hpca.se$cell.types=gsub("T cells CD4 memory resting","T CD4 naive_memory",hpca.se$cell.types)
hpca.se$cell.types=gsub("T cells CD4 memory activated","T CD4 naive_memory",hpca.se$cell.types)
hpca.se$cell.types=gsub("T cells CD8","T CD8",hpca.se$cell.types)
hpca.se$cell.types=gsub("Monocytes","Mono",hpca.se$cell.types)

table(hpca.se$cell.types)

# normalize reference
hpca.se <- NormalizeData(hpca.se, normalization.method = "LogNormalize", scale.factor = 10000)
# create single cell experiment for reference data
ref_sce <- as.SingleCellExperiment(hpca.se)

rm(hpca.se)


# prediction
library(SingleR)
start.time <- Sys.time()

pred.hesc <- SingleR(test = test_assay, ref = ref_sce, assay.type.test = "logcounts",
                     assay.type.ref = "logcounts",
                     labels = ref_sce$cell.types,fine.tune = F)
# save prediction
saveRDS(pred.hesc,file = paste0(dir_analysis,"singleR_prediction_ref_zil_tumor_kotliarov.RDS"))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken #3min?

pred.hesc=readRDS(file = paste0(dir_analysis,"singleR_prediction_ref_zil_tumor_kotliarov.RDS"))
#accuracy harmonize

general_harmonnize_pred_stat(new_clusters=pred.hesc$pruned.labels,test_assay)
# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1      40.54 58.8 80.11       85.26     74.1    47.99