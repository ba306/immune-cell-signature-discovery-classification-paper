# Cell type labeling of Kartha IFNg-PBMC dataset using Seurat, scType and
# RF+our gene sets
# Reference Hao
# Query Kartha IFNg-PBMC 


library(foreach)
library(Seurat)
library(qs)
set.seed(42)

########## Seurat and scType ############

# Input: integrated dataset for seurat and sctype

reference="./data/Hao/pbmc_hao_ref_up.qs"
reference <- qread(reference)
reference=reference[,reference$harmonized_celltype!="Other cells"]

query="./data/Kartha/Kartha_integrated.qs"

# Load query data
query <- qs::qread(file = query)

# Same parameters for all methods
# Colname for the celltypes in the reference
ref_col="harmonized_celltype"
# How to name the results 
ref_name="Hao"
query_name="Kartha"

# # Default  for singleR and chetah
# 2k HVGs

geneset="HVGs"
no_hvgs=2000

source("./functions/benchmarking_classifiers.R")

sctype_class()
seurat_class()

########## RF using our signatures ############

source("./functions/benchmarking_RF_training.R", local = TRUE, chdir = TRUE, keep.source = TRUE)

#### Our genes
our_genes= read.table("./discovery/dbscan_cluster_genes_param_GO_annot.tsv",sep = "\t",
                      stringsAsFactors = F,header = T)
table(our_genes$Annotation)

select=c(                    "Plasma"   ,
                             "NK"                   ,
                             "T CD8"              ,
                             "B"             ,
                             "Monocytes"               ,
                             "DC",
                             "T CD4"   ,
                             "pDC"                           
)

our_genes= our_genes %>%
  filter(!ID%in% c("S_6"))%>%
  subset(Annotation %in%select)
our_genes= our_genes$genes
length(our_genes) #167

# Input: non-integrated dataset 
query="./data/Kartha/Kartha_sobj.qs"

# Load query data
# For RF not integrated also works
# Input only to train the model using common genes!

query <- qs::qread(file = query)

RF_signatures(geneset="ourgenes",signature_genes = our_genes)
# [1] "ourgenes"
# [1] "Query dataset genes / Signature genes"
# [1] "164 / 167"
# [1] "Reference dataset genes / Signature genes"
# [1] "167 / 167"
# [1] "Common genes: 164"

# [1] "Training evaluation"
# [1] "B"     "DC"    "Mono"  "NK"    "T CD4" "T CD8"
# [1] "B"     "DC"    "Mono"  "NK"    "T CD4" "T CD8"
# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1       97.42 98.13 99.62        99.6    99.38    97.77
# [1] "Training test data evaluation"
# [1] "B"     "DC"    "Mono"  "NK"    "T CD4" "T CD8"
# [1] "B"     "DC"    "Mono"  "NK"    "T CD4" "T CD8"
# Sensitivity   PPV   NPV Specificity Accuracy f1_score
# 1       97.54 98.22 99.63       99.62    99.41    97.88