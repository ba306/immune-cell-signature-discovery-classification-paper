# Cell type labeling using Seurat, singleR, chetah, scType
# Reference Hao
# Query Zheng sorted

library(parallel)
library(foreach)
library(Seurat)
library(qs)
set.seed(42)


reference="./data/Hao/pbmc_hao_ref_up.qs"
query="./data/Zheng/pbmc_zheng_sorted_2k.qs"

# Load query data
query <- qread(query)

# Load reference data
reference <- qread(reference)
reference=reference[,reference$harmonized_celltype!="Other cells"]

# Same parameters for all methods
# Colname for the celltypes in the reference
ref_col="harmonized_celltype"
# How to name the results 
ref_name="Hao"
query_name="Zheng"

##### SingleR and chetah ##### 
# Default -> provide all genes
geneset="All"
genes=NA
no_hvgs=NA

#load the function to take the genes 
source("./functions/benchmarking_classifiers.R", local = TRUE, chdir = TRUE, keep.source = TRUE)

# Set the number of cores to use
num_cores <- 4 # Specify the number of cores you want to use

# Initialize the cluster with the desired number of cores
cl <- makeCluster(num_cores)

clusterExport(cl, c("reference", "no_hvgs","genes","file_extension","ref_col","query","geneset","genes",
                    "ref_name","query_name"),envir=environment())
functions_list <- list(
  singleR_class,
  chetah_class)

# Run each function in parallel
parLapply(cl, functions_list, function(f) f())

##### seurat and scType ##### 
# Run using varying HVGs
# Default is 2k HVGs

geneset="HVGs"

cl <- makeCluster(num_cores)

foreach(i=seq(100,2000,100)) %do%{
  no_hvgs=i
  print(i)
  #load the function to take appropriate HVGs
  source("./functions/benchmarking_classifiers.R",local = TRUE, chdir = TRUE, keep.source = TRUE)
  clusterExport(cl, c("reference", "no_hvgs","genes","file_extension","ref_col","query","geneset","genes","cl",
                      "ref_name","query_name"),envir=environment())
  functions_list <- list(seurat_class,
                         sctype_class)

  parLapply(cl, functions_list, function(f) f())
  
}

# Close the cluster
stopCluster(cl)
