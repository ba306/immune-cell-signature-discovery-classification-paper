library(Seurat)
library(foreach)
library(qs)

# Celltyping using CellTypist
############ Run the commands ############
# Define the Python script file and parameters
python_script <- "functions/benchmarking_celltypist.py"

# Load reference data
reference <- qread("./data/Hao/pbmc_hao_ref_up.qs")
reference=reference[,reference$harmonized_celltype!="Other cells"]

# Parameters 
# ref_name=sys.argv[1] # how to name the model "hao"
# reference=sys.argv[2] # name of the reference data with full directory './data/Hao/pbmc_hao_ref_up.h5ad'
# ref_col=sys.argv[3] # which columns are cell type metadata "harmonized_celltype"
# 
# query_name=sys.argv[4] # how to name the predicted query data "zheng"
# query=sc.read(sys.argv[5]) # # name of the query data with full directory "./data/Zheng/pbmc_zheng_sorted_2k.h5ad"
# 
# model_name=sys.argv[6] # if Immune_All_Low use immune_all_low model from the package if not it will print this name as model name 
# genes = sys.argv[7] # list of genes to use c("X,Y,Z")
# genes = genes.split(',') 
# geneset=sys.argv[8] # which genesets are used "HVGs" or gene signatures e.g. "Our_genes"


# Construct the command to run the Python script with parameters
# Training Hao with default setting and classifying query datasets

# Zheng benchmarking
command <- paste("python", python_script, 
                 "Hao", 
                 "./data/Hao/pbmc_hao_ref_up.h5ad",
                 "harmonized_celltype",
                 "Zheng",
                 "./data/Zheng/pbmc_zheng_sorted_2k.h5ad",
                 "default",
                 "default",
                 "default")
system(command)

# Kotliarov benchmarking

command <- paste("python", python_script, 
                 "Hao", 
                 "./data/Hao/pbmc_hao_ref_up.h5ad",
                 "harmonized_celltype",
                 "Kotliarov",
                 "./data/Kotliarov/kotliarov_pbmc.h5ad",
                 "default",
                 "default",
                 "default")
system(command)

# Kartha
command <- paste("python", python_script, 
                 "Hao", 
                 "./data/Hao/pbmc_hao_ref_up.h5ad",
                 "harmonized_celltype",
                 "Kartha",
                 "./data/Kartha/kartha.h5ad",
                 "default",
                 "default",
                 "default")
system(command)


# Training Hao with differing HVGs and classifying query datasets

# Kotliarov and Zheng benchmarking

foreach(i=seq(100,2000,100)) %do%{
  print(i)
  
HVGs_grab= FindVariableFeatures(object = reference, 
                                selection.method = 'vst', nfeatures = i)
genes=VariableFeatures(HVGs_grab)
genes=paste(genes, collapse = ",")

command <- paste("python", python_script, 
                 "Hao", 
                 "./data/Hao/pbmc_hao_ref_up.h5ad",
                 "harmonized_celltype",
                 "Zheng",
                 "./data/Zheng/pbmc_zheng_sorted_2k.h5ad",
                 "HVGs",
                 genes,
                 "HVGs")
system(command)

command <- paste("python", "functions/benchmarking_celltypist.py", 
                 "Hao", 
                 "./data/Hao/pbmc_hao_ref_up.h5ad",
                 "harmonized_celltype",
                 "Kotliarov",
                 "./data/Kotliarov/kotliarov_pbmc.h5ad",
                 "HVGs",
                 genes,
                 "HVGs")
# Execute the command
system(command)
}

