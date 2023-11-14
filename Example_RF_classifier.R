# Examplary script to classify cell types using Random Forest RF classifier based on cell type specific gene signatures

# Cell type labeling using RF trained on our immune cell type genes
# Gene list: our immune cell type genes 
# Reference: Hao
# Query: Zheng sorted 

library(foreach)
library(qs)
library(dplyr)
set.seed(42)


# Call the function for the training and classifier
source("./functions/RF_celltype_classifier.R")
source("./functions/benchmarking_prediction_metrics.R")

#### Define labeling of the outpus ####

# How to name the results 

ref_name="Example_Hao"
query_name="Zheng"
save_dir= "./"

# Colname for the celltype metadata in the reference

ref_col="harmonized_celltype"

##### Load reference and query datasets ####
reference="./data/Hao/pbmc_hao_ref_up.qs"
query="./data/Zheng/pbmc_zheng_sorted_2k.qs"

# Load query data
query <- qread(query)

# Load reference data
reference <- qread(reference)
reference=reference[,reference$harmonized_celltype!="Other cells"]

########## Prepare the gene signatures ############
#### Our genes ####
our_genes= read.table("./discovery/immdisc_aybey_final_list.tsv",sep = "\t",
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
  subset(Annotation %in%select)
our_genes= our_genes$genes
length(our_genes) #167

print("Geneset ready")

########## RF training and prediction ############

print("RF training and prediction starts")

RF_signatures(geneset="ourgenes",
              signature_genes = our_genes,
              assay = "RNA",
              reference =reference ,
              ref_col=ref_col,
              query =query,
              ref_name=ref_name,
              query_name=query_name,
              save_dir=save_dir)

print("Prediction model and results are saved")
print("RF training and prediction ends")

