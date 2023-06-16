# Cell type labeling using RF trained on different immune cell type gene sets
# Reference Hao
# Query Kotliarov

library(parallel)  
library(foreach)
library(readxl)
source("./functions/benchmarking_RF_training.R", local = TRUE, chdir = TRUE, keep.source = TRUE)
set.seed(42)

reference="./data/Hao/pbmc_hao_ref_up.qs"
query="./data/Kotliarov/kotliarov_pbmc.qs"

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
query_name="Kotliarov"

########## Prepare the gene signatures ############
#### Our genes ####
our_genes= read.table("./discovery/dbscan_cluster_genes_param_0_annot.tsv",sep = "\t",
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

#### Abbas genes ####
Abbas= read_excel(path ="./gene_lists_papers/abbas_2005.xlsx",sheet = 1)
# Select relevant immune cell type signatures
Abbas= Abbas %>%
  subset(Cell_type %in%unique(Abbas$Cell_type)[c(3,6,7,8,9)] )
#[1] "Dendritic Cell" "Monocyte"       "B Cell"         "T Cell"         "NK Cell"
abbas_genes=unique(Abbas$Gene)


#### Charoentong genes ####
Charoentong= read_excel(path = paste0("./gene_lists_papers/",
                                      "/Charoentong.xlsx"),sheet = 1)
Charoentong_f= foreach(i=seq_along(rownames(Charoentong)),.combine = rbind) %do% {
  genes= Charoentong[i,-c(1,2)]
  genes=genes[!is.na(genes)]
  cell_type= Charoentong[i,1]
  data.frame(Cell_type=cell_type,Gene=genes)
}
Charoentong_f=cbind(Charoentong_f,data.frame(source="Charoentong_2017"))
colnames(Charoentong_f)[1]="Cell_type"

# Select relevant immune cell type signatures
select=c(
  "B-cells" ,       "CD4+ T-cells" ,  "CD4+ Tcm"    ,   "CD4+ Tem"  ,    
  "CD8+ T-cells"   ,"CD8+ Tcm"  ,     "CD8+ Tem"     ,      "iDC"   ,        
  "Memory B-cells", "Monocytes" ,      
  "NK cells"       ,        "pDC" ,          
  "Tregs"  
  
)
Charoentong_f= Charoentong_f %>% 
  subset(Cell_type %in%select )
Charoentong_genes=unique(Charoentong_f$Gene)

#### Angelova genes ####
published=read.table(paste0("./gene_lists_papers/","published_sign.txt"),sep = "\t",
                     stringsAsFactors = F,header=T)
angelova=published%>%
  filter(meta_publication %in% c("https://doi.org/10.1186/s13059-015-0620-6")) %>%
  select(id_geneset_name,gene_symbol)

cell_types=gsub(pattern = "angelova2015[.]",x = angelova$id_geneset_name,replacement = "")
angelova_gene=data.frame(Cell_type=cell_types,Gene=angelova$gene_symbol,source="Angelova_2015")
# Select relevant immune cell type signatures
select=c( "Central.memory.CD4" ,
          "Central.memory.CD8" ,
          "DC"                ,
          "Effector.memory.CD4",
          "Effector.memory.CD8",
          "Immature.B.cells"  ,        
          "Memory.B.cells"   ,
          "Monocytes"         , 
          "NK"         ,         "NK56.bright"    ,     "NK56.dim"   ,                
          "Treg",               
          "iDC" ,                "mDC"      ,           "pDC"    )  

angelova_gene= angelova_gene %>%
  subset(Cell_type %in%select )

angelova_genes=unique(angelova_gene$Gene)

#### Nieto genes ####
Nieto= read_excel(path = paste0("./gene_lists_papers/",
                                "/nieto_immune_2021.xlsx"),sheet = 1)
selected=c("B cells"       ,  "Plasma B cells"   ,    
           "T cells naive"   ,"T cells regulatory", 
           "CD4 effector memory"    ,"CD4 transitional memory" ,
           "CD4 naïve-memory"    ,       "CD8 effector memory"  ,"NK"     ,"CD8 cytotoxic",
           "Monocytes"      ,          "cDC"        ,              "pDC"      ,                "mDC" )
nieto_genes=Nieto[Nieto$`Cell type` %in% selected,] %>%.$Markers %>%
  strsplit(", ") %>%unlist() %>% unique()

print("Genesets ready")

########## RF training and prediction ############
print("RF training and prediction starts")

# Set the number of cores to use
num_cores <- 6 # Specify the number of cores you want to use

# Initialize the cluster with the desired number of cores
cl <- makeCluster(num_cores)

# Load required packages and functions on the worker nodes
clusterEvalQ(cl, {
  library(Seurat)
  source("./functions/benchmarking_RF_training_prediction_metrics.R", local = TRUE, chdir = TRUE, keep.source = TRUE)
  print <- function(x) { base::print(x) }
})

# Define the parameter combinations
params <- list(
  list(func = RF_signatures,geneset="ourgenes",signature_genes = our_genes),
  list(func = RF_signatures,geneset="random167genes",signature_genes = 167),
  list(func = RF_signatures,geneset="abbas",signature_genes = abbas_genes),
  list(func = RF_signatures,geneset="charoentong",signature_genes = Charoentong_genes),
  list(func = RF_signatures,geneset="angelova",signature_genes = angelova_genes),
  list(func = RF_signatures,geneset="nieto",signature_genes = nieto_genes)
)

clusterExport(cl, c("reference", "query","ref_name","query_name","ref_col"),envir=environment())

# Submit each function as an individual job to the cluster
results <- parLapply(cl, params, function(x) {
  captured_output <- capture.output({
    do.call(x$func, list(x$geneset, x$signature_genes))
  })
  paste(captured_output, collapse = "")  # Concatenate output lines into a single line
  
})

# Close the cluster
stopCluster(cl)


# Save the captured output from all functions to a single text file
output_file <- "./benchmarking/models/RF_outputs_kotliarov.txt"
file_conn <- file(output_file, open = "w")

for (i in seq_along(results)) {
  cat("Output from Function", i, ":\n", file = file_conn)
  cat(results[[i]], "\n", file = file_conn, append = TRUE)
}

close(file_conn)
