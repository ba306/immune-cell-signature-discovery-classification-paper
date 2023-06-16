# Benchmarking functions/wrappers for different cell type classification tools


# General Inputs 
# Reference dataset
# name of the reference
# Query dataset
# name of the query
# no of HVGs (if applicable) or in general no of genes
# prediction ouotput format:
# method_reference_query_geneset_nogenes

# take as gene space HVGs
# calculate in the reference
if(geneset=="HVGs"){
  
HVGs_grab= FindVariableFeatures(object = reference, 
                                selection.method = 'vst', nfeatures = no_hvgs)
genes=VariableFeatures(HVGs_grab)
# file extension
file_extension=paste0(ref_name, "_", query_name, "_", geneset, "_", length(genes))

# take all genes from the reference
}else if (geneset=="All"){
  genes=rownames(reference)
  # file extension
  file_extension=paste0(ref_name, "_", query_name, "_", geneset, "_", geneset)
  
#or e.g. in the case of geneset genes, take common genes from geneset and references
}else{
  genes=intersect(geneset_genes,rownames(reference))
  file_extension=paste0(ref_name, "_", query_name, "_", geneset, "_", length(genes))
  
}

print(paste0(length(genes)," genes used"))

# Define slot for the query
default_slot=DefaultAssay(query)

##### Seurat #####

seurat_class=function(){
  library(Seurat)
  
    reference <- ScaleData(reference,features =genes)
    sim.anchors <- FindTransferAnchors(reference = reference, query = query,
                                       dims = 1:30,features = genes)
    print("Classification starts")
    
  predictions <- TransferData(anchorset = sim.anchors, refdata = reference@meta.data[,ref_col],
                              dims = 1:30)
  predictions=predictions[,"predicted.id",drop=F]
  colnames(predictions)="Predictions"
  write.csv(predictions,file = 
    paste0("./benchmarking/prediction_results/seurat/Seurat_",file_extension,".csv")
  )
  
}

##### singleR #####

singleR_class=function(){
  library(Seurat)
  library(SingleR)
  
  # create single cell experiment for query data
  query_sce <- as.SingleCellExperiment(query)
  # create single cell experiment for reference data
  ref_sce <- as.SingleCellExperiment(reference)
  
  print("Classification starts")
  pred.hesc <- SingleR(test = query_sce, ref = ref_sce,   assay.type.test = "logcounts",
                         assay.type.ref = "logcounts",
                         labels = ref_sce[[ref_col]],fine.tune =T,restrict =genes )
  pred.hesc=pred.hesc[,"pruned.labels",drop=F]
  colnames(pred.hesc)="Predictions"
  rownames(pred.hesc)=rownames(pred.hesc)
  
  write.csv(pred.hesc,file = 
    paste0("./benchmarking/prediction_results/singleR/singleR_",file_extension,".csv")
  )
  
}

##### scType #####

# for sctype no reference is needed
sctype_class=function(){
  lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  # get cell-type-specific gene sets from our in-built database (DB)
  gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

  query <- ScaleData(query,features = genes)
  query <- RunPCA(query, features = genes)
  
  # cluster and visualize
  query <- FindNeighbors(query, dims = 1:10)
  query <- FindClusters(query, resolution = 0.8)
  query <- RunUMAP(query, dims = 1:10)
  
  # assign cell types
  print("Classification starts")
  
  es.max = sctype_score(scRNAseqData = query@assays[[default_slot]]@scale.data, scaled = T, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
  # In case Seurat is used, it is either reference[["RNA"]]@scale.data (default), query[["SCT"]]@scale.data, in case sctransform is used for normalization,
  # or query[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
  
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(query@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(query@meta.data[query@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(query@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  query@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    query@meta.data$customclassif[query@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  predictions=as.data.frame(query$customclassif)
  colnames(predictions)="Predictions"
  write.csv(predictions,file = 
    paste0("./benchmarking/prediction_results/scType/scType_",file_extension,".csv")
  )
  
}

##### CHETAH #####

chetah_class=function(){
  library(CHETAH)
  library(Seurat)
  
  reference=reference[genes,]
  
  ## Make SingleCellExperiments
  # normalized expression used!
  
  ref_sce <- SingleCellExperiment(assays = list(counts = reference@assays[[default_slot]]@data),
                                    colData = DataFrame(celltypes =reference@meta.data[,ref_col]))
  ###optional normalized for query
  query_sce <- SingleCellExperiment(assays = list(counts = query@assays[[default_slot]]@data))

  print("Classification starts")
  
  ## Run CHETAH, minimize unassign by set threshold as -Inf by thresh = -Inf, otherwise use default parameter
  if(geneset=="All"){
  # default n_genes is 200 !
    query_sce <- CHETAHclassifier(input = query_sce, ref_cells = ref_sce,thresh = -Inf,n_genes = 200 )
    ## Extract celltypes
    CHETAH_pred <- as.data.frame(query_sce$celltype_CHETAH)
    colnames(CHETAH_pred) <- 'Predictions'
    
    write.csv(CHETAH_pred,file = 
                paste0("./benchmarking/prediction_results/chetah/CHETAH_",
                       ref_name, "_", query_name,"_","default","_","200",
                       ".csv")
    )
  } else{
    query_sce <- CHETAHclassifier(input = query_sce, ref_cells = ref_sce,thresh = -Inf,n_genes = nrow(ref_sce) )
    ## Extract celltypes
    CHETAH_pred <- as.data.frame(query_sce$celltype_CHETAH)
    colnames(CHETAH_pred) <- 'Predictions'
    
    write.csv(CHETAH_pred,file = 
                paste0("./benchmarking/prediction_results/chetah/CHETAH_",file_extension,
                       ".csv")
    )
  }
  

}

