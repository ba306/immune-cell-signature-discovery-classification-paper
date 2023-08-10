# External validation of the immune cell type signatures 
# Nieto TME immune cell atlas

library(Seurat)
library(dplyr)
library(varhandle)

# Functions for discovery workflow
source("./functions/discovery_functions.R")
# Analysis directory
dir="./discovery/"

######## Prepare validation dataset ########
#read immune cell type atlas 
# downloaded from https://zenodo.org/record/4263972
sobj=readRDS("./data/TME_datasets/TICAtlas_RNA.rds")

# remove liver 2 and melonama 2 since we used those in the discovery datasets
take=colnames(sobj)[!sobj$source %in%c("liver2","melanoma2")]
sobj=sobj[,take]

#only focus on clean annotated cell types
take=colnames(sobj)[!sobj$cell_type %in%c("Naive T cells",
                        "Proliferative T cells",
                  "Proliferative monocytes and macrophages")]
sobj=sobj[,take]

VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,pt.size=0)
sobj <- subset(sobj, subset = nFeature_RNA < 5000 & percent.mt < 15)

dim(sobj)
#[1]  87659 229753

# average number of genes per cell
mean(sobj@meta.data$nFeature_RNA)
# [1] 1265.42

# Meta data harmonization
meta=sobj@meta.data
meta$cell_type=unfactor(meta$cell_type)

sobj$major_celltypes=case_when(
  meta$cell_type== "B cells"~ "B",
  meta$cell_type== "Proliferative B cells"~ "B",
  meta$cell_type== "Plasma B cells"~ "Plasma",
  meta$cell_type== "Regulatory T cells"~ "T CD4",
  meta$cell_type== "T helper cells"~ "T CD4",
  meta$cell_type== "Th17 cells"~ "T CD4",
  meta$cell_type== "Proliferative T cells"~ "T CD4",
  meta$cell_type== "Recently activated CD4 T cells"~ "T CD4",
  meta$cell_type== "Naive-memory CD4 T cells"~ "T CD4",
  meta$cell_type== "Transitional memory CD4 T cells"~ "T CD4",
  meta$cell_type== "Pre-exhausted CD8 T cells"~ "T CD8",
  meta$cell_type== "Cytotoxic CD8 T cells"~ "T CD8",
  meta$cell_type== "Effector memory CD8 T cells"~ "T CD8",
  meta$cell_type== "Terminally exhausted CD8 T cells"~ "T CD8",
  meta$cell_type== "SPP1 TAMs"~ "Macrophage",
  meta$cell_type== "M2 TAMs"~ "Macrophage",
  meta$cell_type== "Proinflamatory TAMs"~ "Macrophage",
  meta$cell_type== "cDC"~ "DC",
  meta$cell_type== "mDC"~ "DC",
  meta$cell_type== "Mast cells"~ "Mast",
  T~meta$cell_type
  
)

sobj$minor_celltypes=case_when(
  meta$cell_type== "B cells"~ "B",
  meta$cell_type== "Proliferative B cells"~ "B proliferative",
  meta$cell_type== "Plasma B cells"~ "Plasma",
  meta$cell_type== "Regulatory T cells"~ "T CD4 reg",
  meta$cell_type== "T helper cells"~ "T CD4 Th",
  meta$cell_type== "Th17 cells"~ "T CD4 Th17",
  meta$cell_type== "Proliferative T cells"~ "T CD4 proliferative",
  meta$cell_type== "Recently activated CD4 T cells"~ "T CD4 activated",
  meta$cell_type== "Naive-memory CD4 T cells"~ "T CD4 naive/memory",
  meta$cell_type== "Transitional memory CD4 T cells"~ "T CD4 transitional memory",
  meta$cell_type== "Pre-exhausted CD8 T cells"~ "T CD8 pre-exhausted",
  meta$cell_type== "Cytotoxic CD8 T cells"~ "T CD8 cytotoxic",
  meta$cell_type== "Effector memory CD8 T cells"~ "T CD8 effector memory",
  meta$cell_type== "Terminally exhausted CD8 T cells"~ "T CD8 terminally exhausted",
  meta$cell_type== "SPP1 TAMs"~ "Macrophage SPP1 TAMs",
  meta$cell_type== "M2 TAMs"~ "Macrophage M2 TAMs",
  meta$cell_type== "Proinflamatory TAMs"~ "Macrophage proinf",
  meta$cell_type== "cDC"~ "DC cDC",
  meta$cell_type== "mDC"~ "DC mDC",
  meta$cell_type== "Mast cells"~ "Mast",
  T~meta$cell_type
  
)

####### Load signatures ##########
# cell type signatures
refined_df_max_median=read.table(paste0(dir,"dbscan_cluster_genes_param_0_annot.tsv"),
                                 sep = "\t",header = T)
rownames(refined_df_max_median)=refined_df_max_median$genes

# Exclude lineage signatures
refined_df_max_median=refined_df_max_median[!refined_df_max_median$Annotation %in% 
                                              c("Myeloid cells",
                                                "Lymphoid cells"),]
#focus on our genes
genes_take=unfactor(refined_df_max_median$genes)
genes_take_intersect=intersect(genes_take,rownames(sobj))
refined_df_max_median=unfactor(refined_df_max_median)

sobj=sobj[genes_take_intersect,]

# scale gene wise expression
sobj <- ScaleData(sobj,features =rownames(sobj))
sobj$dataset="TME_atlas_Nieto"

######## Enrichment violin plots ####################
# Violin plots with statistical testing for given signature and corresponding celltype
valid_mean_cell=mean_signature_scores(refined_df_max_median,
                                      sobj,
                                      "cell",
                                      "RNA",
                                      "major_celltypes") 
valid_mean_cell$ID_Annotation=paste0(valid_mean_cell$ID,
                                     "_",valid_mean_cell$Annotation)

refined_df_max_median$ID_Annotation=paste0(refined_df_max_median$ID,"_",
                                           refined_df_max_median$Annotation)


pdf(paste0(dir,"violin_plots_validation.pdf"),width = 15,height = 10)

refined_df_max_median %>%
  violin_plots_meanscores_updated(gene_list=.,
                                  mean_sig_cell_df=valid_mean_cell,
                                  annot_mode="ID_Annotation")

dev.off()

######## Heatmap for the signatures ####################

# Calculate mean signature score per celltype
valid_mean=refined_df_max_median %>%
mean_signature_scores(.,
                                 sobj,
                                 "celltype",
                                 "RNA",
                                 "major_celltypes") 

valid_mean$Annotation_ID=paste0(valid_mean$Annotation," (",valid_mean$ID,")")

png(paste0(dir,"heatmaps_validation.png"),width = 800,
    height=500)
heatmap_meansig_scores(valid_mean,row_name_mode = "Annotation_ID",
                       sort=T,scale=F)
dev.off()

######## Final gene lists ####################

# Final table summarizing cluster number, cluster annotation
# genes in the cluster and size of the clusters
# only selected final genes are included
# remove S_6 since couldnt be validated
remove= 
  c(
    "S_6")
refined_df_max_median=refined_df_max_median[!refined_df_max_median$ID %in% remove,]

write.table(refined_df_max_median,
  paste0(dir,"immdisc_aybey_final_list.tsv"),sep = "\t",row.names = F)

# summarize the data frame by ID and Annotation, and combine the genes values separated by comma
refined_df_max_median %>%
  group_by(ID, Annotation) %>%
  summarize(genes = paste(genes, collapse = ", "),
            num_genes = n()) %>%
  ungroup()%>%
  arrange(Annotation) %>%
  rename('Cluster number' = ID,'Cluster annotation' = Annotation,
         Genes=genes,'Number of genes'=num_genes) %>%
  write.table(
    paste0(dir,"final_gene_list.tsv"),sep = "\t",row.names = F)




