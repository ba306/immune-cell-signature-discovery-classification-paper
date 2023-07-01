# Script for discovery of signatures from TME scRNAseq datasets
# DBSCAN workflow for reproducing final text files and plots 
# Input: scaled integrated expression matrix
library(Seurat)
library(qs)
library(uwot)
library(ggplot2)
library(varhandle)
library(reshape2)
library(tibble)
library(fpc)
library(dbscan)
set.seed(42)

# Functions for discovery worfklow
source("./functions/discovery_functions.R")

# Directory
dir="./discovery/"

########### Load and create inputs ####################

# Integrated seurat object containing cell types
# dataset was obtained after Seurat integration script
integrated=qread( paste0(dir,"integrated_discovery.qs"))
# Scaled expression matrix for DBSCAN
all_exp=integrated@assays$integrated@scale.data

# Scaled matrices for each individual datasets 
data_list <- qs::qread("./data/TME/discovery_datalist.qs")
# Scale each invidivual datasets individually
data_list <- lapply(data_list, function(x) {
  x <- ScaleData(x, features = rownames(all_exp), verbose = T)
  
})

# In order to refine the clusters silhouette coefficient-based filtering was applied 
# For the calculation of silhouette coefficent, 1-pearson correlation was used as 
# distance matrix
# Here pearson correlation was calculated based on integrated expression matrix

pearson_cor=cor(t(all_exp), method = "pearson")

#1-correlation
correlations = 1-pearson_cor
# Correlation as dist for silhouette calculations
cor_dist=correlations %>%
  as.dist()


############### Density-based clustering #########################

# Dimension reduction using UMAP

umap_spec_100= uwot::umap(all_exp,verbose = T,pca=20)

um_1 =as.data.frame(umap_spec_100)

# UMAP plot 
ggplot(as.data.frame(umap_spec_100),aes(umap_spec_100[,1],umap_spec_100[,2])) +
  geom_point()+ xlab(label = "UMAP_1")+ylab(label = "UMAP_2")

# kNN plot  for DBSCAN
# minPts=10 aim for min. 10 genes per cluster
dbscan::kNNdistplot(um_1, k = 10)
# Optimal seems to be around 0.25-0.3
abline(h = 0.3, lty = 2)

# But we want to have more clusters and higher resolution so we will use a lower epsilon 
# We tried a range between 0.15 and 0.25 by 0.01 interval to find the better resolution
# We ran the gene refinement workflow from the following section incl. 
# filtering out genes or cluster with negative silhouette scores and taking 
# only top 50 genes with highest sillhouette scores 
# Then we examined the heatmaps showing the mean signature scores for each cell type in each discovery cohort
# Epsilon 0.18 captured all cell types -> higher resolution

pdf(paste0(dir,"epsilon_searchheatmap_15_20_by0_01_minpts10.pdf"),width = 30,height = 10)
for(e in seq(0.15, 0.20, by = 0.005)){
  print(e)
  # DBSCAN
  db <- fpc::dbscan(um_1, eps = e, MinPts =10)

  #saving the clusters
  cut_clusters=db$cluster
  names(cut_clusters)=rownames(um_1)
  
  #remove cluster 0 (noise)
  cut_clusters= cut_clusters[cut_clusters!=0]

  # Refining clusters and analysis
  # Refinement helper function with given parameters
  # filter out negative silhouette scored- genes
  # take only clusters with min. 10 genes
  # take only max 50 top silhouette scored genes
  sil_gene_threshold=0
  min_no_gene=10
  max_no_gene=50
  
  refined_df=refinement_filtering(cut_clusters,sil_gene_threshold,
                                  min_no_gene,max_no_gene    )
  
  ### Mean signature scores for each cell type
  
  mean_celltypedf=foreach(i=seq_along(data_list),.combine=rbind)%do%{
    print(i)
    mean_signature_scores(refined_df,
                          data_list[[i]],
                          "celltype",
                          "RNA",
                          "major_celltypes")
  }
  
  all_plots=
    heatmap_meansig_scores_list(mean_celltypedf,
                                row_name_mode = "cluster_gr",
                                row_annot = F,ncol = 7)
  
  
  combined_heatmap <- Reduce("+", all_plots)
  
  draw(combined_heatmap,heatmap_legend_side="top")
  
}
dev.off()

# DBSCAN
db <- fpc::dbscan(um_1, eps = 0.18, MinPts = 10)

#saving the clusters
cut_clusters=db$cluster
names(cut_clusters)=rownames(um_1)

#remove cluster 0 (noise)
cut_clusters= cut_clusters[cut_clusters!=0]

################## Refining clusters and analysis #############################

# Refinement helper function with given parameters
sil_gene_threshold=0
min_no_gene=10
max_no_gene=50

refined_df=refinement_filtering(cut_clusters,sil_gene_threshold,
                                min_no_gene,max_no_gene    )
# [1] "Pre-filtering"
# [1] "No_cluster: 57; No_genes: 2737"
# [1] "Gene-filtering"
# [1] "No_cluster: 30; No_genes: 634"
# [1] "Min gene-filtering"
# [1] "No_cluster: 26; No_genes: 608"
# [1] "Max gene-filtering"
# [1] "No_cluster: 26; No_genes: 602"

### Mean signature scores for each cell type

mean_celltypedf=foreach(i=seq_along(data_list),.combine=rbind)%do%{
  print(i)
  mean_signature_scores(refined_df,
                        data_list[[i]],
                        "celltype",
                        "RNA",
                        "major_celltypes")
}

#Filtering based on maximum-median
# For each signature, max-median difference calculated for each dataset
# Signatures are kept if max-median difference is greater than threshold 
#in at least 3 discovery datasets

refined_df_max_median=max_media_filter(refined_df,mean_celltypedf,
                                       max_median_param=0.6,
                                       min_sign_dataset=3  )

# [1] "Pre- max_median filtering"
# [1] "No_cluster: 26; No_genes: 602"
# [1] "max_median filtering"
# [1] "No_cluster: 14; No_genes: 338"

# ID each cluster with "S_" following sequnetial number 
# Determine the unique values in the vector
unique_nums <- unique(refined_df_max_median$cluster_gr)
refined_df_max_median$ID<- paste0("S_", match(refined_df_max_median$cluster_gr, 
                                              unique_nums))

# Umap plots showing our gene sets in the initial and reconstructed UMAP
# Also we constructed new UMAP only based on our genes 
all_exp_cluster_genes=all_exp[refined_df_max_median$genes,]

umap_clustergenes= umap(all_exp_cluster_genes,pca=20,verbose = T)
umap_clustergenes=as.data.frame(umap_clustergenes)

png(paste0(dir,"umap_selected_genes.png"), 
    width = 35, height = 20, units = "cm", res = 300)
umap_selected_genes(refined_df_max_median,
                    um_1,
                    umap_clustergenes,"ID")
dev.off()


######### Annotate the signatures #############
# Our annotations for our refined gene sets

Annot_df= data.frame(ID=sort(unique(refined_df_max_median$ID)),Annotation=
                       c("Mast",
                         "Macrophage",
                         "T CD4",
                         "Plasma",
                         "B",
                         "DC",
                         "T CD8",
                         "Myeloid cells",
                         "NK",
                         "Macrophage",
                         "Monocytes",
                         "Lymphoid cells",
                         "Monocytes",
                         "pDC"
                         
                         
                       )
)
# Heatmaps for each signature in each discovery dataset 
# Major cell type annotations
# Mean signature scores are drawn

mean_celltypedf_max_median=foreach(i=seq_along(data_list),.combine=rbind)%do%{
  print(i)
  mean_signature_scores(refined_df_max_median,
                        data_list[[i]],
                        "celltype",
                        "RNA",
                        "major_celltypes")
}

# Match annotation to the cluster number
mean_celltypedf_annot=mean_celltypedf_max_median
mean_celltypedf_annot=merge(mean_celltypedf_annot,Annot_df)
refined_df_max_median_annot=merge(refined_df_max_median,Annot_df)

# save final txt files containing refined signatures along with our annotations

write.table(refined_df_max_median_annot,paste0(dir,"dbscan_cluster_genes_param_0_annot.tsv"),
            sep = "\t",row.names = F)
rownames(refined_df_max_median_annot)=refined_df_max_median_annot$genes

#######  Final heatmaps ####### 
# left side Cluster number
# right side signature names

# cell type and cluster ID together
mean_celltypedf_annot$Annotation_ID=paste0(mean_celltypedf_annot$Annotation," (",mean_celltypedf_annot$ID,")")

all_plots=c(
  # these plots are without row annotation labels 
  mean_celltypedf_annot%>%
    filter(dataset!=unique(mean_celltypedf_annot$dataset)[7])%>%
    heatmap_meansig_scores_list(.,
                                row_name_mode = "Annotation_ID",
                                row_annot = F),
  # this is plotted to have row annotations right next to the heatmap
  mean_celltypedf_annot%>%
    filter(dataset==unique(mean_celltypedf_annot$dataset)[7])%>%
    heatmap_meansig_scores_list(.,
                                row_name_mode = "Annotation_ID",
                                row_annot = T)
  
)

combined_heatmap <- Reduce("+", all_plots)

png(paste0(dir,"heatmaps_discovery.png"), 
    width = 55, height = 15, units = "cm", res = 300)
draw(combined_heatmap,heatmap_legend_side="top")
dev.off()

####### Signature enrichments: Violin plots ########
# violin plots to test enrichment of each signature in the corresponding cell type 
# each dataset separately

# calculate mean signature score for each cell type 
disc_mean_cell=foreach(i=seq_along(data_list),.combine=rbind)%do%{
  print(i)
  mean_signature_scores(refined_df_max_median_annot,
                        data_list[[i]],
                        "cell",
                        "RNA",
                        "major_celltypes")
}

#put ID and cell type together e.g. S_1_Mast
disc_mean_cell$ID_Annotation=paste0(disc_mean_cell$ID,"_",disc_mean_cell$Annotation)
refined_df_max_median_annot$ID_Annotation=paste0(refined_df_max_median_annot$ID,"_",
                                                 refined_df_max_median_annot$Annotation)


pdf(paste0(dir,"violin_plots_discovery.pdf"),width = 20,height = 10)
for(i in unique(disc_mean_cell$dataset)){
  print(i)
  disc_mean_cell %>%
    filter(dataset==i) %>%
    violin_plots_meanscores_updated(gene_list=refined_df_max_median_annot,
                                    mean_sig_cell_df=.,
                                    annot_mode="ID_Annotation")
}
dev.off()

