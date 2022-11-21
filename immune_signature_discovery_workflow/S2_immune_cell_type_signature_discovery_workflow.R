
library(dplyr)
library(foreach)
library(Matrix)
library(data.table)
library(cluster)
library(ComplexHeatmap)
library(parallel)
library(kableExtra)
library(matrixStats)
library(uwot)
library("fpc")
library("dbscan")
library("factoextra")
library(gridExtra)
library(RColorBrewer)
library(gprofiler2)
library(dplyr)
library(ggplot2)
library(fst)
library(ggpubr)
library(reshape2)
library(patchwork)
library(Seurat)
library(varhandle)
library(ggrepel)
library(tibble)
library(varhandle)
set.seed(42)

# Script for discovery of signatures from 4 cancer scRNAseq datasets
# DBSCAN workflow for reproducing final text files and plots 

dirname="N:/integrate_4/v3/w_seurat_exp/" 
save_dir="N:./immunsig_scripts_manuscript/immune_signature_discovery_workflow/"

##################################################################
########### Load integrated expression matrix ####################
##################################################################

# Integrated expression matrix along with RData containing cell types
# from the discovery datasets and seperate expression matrices for each discovery
# dataset was obtained after Seurat integration script

# Load cell types from discovery datasets
load(paste0(dirname,"cell_types_list.RData"))

# Load seurat integrated expression matrix
all_exp= read.fst(paste0(dirname,"seurat_exp.fst"))
dim(all_exp)
# Features contain 10k variable genes, since in fst no rownames were allowed,
# these were set as rownames 
load(paste0(dirname,"features.RData"))
rownames(all_exp)=features

##################################################################
############### Density-based clustering #########################
##################################################################

# Dimension reduction using UMAP

umap_spec_100= umap(all_exp,pca=100,verbose = T,init = "agspectral")
save(umap_spec_100,file = paste0(dirname,"seurat_umap_spec_100.RData"))

load(file = paste0(dirname,"seurat_umap_spec_100.RData"))

um_1 =as.data.frame(umap_spec_100)
rownames(um_1)=features

# UMAP plot for 10k genes
# ggplot(as.data.frame(umap_spec_100),aes(umap_spec_100[,1],umap_spec_100[,2])) +
#   geom_point()+ xlab(label = "UMAP_1")+ylab(label = "UMAP_2")

# kNN plot was used to define epsilon for DBSCAN
# 0.075 was decided as epsilon

dbscan::kNNdistplot(um_1, k = 3)
abline(h = 0.075, lty = 2)

# DBSCAN
db <- fpc::dbscan(um_1, eps = 0.075, MinPts = 3)
#table(db$cluster)

#saving the clusters
cut_clusters=db$cluster
names(cut_clusters)=rownames(um_1)

#remove cluster 0 (noise)
cut_clusters= cut_clusters[cut_clusters!=0]
#table(cut_clusters)

##################################################################
################## Refining clusters #############################
##################################################################

# In order to refine the clusters silhouette coefficient-based filtering was applied 
# For the calculation of silhouette coefficent, 1-pearson correlation was used as 
# distance matrix
# Here pearson correlation was calculated for 10k genes based on integrated
# expression matrix

pearson_cor=cor(t(all_exp), method = "pearson")
load(paste0(dirname,"pearson_cor.RData")) 
#1-correlation(of average)
correlations = 1-pearson_cor

sil_cor = correlations[names(cut_clusters),names(cut_clusters)] %>%
  as.dist() %>%
  silhouette(cut_clusters,dist = .)

rownames(sil_cor) = names(cut_clusters)

#Cluster average silhoutte scores
sil_map= fviz_silhouette(sil_cor)

## silhoutte scores for genes less than 0 filtered out 
cut_clusters= cut_clusters[which(sil_cor[,3]>0)]

sil_cor = correlations[names(cut_clusters),names(cut_clusters)] %>%
  as.dist() %>%
  silhouette(cut_clusters,dist = .)

rownames(sil_cor) = names(cut_clusters)

#Cluster average silhoutte scores
sil_map= fviz_silhouette(sil_cor)

#av silhoutte less than 0 filter out
sil_cor_names=names(which(sil_map$plot_env$ave>0))

cut_clusters_filtered= cut_clusters%>% 
  subset(. %in%sil_cor_names)
#table(cut_clusters_filtered)

### Cluster filtering based on numbers in clusters

final_clusters= cut_clusters_filtered %>%
  data.frame(cluster_gr=.) %>%
  subset(cluster_gr %in%  names(which(table(.)>3)))
#table(final_clusters)

### If more than 50 genes, top 50 genes selected based on silhouette coefficient
more_clust= names(which(table(final_clusters)>=50))

if (length(more_clust) >0){
  for (i in seq_along(more_clust)) {
    sub_sil= final_clusters %>%
      subset(cluster_gr==more_clust[i]) %>%
      rownames(.) %>%
      sil_cor[.,]
    
    exc_genes= sort(sub_sil[,3],decreasing = T)
    exc_genes= names(exc_genes)[51:length(exc_genes)]
    
    final_clusters= final_clusters[!rownames(final_clusters) %in% exc_genes,,drop=F]
  }
}
#table(final_clusters)
final_clusters$cluster_gr=as.character(final_clusters$cluster_gr)

### Mean signature scores for each cell

final_sig_mean=
  foreach(j=unique(final_clusters$cluster_gr),.combine = cbind) %do%{
    final_clusters %>%
      subset(cluster_gr==j) %>%
      rownames() %>%
      all_exp[.,,drop=F] %>%
      colMeans() %>%
      as.data.frame() %>%
      setNames(nm = j)
  }

### Mean signature scores for each cell type
# Calculated averaging mean signature scores of cells belonging to individual cell types

cell_type_mean_sig=  
  foreach(i=seq_along(cell_types_list)) %do% {
    uniq_cell_type=unique(cell_types_list[[i]])
    uniq_cell_type=uniq_cell_type[order(uniq_cell_type)]
    foreach(j=seq_along(uniq_cell_type),.combine = cbind) %do%{
      type_pat = names(cell_types_list[[i]])[cell_types_list[[i]]==uniq_cell_type[j]]
      final_sig_mean[type_pat,] %>%
        colMeans(na.rm = T) %>%
        as.data.frame() %>%
        setNames(nm = uniq_cell_type[j])
    }
  }

# ID each cluster with "S_" following sequnetial number 
cluster_cell_type= paste0("S_",1:length(unique(final_clusters$cluster_gr)))
final_clusters_annot= final_clusters %>%
  mutate(annot=cluster_cell_type[match(.$cluster_gr, rownames(cell_type_mean_sig[[1]]))])
rownames(final_clusters_annot)=rownames(final_clusters)


#Filtering based on maximum-median
# For each signature, max-median difference calculated for each dataset
top_median=lapply(cell_type_mean_sig,FUN = function(x){
  round(rowMaxs(as.matrix(x))-rowMedians(as.matrix(x)),digits = 2)
  
})

top_median=Reduce(x = top_median,f = cbind)
storage.mode(top_median)="numeric"
rownames(top_median)=unique(final_clusters_annot$annot)

# Signatures are kept if max-median difference is greater than 1.2 in all discovery dataset
filtered= rownames(top_median[rowSums(top_median<1.2)!=4,])

filtered_final_clusters=final_clusters_annot %>%
  filter(annot %in% filtered)

write.table(filtered_final_clusters,paste0(save_dir,"dbscan_cluster_genes_param_0.tsv"),sep = "\t")
filtered_final_clusters=read.table(paste0(save_dir,"dbscan_cluster_genes_param_0.tsv"),sep = "\t")

# UMAP plot only for our selected genes

genes_found= rownames(um_1) %in% rownames(filtered_final_clusters) 
color= rownames(um_1)
color[!genes_found]=NA
color[genes_found]=filtered_final_clusters$annot

um_1$alpha=ifelse(genes_found,1,0.25)

# pdf(paste0(save_dir,"SuppFig2_filtered_genes_umap.pdf"),width = 11,height = 8)
# ggplot(as.data.frame(um_1),
#        aes(um_1[,1],um_1[,2],col=factor(color),label=factor(color))) +
#   geom_point()+ xlab( "UMAP_1")+ylab("UMAP_2")+ labs(color='Signatures') +
#   theme_bw(base_size = 24)
# dev.off()

# svg(paste0(save_dir,"SuppFig2_filtered_genes_umap.svg"),width = 11,height = 8)
# ggplot(as.data.frame(um_1),
#        aes(um_1[,1],um_1[,2],col=factor(color),label=factor(color))) +
#   geom_point(aes(alpha=alpha))+ xlab( "UMAP_1")+ylab("UMAP_2")+ labs(color='Signatures') +
#   theme_bw(base_size = 24)
# dev.off()

um_1$annot=color
mean_groups=data.frame(
  umap1_mean=um_1[genes_found,] %>% group_by(annot) %>%  
    summarise_at(vars(V1), list(name = mean)) %>% .$name,
  umap2_mean=um_1[genes_found,] %>% group_by(annot) %>%  
    summarise_at(vars(V2), list(name = mean))%>% .$name,
  annot=unique(filtered_final_clusters$annot)
)
svg(paste0(save_dir,"Fig2c_filtered_genes_umap.svg"),width = 11,height = 8)

ggplot() +
  geom_point(data = um_1, 
             mapping = aes(x = um_1[,1], 
                           y = um_1[,2], 
                           colour = annot,alpha=alpha),size=0.5) +
  # geom_point(mapping = aes_string(x = mean_groups$umap1_mean, 
  #                                 y = mean_groups$umap2_mean),
  #            color = "black", size = 1)+ xlab( "UMAP_1")+ylab("UMAP_2")+ labs(color='Signatures') +
  theme_bw(base_size = 24)+ theme(legend.position = "none") +
  geom_text_repel(mapping = aes_string(x = mean_groups$umap1_mean, 
                                       y = mean_groups$umap2_mean,
                                       label =mean_groups$annot)
                  , size = 7,    point.padding = 0, # additional padding around each point
                  min.segment.length = 0, # draw all line segments
                  max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
                  box.padding = 0.5 # additional padding around each text label
  ) + xlab( "UMAP_1")+ylab("UMAP_2") + theme(legend.position = "none")

dev.off()

# Also we constructed new UMAP only based on our genes in 25 clusters
all_exp_cluster_genes=all_exp[genes_found,]
umap_clustergenes= umap(all_exp_cluster_genes,pca=100,verbose = T,init = "agspectral")
save(umap_clustergenes,file = paste0(dirname,"seurat_umap_spec_100_cluster_genes.RData"))
load(file = paste0(dirname,"seurat_umap_spec_100_cluster_genes.RData"))

um_clustergenes =as.data.frame(umap_clustergenes)
rownames(um_clustergenes)=rownames(all_exp_cluster_genes)


um_clustergenes$annot=filtered_final_clusters$annot

mean_groups=data.frame(
  umap1_mean=um_clustergenes %>% group_by(annot) %>%  
    summarise_at(vars(V1), list(name = mean)) %>% .$name,
  umap2_mean=um_clustergenes %>% group_by(annot) %>%  
    summarise_at(vars(V2), list(name = mean))%>% .$name,
  annot=unique(um_clustergenes$annot)
)

svg(paste0(save_dir,"Fig2d_umap_refined_gene_clusters.svg"),width=11,height=8)

ggplot() +
  geom_point(data = um_clustergenes, 
             mapping = aes(x = um_clustergenes[,1], 
                           y = um_clustergenes[,2], 
                           colour = annot),size=0.5) +
  # geom_point(mapping = aes_string(x = mean_groups$umap1_mean, 
  #                                 y = mean_groups$umap2_mean),
  #            color = "black", size = 1)+ xlab( "UMAP_1")+ylab("UMAP_2")+ labs(color='Signatures') +
  theme_bw(base_size = 24)+ theme(legend.position = "none") +
  geom_text_repel( max.overlaps = Inf,mapping = aes_string(x = mean_groups$umap1_mean, 
                                       y = mean_groups$umap2_mean,
                                       label =mean_groups$annot)
                  , size = 6,    point.padding = 0, # additional padding around each point
                  min.segment.length = 0, # draw all line segments
                  max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
                  box.padding = 0.5 # additional padding around each text label
  ) + xlab( "UMAP_1")+ylab("UMAP_2")+ theme_bw(base_size = 24) + theme(legend.position = "none")

dev.off()




##################################################################
################## Defining clusters #############################
##################################################################
# GO analysis
unq_annot= unique(filtered_final_clusters$annot)

GO_table=list()

# Create and save plots for GO analysis
plots= foreach(i=seq_along(unq_annot)) %do% {
  print(i)
  
  clus = rownames(filtered_final_clusters)[filtered_final_clusters$annot==unq_annot[i]] %>%
    gost(query =.,
         organism = "hsapiens",significant = T,correction_method = "fdr",ordered_query = F,
         exclude_iea = T,user_threshold = 0.05)
  Go_tab= clus$result %>%
    select(term_name,term_size,intersection_size,p_value,precision)%>%
    arrange(p_value) %>%
    slice(1:20) #%>%
  GO_table[[i]]= Go_tab
  ggplot(Go_tab, aes(x=precision, y= reorder(term_name,precision), 
                     size=intersection_size, color=p_value)) + 
    theme(axis.text.y=element_blank(), axis.title.y=element_blank())+ ylab(label = "term name")+
    geom_point(alpha = 0.8)  + theme_classic() + ggtitle(label = unq_annot[i] )
}

names(GO_table)=unq_annot

glist <- lapply(plots, ggplotGrob)
ggsave(paste0(save_dir,"SuppFig1_db_GO_analysis.pdf"), width = 20,height = 10,scale=1,marrangeGrob(glist, nrow = 2, ncol = 2))

# Run the script literature_comparison_excel.RDS to search for literature
# It looks through gene sets from original publication (discovery datasets)
# --- finds which genes were assigned to which cell type


##### Annotate the signatures #############
###########################################
# Our annotations for our refined gene sets
table(filtered_final_clusters$annot)

Go_annot= c("Mast",
            "Antigen presentation",
            "Macrophage",
            "Plasma",
            "NK",
            "NK/TCD8",
            "B cell",
            "Monocyte",
            "Endothelial",
            "T CD4/B (IL12)",
            "DC immature",
            "CAF endothelial",
            "Cell cycle 1",
            "Myocyte",
            "CAF (integrin/TGFb)",
            "T CD4 memory/naive",
            "DC mature",
            "T CD8 memory 1",
            "pDC",
            "T CD4 reg",
          "Cell cycle 2",
            "Lymphoid cells",
            "Cell cycle/dna repair",
            "T CD8 memory 2",
            "Plasma (identical protein binding)"
)

# Save heatmap for discovery datasets as p_all
# It shows mean expression scores per cell types in each discovery dataset
# Later our annotations will be added to that and also validation dataset will be added
saveRDS(cell_type_mean_sig,paste0(save_dir,"celltypemeansig.RDS"))

  sorted_rownames=sort(Go_annot)

cell_type_mean_sig_sort=lapply(cell_type_mean_sig,function(x){
  rownames(x)=cluster_cell_type
  x=x[filtered,]
  x
  rownames(x)=Go_annot
  x[sorted_rownames,]
})

names(Go_annot)=unique(filtered_final_clusters$annot)
Go_annot=sort(Go_annot)
ha_names= rowAnnotation(gene_number= anno_text(which = "column",names(Go_annot),just="left"))

p_all=Heatmap(cell_type_mean_sig_sort[[1]], row_title = "Gene signatures",
               name="mean signature score (discovery)", column_title="GSE123139",
               row_names_side = "right",
               cluster_rows = F,cluster_columns=F,column_names_rot = 45,
               show_column_names=T,show_row_names = T,
               heatmap_legend_param = list(legend_direction = "horizontal",
                                           legend_width = unit(3, "cm"))
)+Heatmap(cell_type_mean_sig_sort[[2]], row_title = "Gene signatures",
           row_names_side = "right",
           name="value_2", column_title="GSE115978",column_names_rot = 45,
           cluster_rows = F,cluster_columns=F,show_heatmap_legend = F,
           show_column_names=T,show_row_names = T
)+Heatmap(cell_type_mean_sig_sort[[3]], row_title = "Gene signatures",
          name="value_3", column_title="GSE103322",column_names_rot = 45,
          cluster_rows = F,cluster_columns=F,row_names_side = "right",
          show_column_names=T,show_row_names = T,show_heatmap_legend = F
)+Heatmap(cell_type_mean_sig_sort[[4]], row_title = "Gene signatures",
          name="value_4", column_title="GSE140228",column_names_rot = 45,
          cluster_rows = F,cluster_columns=F,row_names_side = "right",
          show_column_names=T,show_row_names = T,show_heatmap_legend = F
)


# save final txt files containing refined signatures along with our annotations
unq_cluster= unique(filtered_final_clusters$cluster_gr)

Go_annot_df= data.frame(clust= unq_cluster, 
                        GO_annot=Go_annot)

final_clusters_GO_annot= filtered_final_clusters %>%
  mutate(Go=Go_annot_df[match(.$cluster_gr, Go_annot_df$clust),"GO_annot"])
# rownames(final_clusters_GO_annot)=rownames(filtered_final_clusters)
write.table(final_clusters_GO_annot,paste0(save_dir,"dbscan_cluster_genes_param_0_GO_annot.tsv"),sep = "\t")
final_clusters_GO_annot=read.table(paste0(save_dir,"dbscan_cluster_genes_param_0_GO_annot.tsv"),sep = "\t")

##################################################################
################## Validation dataset ############################
##################################################################
# Loading validation dataset
# Lung scRNAseq data from Zilionis et al. 2019

pbmc= readRDS("N:/datasets/myeloid/myeloid_raw_tumor_filtered_sobj.RDS")

# normalize, scale
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 10000)
pbmc <- ScaleData(pbmc,features =rownames(final_clusters_GO_annot))

### Gene signatures from our analysis (DBScan)
# calculate mean signature expression scores 

sig_list= foreach(i= seq_along(unq_annot)) %do% {
  tmp= final_clusters_GO_annot[final_clusters_GO_annot$annot==unq_annot[i],]
  rownames(tmp)
}
names(sig_list) =unq_annot

#intersection genes
inter_gen =foreach(i= seq_along(sig_list),.combine = rbind) %do%{
  data.frame(intersection=paste0(length(intersect(sig_list[[i]],rownames(pbmc))),"/",
                                 length(sig_list[[i]])))%>%
    setDF(rownames = names(sig_list)[i])
}

mean_sig =foreach(i= seq_along(sig_list),.combine = rbind) %do%{
  colMeans(as.matrix(pbmc[sig_list[[i]],]@assays[["RNA"]]@scale.data)) 
}
rownames(mean_sig) = unq_annot


cell_types_val= pbmc@meta.data[["cell.types"]]
names(cell_types_val)= colnames(pbmc)

uniq_cell_type=  unique(cell_types_val)

cell_type_mean_sig_valid= foreach(i=seq_along(uniq_cell_type),.combine = cbind) %do%{
  type_pat = names(cell_types_val)[cell_types_val==uniq_cell_type[i]]
  mean_sig[,type_pat] %>%
    rowMeans(na.rm = T) %>%
    as.data.frame() %>%
    setNames(nm = uniq_cell_type[i])
}

rownames(cell_type_mean_sig_valid)=unq_annot

#order the cell types alphabetically to have similar cells together
cell_type_mean_sig_valid=cell_type_mean_sig_valid[names(Go_annot),order(colnames(cell_type_mean_sig_valid),decreasing = F)]


#Intersected genes added as annotation
ha_n_valid= rowAnnotation(inter= anno_text(which = "column",inter_gen$intersection,just="left"))

rownames(cell_type_mean_sig_valid)=Go_annot
#pdf(paste0(dirname,"valid_lung_cancer.pdf"),width = 17,height = 9)
p_valid=Heatmap(cell_type_mean_sig_valid[sorted_rownames,],
                name="mean signature score (validation)", column_title="GSE127465",
                cluster_rows = F,cluster_columns=F,
                show_column_names=T,show_row_names = T,column_names_rot = 45,
                heatmap_legend_param = list(legend_direction = "horizontal",
                                            legend_width = unit(3, "cm"))
)
#draw(p_valid)
#dev.off()

svg(paste0(save_dir,"Fig2e_mean_score_our_signatures_heatmap.svg"),width = 22,height = 9)
a=ha_names+p_all+p_valid

draw(a,heatmap_legend_side = "top") 


dev.off()

#Table 2 

Genes=foreach(i=seq_along(unq_annot),.combine = rbind) %do% {
  b=rownames(final_clusters_GO_annot)[final_clusters_GO_annot$annot==unq_annot[i]]
  b=paste0(unique(b),collapse = ", ")
  data.frame(Genes=b)
}
rownames(Genes)=unfactor(unique(final_clusters_GO_annot$Go))


Genes=add_column(Genes, Annotation = rownames(Genes), .before = 1)
Genes=add_column(Genes, Cluster = unique(final_clusters_GO_annot$annot),
                 .before = 2)

gene_size=foreach(i=seq_along(unq_annot),.combine = c) %do% {
  b=rownames(final_clusters_GO_annot)[final_clusters_GO_annot$annot==unq_annot[i]]
  
  length(b)
}

Genes['Number of genes']=as.character(gene_size)

#selected signatures

order_selected=Go_annot= c(                           "B cell",
                                                      "DC immature",
                                                      "DC mature",
                                                      "Lymphoid cells",
                                                      "Macrophage",
                                                      "Mast",
                                                      "Monocyte",
                                                      "NK",
                                                      "NK/TCD8",
                                                      "pDC",
                                                      "Plasma",
                                                      "T CD4/B (IL12)",
                                                      "T CD4 memory/naive",
                                                      "T CD4 reg",
                                                      "CAF (integrin/TGFb)",
                                                      "CAF endothelial",
                                                      "Endothelial",
                                                      "Myocyte",
                                                      "Antigen presentation",
                                                      "Cell cycle 1",
                                                      "Cell cycle 2",
                                                      "Cell cycle/dna repair"
                                                      
                                                      
)

Genes=Genes[order_selected,]

write.table(Genes,
            paste0(save_dir,"Tab2_final_gene_list.tsv"),sep = "\t")



