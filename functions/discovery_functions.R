
library(foreach)
library(dplyr)
library(ggplot2)
library("factoextra")
library(cluster)
library(tibble)
library(reshape2)
library(tidyr)
library(ComplexHeatmap)
library(ggrepel)
library(gridExtra)
library(gridGraphics)
library(grid)
library(patchwork)
library(ggpubr)
library(circlize)

# color palette 

chr=c("darkblue","brown4","darkgreen","darkorange","darkcyan",
      "darksalmon","deeppink","red2", "goldenrod1", "mediumorchid",
      "#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0","black")


####### Violin plot generation for cell type composition from different datasets #####

#input: list with all datasets
# only for major cell types 
# by="dataset" -> x axis dataset 
# by=major_celltypes -> x axis cell types
violin_celltype_dataperc=function(list_data,by="dataset"){
  
  plot_df=foreach(i=seq_along(list_data),.combine=rbind) %do% {
    data.frame(major_celltypes=list_data[[i]]@meta.data$major_celltypes,
               dataset=list_data[[i]]@meta.data$dataset)
    
  }
  #Percentual 
  
  if(by=="dataset"){
    prop_data <- plot_df %>%
      group_by(major_celltypes,dataset)%>%
      summarise(n = n()) %>%
      mutate(freq = 100*n / sum(n))
    
    print(
      ggplot(prop_data, aes(x = major_celltypes,y=freq, fill = dataset)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=chr) +
        labs(x = "Major celltypes", y = "Cell type percentages [%]", fill = "Dataset") +theme_bw()+
        theme(strip.text = element_text(size = 12), legend.direction = "horizontal",
              legend.position = "top",axis.text.x = element_text(angle = 20, hjust = 1))
      
    )
    
  }else if(by=="major_celltypes"){
    
    prop_data <- plot_df %>%
      group_by(dataset,major_celltypes)%>%
      summarise(n = n()) %>%
      mutate(freq = 100*n / sum(n))
    
    print(
      ggplot(prop_data, aes(x = dataset,y=freq, fill = major_celltypes)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=chr) +
        labs(x = "Dataset", y = "Cell type percentages [%]", fill = "Major celltypes")  +theme_bw()+
        theme(strip.text = element_text(size = 12), legend.direction = "horizontal",
              legend.position = "top"))+coord_flip()
  }
  
  
}

####### Various wrappers for cell cluster refinement and analysis ###############

#cut_clusters: cluster assignments- names gene sand values clusters
# refinement filtering gives dataframe with gene names and clusters
# sil_gene_threshold: silhouette score filter for genes
# min_no_gene= minimum genes for a cluster
# max_no_gene= select top n number of genes for each cluster

# take only clusters with min. 10 genes
# take only max 50 top silhouette scored genes

refinement_filtering=function(cut_clusters,
                              sil_gene_threshold,
                              min_no_gene,
                              max_no_gene){
  
  # Silhouette scores for  each gene
  
  sil_cor_gene = correlations[names(cut_clusters),names(cut_clusters)] %>%
    as.dist() %>%
    silhouette(cut_clusters,dist = .)
  
  rownames(sil_cor_gene) = names(sil_cor_gene)
  
  #Cluster average silhoutte scores
  sil_map_gene= fviz_silhouette(sil_cor_gene)
  
  ## silhoutte scores for genes less than 0 filtered out 
  cut_clusters_gene= cut_clusters[which(sil_cor_gene[,3]>sil_gene_threshold)]

  ### Cluster filtering based on numbers in clusters
  
  final_clusters_1= cut_clusters_gene %>%
    data.frame(cluster_gr=.) %>%
    subset(cluster_gr %in%  names(which(table(.)>=min_no_gene)))
  
  final_clusters_2=final_clusters_1
  
  ### If more than 50 genes, top 50 genes selected based on silhouette coefficient
  more_clust= names(which(table(final_clusters_1)>=max_no_gene))
  
  if (length(more_clust) >0){
    for (i in seq_along(more_clust)) {
      sub_sil= final_clusters_1 %>%
        subset(cluster_gr==more_clust[i]) %>%
        rownames(.) %>%
        sil_cor_clust[.,]
      
      exc_genes= sort(sub_sil[,3],decreasing = T)
      exc_genes= names(exc_genes)[51:length(exc_genes)]
      
      final_clusters_2= final_clusters_2[!rownames(final_clusters_2) %in% exc_genes,,drop=F]
    }
  }
  table(final_clusters_1)
  table(final_clusters_2)
  
  final_clusters_2$cluster_gr=as.character(final_clusters_2$cluster_gr)
  final_clusters_2$genes=rownames(final_clusters_2)
  
  # Result  summary 
  
  print("Pre-filtering")
  
  print(
    paste0("No_cluster: ",length(unique(cut_clusters)),"; No_genes: ", length(cut_clusters))
  )
  print("Gene-filtering")
  
  print(
    paste0("No_cluster: ",length(unique(cut_clusters_gene)),"; No_genes: ", length(cut_clusters_gene))
  )  
  print("Min gene-filtering")
  
  print(
    paste0("No_cluster: ",length(unique(final_clusters_1$cluster_gr)),
           "; No_genes: ", nrow(final_clusters_1))
  )
  
  print("Max gene-filtering")
  
  print(
    paste0("No_cluster: ",length(unique(final_clusters_2$cluster_gr)),
           "; No_genes: ", nrow(final_clusters_2))
  )
  
  final_clusters_2
}


# mean signature score per cell or per cell type

# cell_mode is either "celltype" (per cell type) or  "cell" (per cell)
# assay_mode either integrated or RNA to obtain scaled gene expression

mean_signature_scores=function(gene_list,data,cell_mode,assay_mode="integrated",
                               celltype_granularity="major_celltypes"){
  exp_data=data@assays[[assay_mode]]@scale.data
  gene_list$cluster_gr=as.factor(gene_list$cluster_gr)
  gene_list_matches=unique(gene_list[,!colnames(gene_list)%in%"genes",drop=F])
  
  # mean signature score per cell
  mean_score=function(gene_list,exp_data,data) {
    mean_cells=foreach(j=unique(gene_list$cluster_gr),.combine = cbind) %do%{
      gene_list %>%
        subset(cluster_gr==j) %>%
        rownames() %>%
        intersect(.,rownames(exp_data)) %>%
        exp_data[.,,drop=F] %>%
        colMeans() %>%
        as.data.frame() %>%
        setNames(nm = j) 
    }%>%
      rownames_to_column(var="cell")
    
    meta=data@meta.data[,c("dataset",celltype_granularity)] %>%
      rownames_to_column(var="cell")
    mean_melt=melt(mean_cells)
    # Merge the melted data frame with the second data frame by RowNames
    # merge metadata-cell types
    merged_df=merge(mean_melt, meta, by = "cell")
    colnames(merged_df)=c("cell","cluster_gr","score","dataset","celltypes")
    
    # merge and match signature list to the mean scores
    
    merge(merged_df,gene_list_matches)
  }
  if(cell_mode=="cell"){
    mean_score(gene_list,exp_data,data)
  }else if(cell_mode=="celltype"){
    # mean signature score per cell type
    # Calculated averaging mean signature scores of cells belonging to individual cell types
    
    mean_score(gene_list,exp_data,data)%>%
      group_by(dataset,celltypes,cluster_gr)  %>% 
      summarize(celltype_mean = mean(score)) %>%
      merge(.,gene_list_matches)
    
  }
}



# max-median filtering
# For each signature, max-median difference calculated for each dataset
# Signatures are kept if max-median difference is greater than threshold 
#in at least 3 discovery datasets

max_media_filter=function(gene_list,mean_score_df,
                          max_median_param,
                          min_sign_dataset){
  
  # mean signature score per cell type
  
  max_median_diff_df=mean_score_df %>%
    group_by(dataset,cluster_gr) %>%
    #maximum-median
    summarize(max_median_diff= round(max(celltype_mean)- median(celltype_mean),2))
  
  p=max_median_diff_df%>%
    ggplot( aes(x = max_median_diff)) + 
    geom_histogram(alpha = 0.5, position = "identity", bins = 20) + 
    labs(title = "Distribution Plot by dataset", x = "max_median_diff", y = "Frequency") +
    facet_wrap("dataset")+ geom_vline(xintercept = max_median_param, linetype = "dashed", color = "red")
  print(p)
  
  take=max_median_diff_df%>%
    filter(max_median_diff>max_median_param) %>%
    ungroup(dataset,cluster_gr) %>%
    dplyr::count(cluster_gr)%>%
    filter(n>=min_sign_dataset) %>% .$cluster_gr
  
  new_genelist=gene_list %>%
    filter(cluster_gr %in% take)
  
  print("Pre- max_median filtering")
  
  print(
    paste0("No_cluster: ",length(unique(gene_list$cluster_gr)),"; No_genes: ", nrow(gene_list))
  )
  print("max_median filtering")
  
  print(
    paste0("No_cluster: ",length(unique(new_genelist$cluster_gr)),"; No_genes: ", nrow(new_genelist))
  )
  
  new_genelist
}



####### UMAP creation for final gene sets and re-analysed UMAP #####


# UMAP plot only for our selected genes


umap_selected_genes=function(gene_list, # gene list
                             umap_df, #initial UMAP
                             reUMAP, # re-constructed UMAP based on final gene sets
                             label_var # based on which column from gene list to color and label 
                             ){
  
  #color our final genes
  
  umap_df$genes=rownames(umap_df)
  umap_df=left_join(umap_df,gene_list)
  
  # other genes gray ours with color and not faded
  umap_df$alpha=ifelse(is.na(umap_df$ID)==F,1,0.25)
  
  mean_groups_initial=data.frame(
    umap1_mean=umap_df %>% group_by(ID) %>%  
      summarise_at(vars(V1), list(name = mean)) %>% .$name,
    umap2_mean=umap_df%>% group_by(ID) %>%  
      summarise_at(vars(V2), list(name = mean))%>% .$name,
    ID=umap_df%>% group_by(ID) %>%  
      summarise_at(vars(V2), list(name = mean)) %>% .$ID
  ) %>%na.omit()

  umap_df_plotting=merge(mean_groups_initial,umap_df[!duplicated(umap_df$ID), ] %>%na.omit())


selected_genes_initialumap= ggplot() +
    geom_point(data = umap_df, 
               mapping = aes(x = umap_df[,1], 
                             y = umap_df[,2], 
                             colour =  !!as.name(label_var),alpha=alpha),size=0.5) +
    # geom_point(mapping = aes_string(x = mean_groups$umap1_mean, 
    #                                 y = mean_groups$umap2_mean),
    #            color = "black", size = 1)+ xlab( "UMAP_1")+ylab("UMAP_2")+ labs(color='Signatures') +
    theme_bw(base_size = 24)+ theme(legend.position = "none") +
    geom_text_repel(mapping = aes(x = umap_df_plotting$umap1_mean, 
                                  y = umap_df_plotting$umap2_mean,
                                  label =umap_df_plotting[,label_var])
                    , size = 4,    point.padding = 0, # additional padding around each point
                    min.segment.length = 0, # draw all line segments
                    max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
                    box.padding = 0.5 # additional padding around each text label
    ) + xlab( "UMAP1")+ylab("UMAP2") + theme(legend.position = "none")
  
  
  # Also we constructed new UMAP only based on our genes 
  # provided to the function
  um_clustergenes =as.data.frame(reUMAP)
  um_clustergenes$genes=rownames(um_clustergenes)
  
  um_clustergenes=merge(um_clustergenes,gene_list)
  
  mean_groups=data.frame(
    umap1_mean=um_clustergenes %>% group_by(ID) %>%  
      summarise_at(vars(V1), list(name = mean)) %>% .$name,
    umap2_mean=um_clustergenes %>% group_by(ID) %>%  
      summarise_at(vars(V2), list(name = mean))%>% .$name,
    ID=um_clustergenes%>% group_by(ID) %>%  
      summarise_at(vars(V2), list(name = mean)) %>% .$ID
  )
  
  umap_df_plotting_reumap=merge(mean_groups,um_clustergenes[!duplicated(um_clustergenes$ID), ] %>%na.omit())
  
  selected_genes_reumap=ggplot() +
    geom_point(data = um_clustergenes, 
               mapping = aes(x = um_clustergenes[,2], 
                             y = um_clustergenes[,3], 
                             colour = ID),size=0.5) +
    # geom_point(mapping = aes_string(x = mean_groups$umap1_mean, 
    #                                 y = mean_groups$umap2_mean),
    #            color = "black", size = 1)+ xlab( "UMAP_1")+ylab("UMAP_2")+ labs(color='Signatures') +
    theme_bw(base_size = 24)+ theme(legend.position = "none") +
    geom_text_repel( max.overlaps = Inf,mapping = aes(x = umap_df_plotting_reumap$umap1_mean, 
                                                             y = umap_df_plotting_reumap$umap2_mean,
                                                             label =umap_df_plotting_reumap[,label_var])
                     , size = 4,    point.padding = 0, # additional padding around each point
                     min.segment.length = 0, # draw all line segments
                     max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
                     box.padding = 0.5 # additional padding around each text label
    ) + xlab( "UMAP1")+ylab("UMAP2")+ theme_bw(base_size = 24) + theme(legend.position = "none")
  
  
  selected_genes_initialumap +selected_genes_reumap
}


######### Heatmap creation for gene signatures ###################

heatmap_meansig_scores_list=function(mean_sig_df, #mean signature score dataframe
                                     row_name_mode="cluster_gr", #row names 
                                row_annot=T, #show row names or not 
                                sort=T, # if T automatically sort row and column names
                                row_ord, #specific row order , character vector
                                col_ord # specific column order, character vector
                                # row_ord and col_ord should be set at the same time, if not doesnt work
                                ){
  
  unq_dataset=unique(mean_sig_df$dataset)
  
 foreach(data=unq_dataset) %do% {
    unq_celltype=mean_sig_df %>% 
      filter(dataset==data) 
    
    heatmap_data <- mean_sig_df %>% 
      filter(dataset==data) %>%
      pivot_wider(names_from = celltypes, values_from = celltype_mean)%>%
      ungroup()%>%
      select(-dataset)%>%
      column_to_rownames(var=row_name_mode)
    
    
    heatmap_data=heatmap_data[,  unique(unq_celltype$celltypes)]
    if(sort==T){
      heatmap_data=heatmap_data[sort(rownames(heatmap_data)),
                                sort(colnames(heatmap_data))]
    }else if(length(row_ord)>0 &length(col_ord)>0 ){
      heatmap_data=heatmap_data[row_ord,col_ord]
    }
   Heatmap(heatmap_data,
                 name="Mean signature score", column_title=data,
           row_names_side = "right", col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
           cluster_rows = F,cluster_columns=F,column_names_rot = 45,
                 show_column_names=T,show_row_names = row_annot,
                 heatmap_legend_param = list(legend_direction = "horizontal",
                                             legend_width = unit(3, "cm"))
    )

 }
}


violin_plots_meanscores_updated=function(gene_list,mean_sig_cell_df,annot_mode="ID_Annotation"){
  # Annotation
  unq_annot=unique(gene_list[,annot_mode])
  uniq_cell_type=unique(mean_sig_cell_df$celltypes)
  
  plot_list=foreach(i = seq_along(unq_annot))%do%{
    signature_celltype=gsub(pattern = "S_(.*)_(.*)",x = unq_annot[i],replacement = "\\2")
    pos=grep(paste0("^",signature_celltype,"$"),uniq_cell_type)
    
  if(signature_celltype=="Myeloid cells"){
      mean_sig_cell_df %>%
        mutate(celltypes=
                 case_when(
                   mean_sig_cell_df$celltypes=="DC" ~ "Myeloid cells",
                   mean_sig_cell_df$celltypes=="Monocytes" ~ "Myeloid cells",
                     mean_sig_cell_df$celltypes=="Macrophage" ~ "Myeloid cells",
                   T~mean_sig_cell_df$celltype
                 )
        )%>%
        filter(ID_Annotation==unq_annot[i]) %>%
        # Create violin plots
        ggplot(
          aes(x = reorder(celltypes, score, FUN = mean), y = score)) +
        geom_violin(trim = FALSE,scale = "width") +
        ggtitle(unq_annot[i])+
        theme_bw()+coord_flip() + theme(axis.title.x = element_blank(),
                                        axis.title.y = element_blank())+
        stat_compare_means(label = "p.signif", method = "wilcox",
                           ref.group = "Myeloid cells")
      
      
    }else if(signature_celltype=="Lymphoid cells"){
      mean_sig_cell_df %>%
        mutate(celltypes=
                 case_when(
                   mean_sig_cell_df$celltypes=="T CD4" ~ "Lymphoid cells",
                   mean_sig_cell_df$celltypes=="T CD8" ~ "Lymphoid cells",
                   T~mean_sig_cell_df$celltype
                 )
        )%>%
        filter(ID_Annotation==unq_annot[i]) %>%
        # Create violin plots
        ggplot(
          aes(x = reorder(celltypes, score, FUN = mean), y = score)) +
        geom_violin(trim = FALSE,scale = "width") +
        ggtitle(unq_annot[i])+
        theme_bw()+coord_flip() + theme(axis.title.x = element_blank(),
                                        axis.title.y = element_blank())+
        stat_compare_means(label = "p.signif", method = "wilcox",
                           ref.group = "Lymphoid cells")
      
      
    }else{
      
      mean_sig_cell_df %>%
        filter(ID_Annotation==unq_annot[i]) %>%
        # Create violin plots
        ggplot(
          aes(x = reorder(celltypes, score, FUN = mean), y = score)) +
        geom_violin(trim = FALSE,scale = "width") +
        ggtitle(unq_annot[i])+
        theme_bw()+coord_flip() + theme(axis.title.x = element_blank(),
                                        axis.title.y = element_blank())+
        stat_compare_means(label = "p.signif", method = "wilcox",
                           ref.group = uniq_cell_type[pos])
      
    } 
  }
  
  vln_arranged= grid.arrange(
    grobs = plot_list,
    ncol = 5,
    bottom = "Signature score",
    left = "Cell type",
    top = paste0("Dataset: ",unique(mean_sig_cell_df$dataset))
  )
  print(vln_arranged)
}


