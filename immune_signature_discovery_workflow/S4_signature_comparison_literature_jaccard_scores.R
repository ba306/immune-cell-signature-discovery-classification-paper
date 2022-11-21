
library(foreach)
library(ComplexHeatmap)
library(bayesbio)
library(readxl)
library(dplyr)
library(matrixStats)
library(reshape2)

#Script for calculating jaccard index between our signatures and published signatures

lit_dirname= "N:/immunsig_scripts_manuscript/gene_lists_papers/"
dirname="N:./immunsig_scripts_manuscript/immune_signature_discovery_workflow/"

db_genes_0= read.table(paste0(dirname,"dbscan_cluster_genes_param_0_GO_annot.tsv"),sep = "\t",stringsAsFactors = F)
immune_annot= c("Mast",
            # "Antigen presentation",
            "Macrophage",
            "Plasma",
            "NK",
            "NK/TCD8",
            "B cell",
            "Monocyte",
            # "Endothelial",
            "T CD4/B (IL12)",
            "DC immature",
            # "CAF endothelial",
            # "Cell cycle 1",
            # "Myocyte",
            # "CAF (integrin/TGFb)",
            "T CD4 memory/naive",
            "DC mature",
            # "T CD8 memory 1",
            "pDC",
            "T CD4 reg",
            # "Cell cycle 2",
            "Lymphoid cells"
            # "Cell cycle/dna repair",
            # "T CD8 memory 2",
            # "Plasma (identical protein binding)"
)

 db_genes_0=db_genes_0[db_genes_0$Go %in% immune_annot,]
# table(db_genes_0$Go)
unq_annot= unique(db_genes_0$annot)

lit_gene_list=read.table( paste0(lit_dirname,"all_immune_cell_signatures_curated.txt"),sep = "\t",header = T)

# cell types and source are mixed into one name 
cell_source= paste0(lit_gene_list$Cell_type,"_",gsub("(.*)_(.*)",x = lit_gene_list$source,replacement = "\\1"))
lit_gene_list$cell_source= cell_source
lit_annot= unique(lit_gene_list$cell_source)
lit_annot=sort(lit_annot)

########### Jaccard score calculation ###########################

lit_cell= foreach(i=seq_along(lit_annot),.combine = rbind) %do% {
  cell_f= lit_gene_list %>%
    filter(cell_source ==lit_annot[i])
  genes=as.character(cell_f$Gene)
  
  
  cell_col= foreach(j=seq_along(unq_annot),.combine = cbind) %do% {
    db_genes= rownames(db_genes_0)[db_genes_0$annot==unq_annot[j] ]
    round(jaccardSets(db_genes,genes),digits = 2) 
  }
  cell_col
}
rownames(lit_cell)=lit_annot
colnames(lit_cell)=unq_annot

lit_cell=t(lit_cell)


lit_cell=data.frame(lit_cell,annot=unique(db_genes_0$Go))

#length(Go_annot)
# Save the heatmap
# Heatmap is alphabetically ordered in rows and columns!

# pdf(paste0(dirname,"SuppFig2_literature_gene_sign_comparison.pdf"),width = 32,height = 8)
plot1=Heatmap(lit_cell[order(lit_cell$annot),-ncol(lit_cell)],column_title_side = "bottom",
             row_title = "cluster names",row_title_side = "left",show_row_names = T,cluster_rows = F,cluster_columns = F,
             column_names_rot = 45,name = "jaccard_index",row_names_side = "left", 
             heatmap_legend_param = list(legend_direction = "horizontal",
                                         legend_width = unit(3, "cm"), title_position = "lefttop"))+ 
  rowAnnotation(annotation= anno_text(which = "column",sort(lit_cell$annot),just="left"))
draw(plot1, heatmap_legend_side = "top")

# dev.off()

# Max jaccard scores for each signature :
score_df=data.frame(signature=immune_annot,
           maxs=rowMaxs(as.matrix(lit_cell[,-ncol(lit_cell)]))) %>% arrange(desc(maxs))
##signatures above 0.1
score_df[score_df$maxs>0.1,]

# 1       B cell 0.31
# 2      NK/TCD8 0.17
# 3     Monocyte 0.17
# 4  DC immature 0.17
# 5 Cell cycle 1 0.15
# 6           NK 0.13
# 7         Mast 0.11

#7/14 ~ 50% above 0.1 Jaccard score 
nrow(score_df)
median(score_df$maxs) #~0.09

max_jaccard_coeff=foreach(i=seq_along(rownames(lit_cell)),.combine = rbind) %do% {
  sort(lit_cell[i,],decreasing = T) [1:4] %>% melt()
  
}
max_jaccard_coeff %>% arrange(desc(value)) %>%head(10)

# annot        variable value
# 1      B cell        B_Nirmal  0.31
# 2      B cell  B_naive_Newman  0.24
# 3      B cell         B_Becht  0.20
# 4     NK/TCD8    T_CD8_Newman  0.17
# 5    Monocyte Monocyte_Newman  0.17
# 6 DC immature    DC_mDC_Becht  0.17

################

# overlap coefficient

overlap_coeff=function(x,y){
  length(intersect(x,y))/min(length(x),length(y))
}



overlap_lit_cell= foreach(i=seq_along(lit_annot),.combine = rbind) %do% {
  cell_f= lit_gene_list %>%
    filter(cell_source ==lit_annot[i])
  genes=as.character(cell_f$Gene)
  
  
  cell_col= foreach(j=seq_along(unq_annot),.combine = cbind) %do% {
    db_genes= rownames(db_genes_0)[db_genes_0$annot==unq_annot[j] ]
    round(overlap_coeff(db_genes,genes),digits = 2) 
  }
  cell_col
}
rownames(overlap_lit_cell)=lit_annot
colnames(overlap_lit_cell)=unq_annot

overlap_lit_cell=t(overlap_lit_cell)


overlap_lit_cell=data.frame(overlap_lit_cell,annot=unique(db_genes_0$Go))

#length(Go_annot)
# Save the heatmap
# Heatmap is alphabetically ordered in rows and columns!

# pdf(paste0(dirname,"SuppFig2_literature_gene_sign_comparison_Overlap_coefficient.pdf"),width = 32,height = 8)
plot2=Heatmap(overlap_lit_cell[order(overlap_lit_cell$annot),-ncol(overlap_lit_cell)],column_title = "Literature signatures",column_title_side = "bottom",
             row_title = "cluster names",row_title_side = "left",show_row_names = T,cluster_rows = F,cluster_columns = F,
             column_names_rot = 45,name = "overlap coefficient",row_names_side = "left", 
             heatmap_legend_param = list(legend_direction = "horizontal",
                                         legend_width = unit(3, "cm"), title_position = "lefttop"))+ 
  rowAnnotation(annotation= anno_text(which = "column",sort(overlap_lit_cell$annot),just="left"))
draw(plot2, heatmap_legend_side = "top")

# dev.off()

# Max overlap scores for each signature :
score_df=data.frame(signature=immune_annot,
                    maxs=rowMaxs(as.matrix(overlap_lit_cell[,-ncol(overlap_lit_cell)]))) %>% arrange(desc(maxs))
##signatures above 0.1
score_df[score_df$maxs>0.1,]
# 1              B cell 0.89
# 2             NK/TCD8 0.80
# 3      Lymphoid cells 0.80
# 4          Macrophage 0.75
# 5  T CD4 memory/naive 0.57
# 6                  NK 0.56
# 7      T CD4/B (IL12) 0.43
# 8         DC immature 0.38
# 9              Plasma 0.30
# 10           Monocyte 0.30
# 11          T CD4 reg 0.29
# 12               Mast 0.27


#12/14 ~ above 0.2 
# only 2 below 0.2
# 13                pDC 0.07
# 14          DC mature 0.06

median(score_df$maxs) #~0.4


max_overlap_coeff=foreach(i=seq_along(rownames(overlap_lit_cell)),.combine = rbind) %do% {
  sort(overlap_lit_cell[i,],decreasing = T) [1:4] %>% melt()
  
}
max_overlap_coeff %>% arrange(desc(value)) %>%head()
# annot              variable value
# 1         B cell               B_Becht  0.89
# 2        NK/TCD8          T_CD8_Newman  0.80
# 3 Lymphoid cells              T_Nirmal  0.80
# 4     Macrophage Macrophages_DC_Nirmal  0.75
# 5        NK/TCD8   NK_activated_Newman  0.67
# 6        NK/TCD8     NK_resting_Newman  0.67

####
library(gridExtra)
svg(paste0(dirname,"SuppFig2_overlap_and_jaccard.svg"),width = 32,height = 16)

grob1 = grid.grabExpr(draw(plot1, heatmap_legend_side = "top")) 
grob2 = grid.grabExpr(draw(plot2, heatmap_legend_side = "top")) 

#grid.arrange(grob1,grob2) 
plot_grid(grob1, grob2, labels = c('A', 'B'),ncol = 1,label_size = 18)

dev.off()
