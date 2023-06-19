# Calculate jaccard scores and Szymkiewicz–Simpson coefficient (overlap coeffient)
# between our gene sets and published immune cell type signatures

library(foreach)
library(ComplexHeatmap)
library(bayesbio)
library(readxl)
library(dplyr)
library(matrixStats)
library(reshape2)
library(gridExtra)
library(cowplot)
#Script for calculating jaccard index between our signatures and published signatures

lit_dirname= "./gene_lists_papers/"
dirname="./discovery/"

db_genes_0= read.table(paste0(dirname,
                              "dbscan_cluster_genes_param_0_annot.tsv"),
                       sep = "\t",stringsAsFactors = F,header=T)
rownames(db_genes_0)=db_genes_0$genes

# Remove lineage signatures
# Remove invalidated S_6 monocytes signauture
immune_annot= 
  c("S_3",
    "S_7",
    "S_6")
db_genes_0=db_genes_0[!db_genes_0$ID %in% immune_annot,]

db_genes_0$Annotation=paste0(
  db_genes_0$Annotation," (",db_genes_0$ID,")"
)

unq_annot= unique(db_genes_0$Annotation)

# Load published immune cell signatures
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
    db_genes= rownames(db_genes_0)[db_genes_0$Annotation==unq_annot[j] ]
    round(jaccardSets(db_genes,genes),digits = 2) 
  }
  cell_col
}
rownames(lit_cell)=lit_annot
colnames(lit_cell)=unq_annot

lit_cell=t(lit_cell)

lit_cell=data.frame(lit_cell,annot=unique(db_genes_0$Annotation))

# Heatmap is alphabetically ordered in rows and columns!
# no of genes in brackets
lit_cell$annot=paste0(lit_cell$annot," (",table(db_genes_0$Annotation)[rownames(lit_cell)],")")
colnames(lit_cell)[-ncol(lit_cell)]=paste0(colnames(lit_cell)[-ncol(lit_cell)],
                                           " (",table(lit_gene_list$cell_source)[colnames(lit_cell)[-ncol(lit_cell)]],")")

plot1=Heatmap(lit_cell[order(lit_cell$annot),-ncol(lit_cell)],column_title_side = "bottom",
              row_title = "cluster names",row_title_side = "left",show_row_names = T,cluster_rows = F,cluster_columns = F,
              column_names_rot = 45,name = "jaccard_index",row_names_side = "left", 
              heatmap_legend_param = list(legend_direction = "horizontal",
                                          legend_width = unit(3, "cm"), title_position = "lefttop"))+ 
  rowAnnotation(annotation= anno_text(which = "column",sort(lit_cell$annot),just="left"))

# Max jaccard scores for each signature :
score_df=data.frame(signature=rownames(lit_cell),
                    maxs=rowMaxs(as.matrix(lit_cell[,-ncol(lit_cell)]))) %>% arrange(desc(maxs))
summary(score_df$maxs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0200  0.0950  0.1200  0.1536  0.1900  0.3200
##signatures above 0.1
score_df[score_df$maxs>0.1,]

# signature maxs
# 1        Mast (S_1) 0.32
# 2          B (S_13) 0.31
# 3      T CD4 (S_11) 0.19
# 4         DC (S_14) 0.19
# 5          NK (S_4) 0.16
# 6 Macrophage (S_10) 0.12
# 7       T CD8 (S_2) 0.12

# => 7/11

# What are top matches?

max_jaccard_coeff=foreach(i=seq_along(rownames(lit_cell)),.combine = rbind) %do% {
  sort(lit_cell[i,],decreasing = T) [1:4] %>% melt()
  
}
max_jaccard_coeff %>% arrange(desc(value)) %>%head(10)

# annot                   variable value
# 1    Mast (S_1) (20)           Mast_Bindea (30)  0.32
# 2      B (S_13) (38)              B_Nirmal (42)  0.31
# 3      B (S_13) (38)        B_naive_Newman (60)  0.21
# 4  T CD4 (S_11) (16)               T_Becht (16)  0.19
# 5      B (S_13) (38)       B_memory_Newman (62)  0.19
# 6     DC (S_14) (13)           DC_mDC_Becht (6)  0.19
# 7      NK (S_4) (20)               NK_Becht (9)  0.16
# 8    Mast (S_1) (20) Mast_activated_Newman (39)  0.13
# 9      NK (S_4) (20)   NK_activated_Newman (59)  0.13
# 10     NK (S_4) (20)     NK_resting_Newman (57)  0.13


########### Szymkiewicz–Simpson coefficient calculation ###########################

# Function to calculate Szymkiewicz–Simpson coefficient
overlap_coeff=function(x,y){
  length(intersect(x,y))/min(length(x),length(y))
}

overlap_lit_cell= foreach(i=seq_along(lit_annot),.combine = rbind) %do% {
  cell_f= lit_gene_list %>%
    filter(cell_source ==lit_annot[i])
  genes=as.character(cell_f$Gene)
  
  
  cell_col= foreach(j=seq_along(unq_annot),.combine = cbind) %do% {
    db_genes= rownames(db_genes_0)[db_genes_0$Annotation==unq_annot[j] ]
    round(overlap_coeff(db_genes,genes),digits = 2) 
  }
  cell_col
}
rownames(overlap_lit_cell)=lit_annot
colnames(overlap_lit_cell)=unq_annot

overlap_lit_cell=t(overlap_lit_cell)


overlap_lit_cell=data.frame(overlap_lit_cell,annot=unique(db_genes_0$Annotation))

# no of genes in brackets
overlap_lit_cell$annot=paste0(overlap_lit_cell$annot," (",table(db_genes_0$Annotation)[rownames(overlap_lit_cell)],")")
colnames(overlap_lit_cell)[-ncol(overlap_lit_cell)]=paste0(colnames(overlap_lit_cell)[-ncol(overlap_lit_cell)],
                                                           " (",table(lit_gene_list$cell_source)[colnames(overlap_lit_cell)[-ncol(overlap_lit_cell)]],")")
plot2=Heatmap(overlap_lit_cell[order(overlap_lit_cell$annot),-ncol(overlap_lit_cell)],column_title = "Literature signatures",column_title_side = "bottom",
              row_title = "cluster names",row_title_side = "left",show_row_names = T,cluster_rows = F,cluster_columns = F,
              column_names_rot = 45,name = "overlap coefficient",row_names_side = "left", 
              heatmap_legend_param = list(legend_direction = "horizontal",
                                          legend_width = unit(3, "cm"), title_position = "lefttop"))+ 
  rowAnnotation(annotation= anno_text(which = "column",sort(overlap_lit_cell$annot),just="left"))
draw(plot2, heatmap_legend_side = "top")


# Max overlap scores for each signature :
score_df=data.frame(signature=rownames(overlap_lit_cell),
                    maxs=rowMaxs(as.matrix(overlap_lit_cell[,-ncol(overlap_lit_cell)]))) %>% arrange(desc(maxs))
summary(score_df$maxs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0800  0.2500  0.4000  0.3709  0.4750  0.6700 
##signatures above 0.1
score_df[score_df$maxs>0.1,]

# signature maxs
# 1           B (S_13) 0.67
# 2         Mast (S_1) 0.60
# 3          DC (S_14) 0.50
# 4           NK (S_4) 0.45
# 5        T CD8 (S_2) 0.43
# 6  Macrophage (S_10) 0.40
# 7       T CD4 (S_11) 0.31
# 8    Monocytes (S_8) 0.31
# 9   Macrophage (S_5) 0.19
# 10     Plasma (S_12) 0.14



# What are top matches?

max_overlap_coeff=foreach(i=seq_along(rownames(overlap_lit_cell)),.combine = rbind) %do% {
  sort(overlap_lit_cell[i,],decreasing = T) [1:4] %>% melt()
  
}
max_overlap_coeff %>% arrange(desc(value)) %>%head()
# annot              variable value
# annot                 variable value
# 1   B (S_13) (38)              B_Becht (9)  0.67
# 2 Mast (S_1) (20)         Mast_Bindea (30)  0.60
# 3   B (S_13) (38)            B_Nirmal (42)  0.50
# 4  DC (S_14) (13)         DC_mDC_Becht (6)  0.50
# 5   B (S_13) (38)      B_naive_Newman (60)  0.45
# 6   NK (S_4) (20) NK_activated_Newman (59)  0.45


# Draw both heatmaps together
svg(paste0(dirname,"overlap_and_jaccard.svg"),width = 46,height = 16)

grob1 = grid.grabExpr(draw(plot1, heatmap_legend_side = "top")) 
grob2 = grid.grabExpr(draw(plot2, heatmap_legend_side = "top")) 

plot_grid(grob1, grob2, labels = c('A', 'B'),ncol = 1,label_size = 18)

dev.off()

