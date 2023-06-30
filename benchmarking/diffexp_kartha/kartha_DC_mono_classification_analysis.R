# Show downstream analysis bias based on IFN-g effect on DC 
# Dataset: IFNg PBMC Kartha scRNA-seq
# Cell labels from celltypist, seurat, sctype and our RF
# DC classifications checked in each cell typing method
# Misclassified DCs detected --> in reality monocytes with high IFNg response

library(qs)
library(Seurat)
library(msigdbr)
library(ggsankey)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(gtable)
library(ggbeeswarm) # to draw violin plots together with jitter
set.seed(42)

dir_analysis="./benchmarking/diffexp_kartha/"

##### Prepare and load data ##### 
# Load Kartha data
query="./data/Kartha/Kartha_sobj.qs"
query <- qs::qread(file = query)

# Only focus on IFN and control for the demonstration of downstream analysis bias

take=query$StimType=="Control" | query$StimType=="IFN"
query=query[,take]

# Remove ControlGolgiPlug_6h and IFNGolgiPlug_6h
# remove donor 4 since this only has IFNg treatment
remove=query$Condition=="ControlGolgiPlug_6h" | query$Condition=="IFNGolgiPlug_6h" | query$Donor =="Donor4"
query=query[,!remove]

# Load different cell type classification results
seurat_annot=read.csv("./benchmarking/prediction_results/seurat/Seurat_Hao_Kartha_HVGs_2000.csv",header = T,row.names = 1)
sctype_annot=read.csv("./benchmarking/prediction_results/scType/scType_Hao_Kartha_HVGs_2000.csv",header = T,row.names = 1)
celltypist_annot=read.csv("./benchmarking/prediction_results/celltypist/celltypist_Hao_Kartha_default_1060.csv",header = T,row.names = 1)
rf_annot=read.csv("./benchmarking/prediction_results/RF/rf_Hao_Kartha_ourgenes_164.csv",header = T,row.names = 1)

# harmonize scType annotations
source("./functions/benchmarking_prediction_metrics.R")
sctype_annot$Predictions =sctype_harmonization(sctype_annot$Predictions)


query$seurat_annot=seurat_annot[colnames(query),"Predictions"]
query$sctype_annot=sctype_annot[colnames(query),"Predictions"]
query$celltypist_annot=celltypist_annot[colnames(query),"Predictions"]
query$rf_annot=rf_annot[colnames(query),"Predictions"]

table(query$seurat_annot)
# B    DC  Mono    NK T CD4 T CD8 
# 447   126  1161   376  3264  1128 
table(query$sctype_annot)
# B                          DC ISG expressing immune cells 
# 473                         105                         279 
# Mono                       T CD4                    T unconv 
# 1240                        2856                        1549 
table(query$celltypist_annot)
# B    DC  Mono    NK T CD4 T CD8 
# 452   115  1364  2119  1262  1190 
table(query$rf_annot)
# B    DC  Mono    NK T CD4 T CD8 
# 411    66  1269   338  3124  1294 


##### Load DC, mono and IFNg Hallmark signatures ##### 
# DC markers 

DC_markers=c(
  "FCER1A","CD1C",
  "FLT3","CD1E")
monocyte_markers =c("CD14","FCGR3A",
                    #https://www.nature.com/articles/s41467-018-04985-0
                    "CTSS", "FCN1", "S100A9", "LYZ", "VCAN"   ,
                    #https://www.science.org/doi/10.1126/science.aah4573
                    "TLR2","ITGB2","ITGAM","CTSD","CTSA","NLRP3")

# IFNg Hallmark score 
hallmark_gene_sets <- msigdbr(species="human",category = "H") %>%
  filter(gs_name %in%"HALLMARK_INTERFERON_GAMMA_RESPONSE" )
IFNg_Hallmark=hallmark_gene_sets$human_gene_symbol

##### 
# Compare DCs annotated by RF (DC_RF) against 
# monocytes annotated by RF but annotated as DCs by other methods (Mono_RF)
# Do the comparison for each method
##### 

#####  Seurat vs RF ##### 
# Get DCs by RF or Seurat

take_dcs=query$seurat_annot=="DC" |
  query$rf_annot=="DC"
query_alldcs=query[,take_dcs]

comp_df=data.frame(RF=query_alldcs$rf_annot,
                   Seurat=query_alldcs$seurat_annot)

# Sankey plot
sankey_seurat=comp_df %>%
  make_long(RF, Seurat) %>%
  ggplot( aes(x = x, next_x = next_x,
              node = node, next_node = next_node,
              fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))+ggtitle("")

# Mono and DC markers in DC_RF vs Mono_RF
query_alldcs$celltype_comparison=case_when(query_alldcs$rf_annot=="DC" ~"DC_RF",
                                           query_alldcs$rf_annot=="Mono"&
                                             query_alldcs$seurat_annot=="DC" ~ "Mono_RF",
                                           T~ "Others")
query_alldcs=query_alldcs[,query_alldcs$celltype_comparison!="Others"]

# Dot plots- Mono and DC markers in DC_RF vs Mono_RF
mono_dot=DotPlot(query_alldcs, features = monocyte_markers,
                            group.by = "celltype_comparison", cols = "RdYlBu")+ 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.6, "cm"),# Adjust the size of the legend key
        legend.box = "vertical"  # Display the legend in a vertical layout
        )+coord_flip()+xlab("Monocyte markers")+ 
  ggtitle("Seurat vs RF")

dc_dot=DotPlot(query_alldcs, features = DC_markers,
                      group.by = "celltype_comparison", cols = "RdYlBu")+
  coord_flip()+
  xlab("DC markers")+
  theme(
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),# Adjust the size of the legend key
    legend.box = "vertical",  # Display the legend in a vertical layout
    legend.title = element_text(size = 10)
  )

seurat_dot=grid.arrange(mono_dot,dc_dot,nrow=2)

# IFNg scores in mono_RF and DC_RF

query_alldcs=ScaleData(query_alldcs,features = IFNg_Hallmark)
query_alldcs$IFNg_hallmark=colMeans(query_alldcs@assays$RNA@scale.data)

seurat_ifng=data.frame(celltype=query_alldcs$celltype_comparison,
                       sample=query_alldcs$Condition,
                       IFNg_hallmark=query_alldcs$IFNg_hallmark) %>%
  ggplot(aes(x = sample, y = IFNg_hallmark, fill = celltype)) +
  geom_violin(alpha = 0.5, scale = "width",position = position_dodge(width = 1)) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = 1, 
                               color="black",alpha=.5,show.legend = F)+
  stat_compare_means(aes(group = celltype), label = "p.signif",method = "wilcox.test")+
  theme_bw(base_size = 16)+
  ggtitle("Seurat vs RF") +
  theme(
    legend.position = "none"
  )+xlab("")


##### CellTypist vs RF ##### 

# Get DCs by RF or CellTypist
take_dcs=query$rf_annot=="DC" |
  query$celltypist_annot=="DC"
query_alldcs=query[,take_dcs]

comp_df=data.frame(RF=query_alldcs$rf_annot,
                   CellTypist=query_alldcs$celltypist_annot)

# Sankey plot
sankey_celltypist=comp_df %>%
  make_long(RF, CellTypist) %>%
  ggplot( aes(x = x, next_x = next_x,
              node = node, next_node = next_node,
              fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))+ggtitle("")

# Mono and DC markers in DC_RF vs Mono_RF
query_alldcs$celltype_comparison=case_when(query_alldcs$rf_annot=="DC" ~"DC_RF",
                                           query_alldcs$rf_annot=="Mono"&
                                             query_alldcs$celltypist_annot=="DC" ~ "Mono_RF",
                                           T~ "Others")
query_alldcs=query_alldcs[,query_alldcs$celltype_comparison!="Others"]

# Dot plots- Mono and DC markers in DC_RF vs Mono_RF
mono_dot=DotPlot(query_alldcs, features = monocyte_markers,
                 group.by = "celltype_comparison", cols = "RdYlBu")+ 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.6, "cm"),# Adjust the size of the legend key
        legend.box = "vertical"  # Display the legend in a vertical layout
  )+coord_flip()+xlab("Monocyte markers")+ 
  ggtitle("CellTypist vs RF")
dc_dot=DotPlot(query_alldcs, features = DC_markers,
               group.by = "celltype_comparison", cols = "RdYlBu")+
  coord_flip()+
  xlab("DC markers")+
  theme(        axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),# Adjust the size of the legend key
    legend.box = "vertical",  # Display the legend in a vertical layout
    legend.title = element_text(size = 10)
  )

celltypist_dot=grid.arrange(mono_dot,dc_dot,nrow=2)

# IFNg scores in mono_RF and DC_RF

query_alldcs=ScaleData(query_alldcs,features = IFNg_Hallmark)
query_alldcs$IFNg_hallmark=colMeans(query_alldcs@assays$RNA@scale.data)


celltypist_ifng=data.frame(celltype=query_alldcs$celltype_comparison,
                           sample=query_alldcs$Condition,
                           IFNg_hallmark=query_alldcs$IFNg_hallmark) %>%
  ggplot(aes(x = sample, y = IFNg_hallmark, fill = celltype)) +
  geom_violin(alpha = 0.5, scale = "width",position = position_dodge(width = 1)) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = 1, 
                               color="black",alpha=.5,show.legend = F)+
  stat_compare_means(aes(group = celltype), label = "p.signif",method = "wilcox.test")+
  theme_bw(base_size = 16)+
  ggtitle("CellTypist vs RF") +
  theme(
    legend.position = "top"
  )+ labs(fill = "RF myeloid celltypes")+
  ylab("")

# get the legend to plot in the middle only
legend_plot=gtable_filter(ggplot_gtable(ggplot_build(celltypist_ifng)), "guide-box")

celltypist_ifng=data.frame(celltype=query_alldcs$celltype_comparison,
                           sample=query_alldcs$Condition,
                           IFNg_hallmark=query_alldcs$IFNg_hallmark) %>%
  ggplot(aes(x = sample, y = IFNg_hallmark, fill = celltype)) +
  geom_violin(alpha = 0.5, scale = "width",position = position_dodge(width = 1)) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = 1, 
                               color="black",alpha=.5,show.legend = F)+
  stat_compare_means(aes(group = celltype), method = "wilcox.test",label = "p.signif")+
  theme_bw(base_size = 16)+
  ggtitle("CellTypist vs RF") +
  theme(
    legend.position = "none"
  )+ labs(fill = "RF myeloid celltypes")+
  ylab("")


##### sctype vs RF ##### 
# Get DCs by RF or sctype

take_dcs=query$rf_annot=="DC" |
  query$sctype_annot=="DC"
query_alldcs=query[,take_dcs]

comp_df=data.frame(RF=query_alldcs$rf_annot,
                   scType=query_alldcs$sctype_annot)

# Sankey plot
sankey_sctype=comp_df %>%
  make_long(RF, scType) %>%
  ggplot( aes(x = x, next_x = next_x,
              node = node, next_node = next_node,
              fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))+ggtitle("")

# Mono and DC markers in DC_RF vs Mono_RF
query_alldcs$celltype_comparison=case_when(query_alldcs$rf_annot=="DC" ~"DC_RF",
                                           query_alldcs$rf_annot=="Mono" &
                                             query_alldcs$sctype_annot=="DC"~ "Mono_RF",
                                           T~ "Others")
query_alldcs=query_alldcs[,query_alldcs$celltype_comparison!="Others"]

# Dot plots- Mono and DC markers in DC_RF vs Mono_RF
mono_dot=DotPlot(query_alldcs, features = monocyte_markers,
                 group.by = "celltype_comparison", cols = "RdYlBu")+ 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),        
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.6, "cm"),# Adjust the size of the legend key
        legend.box = "vertical"  # Display the legend in a vertical layout
  )+coord_flip()+xlab("Monocyte markers")+ 
  ggtitle("scType vs RF")
dc_dot=DotPlot(query_alldcs, features = DC_markers,
               group.by = "celltype_comparison", cols = "RdYlBu")+
  coord_flip()+
  xlab("DC markers")+
  theme(        
    axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),# Adjust the size of the legend key
    legend.box = "vertical",  # Display the legend in a vertical layout
    legend.title = element_text(size = 10)
  )

sctype_dot=grid.arrange(mono_dot,dc_dot,nrow=2)

# IFNg scores in mono_RF and DC_RF

query_alldcs=ScaleData(query_alldcs,features = IFNg_Hallmark)
query_alldcs$IFNg_hallmark=colMeans(query_alldcs@assays$RNA@scale.data)

sctype_ifng=data.frame(celltype=query_alldcs$celltype_comparison,
                       sample=query_alldcs$Condition,
                       IFNg_hallmark=query_alldcs$IFNg_hallmark) %>%
  # filter(celltype=="DC"|celltype=="Mono") %>%
  ggplot(aes(x = sample, y = IFNg_hallmark, fill = celltype)) +
  geom_violin(alpha = 0.5, scale = "width",position = position_dodge(width = 1)) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = 1, 
                               color="black",alpha=.5,show.legend = F)+
  stat_compare_means(aes(group = celltype), label = "p.signif",method = "wilcox.test")+
  theme_bw(base_size = 16)+
  ggtitle("scType vs RF") +
  theme(
    legend.position = "none"
  )+xlab("")

# Plots

# Dot plots
png(paste0(dir_analysis,"dotplots_mono_DC_allcond.png"),width = 1200,height = 500)
grid.arrange(seurat_dot,celltypist_dot,sctype_dot,nrow=1)
dev.off()

# sankey plot
png(paste0(dir_analysis,"sankey_DCs.png"),width = 800,height = 400)
grid.arrange(sankey_seurat,sankey_celltypist,sankey_sctype,nrow=1)
dev.off()

# IFN violin plots
png(paste0(dir_analysis,"IFNg_DCs_monos.png"),width = 1000,height = 400)
grid.arrange(legend_plot,arrangeGrob(
  seurat_ifng,celltypist_ifng,sctype_ifng,nrow=1),nrow=2,
  heights = c(0.1,2)
)
dev.off()


