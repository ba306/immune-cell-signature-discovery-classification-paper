# Benchmarking results
# Create plots from benchmarking results
# Benchmarking datasets: Kotliarov, Zheng
# Methods
#   RF with our gene sets and other immune cell type signatures
#   Seurat
#   singleR
#   CHETAH
#   CellTypist
#   scType

library(Seurat)
library(qs)
library(ggplot2)
library("ggsci")
library(reshape2)
library(cowplot)
source("./functions/benchmarking_prediction_metrics.R") # to calculate prediction metrics


dir_analysis="./benchmarking/"

# Get cell type annotations from Zheng 
zheng=qread("./data/Zheng/pbmc_zheng_sorted_2k.qs")
zheng_annot=zheng$harmonized_celltype
zheng_annot=as.data.frame(zheng_annot)
head(zheng_annot)
rm(zheng)

# Get cell type annotations from Kotliarov
kotliarov=qread("./data/Kotliarov/kotliarov_pbmc.qs")

# table(kotliarov$harmonized_celltype)
# table(kotliarov$K3,kotliarov$harmonized_celltype)

kotliarov_annot=kotliarov$harmonized_celltype
kotliarov_annot=as.data.frame(kotliarov_annot)
head(kotliarov_annot)
rm(kotliarov)

# read every prediction result and compare with the original annotations
all_pred_files=list.files("./benchmarking/prediction_results/",recursive = T,full.names = T)

all_df=data.frame(
method=gsub("./benchmarking/prediction_results//.*/(.*)_(.*)_(.*)_(.*)_(.*).csv","\\1",all_pred_files),
ref=gsub("./benchmarking/prediction_results//.*/(.*)_(.*)_(.*)_(.*)_(.*).csv","\\2",all_pred_files),
query=gsub("./benchmarking/prediction_results//.*/(.*)_(.*)_(.*)_(.*)_(.*).csv","\\3",all_pred_files),
geneset=gsub("./benchmarking/prediction_results//.*/(.*)_(.*)_(.*)_(.*)_(.*).csv","\\4",all_pred_files),
no_genes=gsub("./benchmarking/prediction_results//.*/(.*)_(.*)_(.*)_(.*)_(.*).csv","\\5",all_pred_files),
file=all_pred_files
)


# Cell type harmonization

all_stats=foreach(i=seq_along(rownames(all_df)),.combine=rbind)%do%{
  print(all_df[i,"file"])
pred_res=read.csv(all_df[i,"file"],row.names = 1)

# Take the results from only Zheng or Kotliarov
if(all_df[i,"query"]=="Zheng"){
  pred_res=pred_res[rownames(zheng_annot),,drop=F]
}else if(all_df[i,"query"]=="Kotliarov"){
  pred_res=pred_res[rownames(kotliarov_annot),,drop=F]
  }else{NA}

# For scType harmonization of their database 
if(all_df[i,"method"]=="scType"){
  pred=sctype_harmonization(pred_res$Predictions)
  
}else{
  pred=pred_res$Predictions
}

# Benchmarking- calculating prediction scores

if(all_df[i,"query"]=="Zheng"){
  harmonize_pred_stat_mean(factor(pred),factor(zheng_annot$zheng_annot))
}else if(all_df[i,"query"]=="Kotliarov"){
  harmonize_pred_stat_mean(factor(pred),factor(kotliarov_annot$kotliarov_annot))
  
}else{NA}
}

#add predictions to the info about the method, parameters ...
all_info_stats=cbind(all_df,all_stats)

# For plotting purpose improve the method annotations
all_info_stats$method=case_when(all_info_stats$method=="rf"~"RF",
                                all_info_stats$method=="celltypist"~"CellTypist",
                                T~          all_info_stats$method)

all_info_stats$geneset=case_when(all_info_stats$geneset=="angelova"~"Angelova",
                                all_info_stats$geneset=="abbas"~"Abbas",
                                all_info_stats$geneset=="charoentong"~"Charoentong",
                                all_info_stats$geneset=="nieto"~"Nieto",
                                all_info_stats$geneset=="random163genes"~"Random genes",
                                all_info_stats$geneset=="random167genes"~"Random genes",
                                
                                all_info_stats$geneset=="ourgenes"~"Our geneset",
                                
                                T~          all_info_stats$geneset)

# Represent the data in melted format
all_info_stats_melt=melt(all_info_stats)

# factoring/ordering of the methods and genesets for plotting 
methods=c(  "RF",
            "Seurat", 
            "CellTypist",
            "scType",
            "singleR" ,
            "CHETAH")

published_signatures=c("Our geneset",
                       "Charoentong",
                       "Angelova" ,
                       "Nieto",
                       "Abbas",
                       "Random genes")

# Benchmarking plots based on default settings from different classification methods
# Take only those we want to plot
take=   (all_info_stats$method=="singleR" & all_info_stats$geneset=="All") |
  (all_info_stats$method=="CellTypist" & all_info_stats$geneset=="default") |
  (all_info_stats$method=="CHETAH" & all_info_stats$geneset=="default") |
  (all_info_stats$method=="Seurat" & all_info_stats$no_genes=="2000") |
  (all_info_stats$method=="RF" & all_info_stats$geneset=="Our geneset") |
  (all_info_stats$method=="scType" & all_info_stats$no_genes=="2000")

all_methods=all_info_stats[take,]
all_methods=na.omit(all_methods)

all_methods=melt(all_methods)
# Order the methods
all_methods$method <- factor(all_methods$method, levels = methods)

# Plot comparing prediction scores from different methods in kotliarov and zheng datasets
all_methods_figure=all_methods%>%
ggplot( aes(variable, value,fill=method))+
  geom_bar(position=position_dodge(0.9), colour="black", stat="identity", width=0.9) +
    labs(x="Metrics", y = "Metrics score")+ theme_bw(base_size=22)+coord_flip()+
  guides(shape = guide_legend(override.aes = list(size = 0.5)))+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))+facet_wrap("query",scales = "free_x")+ 
  theme(legend.position="top")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  scale_fill_nejm()

# Plot comparing prediction scores from scType, CellTypist and Seurat in kotliarov and zheng datasets
# Different HVGs are employed in these three top ranking methods
# Our RF approach was shown as red dotted line and point

#take only celltypist, sctype, seurat and our RF
take=   (all_info_stats$method=="CellTypist" |
           all_info_stats$method=="Seurat"|
  all_info_stats$method=="scType" |
    all_info_stats$geneset=="Our geneset") &
  all_info_stats$geneset!="ImmuneAllLow" &
  all_info_stats$geneset!="default"
  
hvgs_methods=all_info_stats[take,]
hvgs_methods=na.omit(hvgs_methods)

hvgs_methods=melt(hvgs_methods)
hvgs_methods$no_genes=as.numeric(hvgs_methods$no_genes)

hvgs_methods$hline=ifelse(
  hvgs_methods$method=="RF" ,hvgs_methods$value,NA
)

hvgs_methods$method <- factor(hvgs_methods$method, levels = methods)


hvgs_plot=ggplot(hvgs_methods,aes(no_genes, value,colour=method)) + geom_point() + geom_line()+ 
  scale_x_continuous(breaks = seq(0, 2000, by = 500))+
  geom_hline(aes(yintercept = hline), linetype = "dashed", color = "red") +
  facet_grid(query~variable,scales = "free")+ 
  xlab("Number of HVGs") + ylab("Metrics score")+ labs(colour = "Metrics")+theme_bw(base_size=22)+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))+ theme(legend.position="top",
                                                      legend.margin=margin(0,0,0,0),
                                                      legend.box.margin=margin(-10,-10,-10,-10))+scale_color_nejm()


# Plots comparing different immune cell type signatures on RF approach
take= all_info_stats$method=="RF"

rf_signatures_pred=all_info_stats[take,]
rf_signatures_pred=na.omit(rf_signatures_pred)
rf_signatures_pred=melt(rf_signatures_pred)

rf_signatures_pred$geneset <- factor(rf_signatures_pred$geneset, levels = published_signatures)

rf_signatures_plot=rf_signatures_pred%>%
  ggplot( aes(variable, value,fill=geneset ))+
  geom_bar(position=position_dodge(0.9), colour="black", stat="identity", width=0.9) +
  labs(x="Metrics", y = "Metrics score")+ theme_bw(base_size=22)+coord_flip()+
  guides(shape = guide_legend(override.aes = list(size = 0.5)))+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))+facet_wrap("query",scales = "free_x")+ 
  theme(legend.position="top")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  scale_fill_nejm()


# Main figure plots from benchmarking part

png(paste0(dir_analysis,"Benchmarking_methods_plots.png"),width = 1400,height = 1200)

plot_grid(all_methods_figure, hvgs_plot, labels = c("A", "B"), nrow = 2, align = "h",label_size = 20,vjust = 1.2,
          hjust = -0.2)
dev.off()

png(paste0(dir_analysis,"Benchmarking_RF_plots.png"),width = 1000,height = 600)

rf_signatures_plot
dev.off()



