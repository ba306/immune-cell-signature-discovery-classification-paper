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


# Calculate score differences between our  gene sets and other gene sets
# Get unique queries
unique_queries <- unique(rf_signatures_pred$query)

# Create an empty data frame to store the results
diff_df <- data.frame(
  query = character(),
  geneset = character(),
  sensitivity_diff = numeric(),
  ppv_diff = numeric(),
  npv_diff = numeric(),
  specificity_diff = numeric(),
  accuracy_diff = numeric(),
  f1_score_diff = numeric(),
  stringsAsFactors = FALSE
)

# Calculate differences for each query separately
for (query in unique_queries) {
  # Filter rows for the current query
  query_rows <- rf_signatures_pred[rf_signatures_pred$query == query, ]
  
  # Calculate differences for "Our geneset" against other genesets
  our_geneset_row <- query_rows[query_rows$geneset == "Our geneset", ]
  other_geneset_rows <- query_rows[query_rows$geneset != "Our geneset", ]
  
  sensitivity_diff <- our_geneset_row$Sensitivity - other_geneset_rows$Sensitivity
  ppv_diff <- our_geneset_row$PPV - other_geneset_rows$PPV
  npv_diff <- our_geneset_row$NPV - other_geneset_rows$NPV
  specificity_diff <- our_geneset_row$Specificity - other_geneset_rows$Specificity
  accuracy_diff <- our_geneset_row$Accuracy - other_geneset_rows$Accuracy
  f1_score_diff <- our_geneset_row$f1_score - other_geneset_rows$f1_score
  
  # Create a data frame with the differences for the current query
  query_diff_df <- data.frame(
    query = rep(query, length(sensitivity_diff)),
    geneset = other_geneset_rows$geneset,
    sensitivity_diff = sensitivity_diff,
    ppv_diff = ppv_diff,
    npv_diff = npv_diff,
    specificity_diff = specificity_diff,
    accuracy_diff = accuracy_diff,
    f1_score_diff = f1_score_diff,
    stringsAsFactors = FALSE
  )
  
  # Append the data frame to the result
  diff_df <- rbind(diff_df, query_diff_df)
}

# View the resulting differences
diff_df

# query      geneset            sensitivity_diff ppv_diff npv_diff specificity_diff accuracy_diff 
# 1  Kotliarov        Abbas             8.23     4.86     0.99             1.44          1.83
# 2  Kotliarov     Angelova             5.87     5.10     1.45             1.66          2.24
# 3  Kotliarov  Charoentong             2.44     2.40     0.58             0.79          1.07
# 4  Kotliarov        Nieto             4.32     3.14     0.92             1.04          1.37
# 5  Kotliarov Random genes            39.99    33.70     6.39             6.45         10.09
# 6      Zheng        Abbas            16.28     9.56     4.52             3.66          5.72
# 7      Zheng     Angelova            14.55    17.83     5.40             5.61          7.41
# 8      Zheng  Charoentong             8.01     5.26     2.07             2.96          3.57
# 9      Zheng        Nieto            10.48     9.70     3.57             3.88          4.98
# 10     Zheng Random genes            56.21    55.43    15.65            15.08         21.49

# f1_score_diff
# 1           6.54
# 2           5.48
# 3           2.42
# 4           3.72
# 5          36.97
# 6          13.16
# 7          16.17
# 8           6.70
# 9          10.10
# 10         55.87

diff_df %>%
  filter(geneset=="Charoentong")

# Create plots 
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



