# Show downstream analysis bias based on IFN-g effect on DC 
# Dataset: IFNg PBMC Kartha scRNA-seq
# Cell labels from celltypist, seurat, sctype and our RF
# For differential gene expression only 2k HVGs are used 

library(cowplot)
library(ggplot2)
library(foreach)
library(qs)
library(Seurat)
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

# average number of genes per cell
mean(query@meta.data$nFeature_RNA)
#943.6993

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


##### Differential expression pvalues from 2k HVGs on DCs ##### 

# Extract 2k HVGs from integrated matrix
integ <- qread("./data/Kartha/Kartha_integrated.qs")
hvgs=rownames(integ)
rm(integ)


# Calculate p-values for each cell type labeling method
# Select only DCs
# IFNg vs control
# Separately for each time point

Idents(query)=query$Stimulation


# RF
diff_test=query[,query$rf_annot=="DC"]

cond=unique(query$time)

diff_ourcells_hvgs=foreach(time=cond,.combine = rbind)%do%{
  print(time)
  diff_test_time=diff_test[,diff_test$time==time]
  
    diff_mat=FindMarkers(diff_test_time,ident.1="IFN",ident.2 = "Control",
                         features = hvgs,
                         logfc.threshold = -Inf,
                         min.pct = -Inf,min.cells.feature = 1,min.cells.group = 1)
    diff_mat$celltype="DC"
    diff_mat$genes=rownames(diff_mat)
    diff_mat$method="RF"
    diff_mat$stim="IFN-g"
    diff_mat$time=time
    diff_mat
  }


# Seurat
diff_test=query[,query$seurat_annot=="DC"]


diff_seurat_hvgs=foreach(time=cond,.combine = rbind)%do%{
  print(time)
  diff_test_time=diff_test[,diff_test$time==time]
  
  diff_mat=FindMarkers(diff_test_time,ident.1="IFN",ident.2 = "Control",
                       features = hvgs,
                       logfc.threshold = -Inf,
                       min.pct = -Inf,min.cells.feature = 1,min.cells.group = 1)
    diff_mat$celltype="DC"
    diff_mat$genes=rownames(diff_mat)
    diff_mat$method="Seurat"
    diff_mat$stim="IFN-g"
    diff_mat$time=time
    
    diff_mat
  }

# CellTypist
diff_test=query[,query$celltypist_annot=="DC"]

diff_celltypist_hvgs=foreach(time=cond,.combine = rbind)%do%{
  print(time)
  diff_test_time=diff_test[,diff_test$time==time]
  
  diff_mat=FindMarkers(diff_test_time,ident.1="IFN",ident.2 = "Control",
                       features = hvgs,
                       logfc.threshold = -Inf,
                       min.pct = -Inf,min.cells.feature = 1,min.cells.group = 1)
    diff_mat$celltype="DC"
    diff_mat$genes=rownames(diff_mat)
    diff_mat$method="CellTypist"
    diff_mat$stim="IFN-g"
    diff_mat$time=time
    diff_mat
  }

# scType
diff_test=query[,query$sctype_annot=="DC"]

diff_sctype_hvgs=foreach(time=cond,.combine = rbind)%do%{
  print(time)
  diff_test_time=diff_test[,diff_test$time==time]
  
  diff_mat=FindMarkers(diff_test_time,ident.1="IFN",ident.2 = "Control",
                       features = hvgs,
                       logfc.threshold = -Inf,
                       min.pct = -Inf,min.cells.feature = 1,min.cells.group = 1)
    diff_mat$celltype="DC"
    diff_mat$genes=rownames(diff_mat)
    diff_mat$method="scType"
    diff_mat$stim="IFN-g"
    diff_mat$time=time
    
    diff_mat
  }

#save the results as list
saveRDS(list(diff_ourcells_hvgs,diff_seurat_hvgs,diff_celltypist_hvgs,diff_sctype_hvgs),
        file = paste0(dir_analysis,"2khvgs_diff_exp_DC_annots.RDS"))

# p-value-to-p-value plots 
# scatter plots

# create dataframe for plotting from the diff exp results 
#    Our_pval celltype     gene  stim time   other_pval

# Compare pvalues from RF vs other methods

# seurat vs ours 
all_df_seurat=
  foreach(t=cond,.combine = rbind)%do%{
    a=diff_ourcells_hvgs %>%
      filter(time==t) 
    
    b=diff_seurat_hvgs %>%
      filter(time==t) 
    
    cbind(
      data.frame(
        Our_pval=
          a
        [a$genes %in% hvgs,][,c("p_val")],
        gene=hvgs,
        celltype="DC",
        stim="IFN-g",
        time=t
      ),
      data.frame(
        other_pval=
          b
        [b$genes %in% hvgs,][,c("p_val")]
        
      )
    )
  }


#celltypist vs ours

all_df_celltypist=
  foreach(t=cond,.combine = rbind)%do%{
    a=diff_ourcells_hvgs %>%
      filter(time==t) 
    
    b=diff_celltypist_hvgs %>%
      filter(time==t) 
    
    cbind(
      data.frame(
        Our_pval=
          a
        [a$genes %in% hvgs,][,c("p_val")],
        gene=hvgs,
        celltype="DC",
        stim="IFN-g",
        time=t
      ),
      data.frame(
        other_pval=
          b
        [b$genes %in% hvgs,][,c("p_val")]
        
      )
    )
  }

# sctype vs ours
all_df_sctype=
  foreach(t=cond,.combine = rbind)%do%{
    a=diff_ourcells_hvgs %>%
      filter(time==t) 
    
    b=diff_sctype_hvgs %>%
      filter(time==t) 
    
    cbind(
      data.frame(
        Our_pval=
          a
        [a$genes %in% hvgs,][,c("p_val")],
        gene=hvgs,
        celltype="DC",
        stim="IFN-g",
        time=t
      ),
      data.frame(
        other_pval=
          b
        [b$genes %in% hvgs,][,c("p_val")]
        
      )
    )
  }


# p value from other method on y axis
# p value from RF on x axis
# trend line added

seurat_plot_dc=all_df_seurat%>% 
ggplot( aes(x = Our_pval, y = other_pval)) +
  geom_point(
    size = 1) +ylab("Seurat p values (p)")+xlab("RF p values (p)")+
  geom_abline(color = "black")+facet_grid(celltype~time)+
  theme_bw(base_size=15)+  theme( legend.position = "top")

celltypist_plot_dc=all_df_celltypist%>% 
  ggplot( aes(x = Our_pval, y = other_pval)) +
  geom_point(
    size = 1) +ylab("CellTypist p values (p)")+xlab("RF p values (p)")+
  geom_abline(color = "black")+facet_grid(celltype~time)+
  theme_bw(base_size=15)+  theme( legend.position = "top")

sctype_plot_dc=all_df_celltypist%>% 
  ggplot( aes(x = Our_pval, y = other_pval)) +
  geom_point(
    size = 1) +ylab("scType p values (p)")+xlab("RF p values (p)")+
  geom_abline(color = "black")+facet_grid(celltype~time)+
  theme_bw(base_size=15)+  theme( legend.position = "top")


png(paste0(dir_analysis,"Kartha_DC_pvalue_plots.png"),width = 1400,height = 400)

plot_grid(seurat_plot_dc, sctype_plot_dc, celltypist_plot_dc, labels = c("A", "B", "C"), nrow = 1, align = "h",label_size = 20,vjust = 1.2,
          hjust = -0.2)


dev.off()


