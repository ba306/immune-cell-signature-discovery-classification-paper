# Curate cell type annotations for the workflow

library(dplyr)
library(varhandle)
dir="./data/TME/"
set.seed(42)

##################  Breast cancer ################
# GSE176078
GSE176078=read.table(paste0(dir,"BRCA_GSE176078_metadata.csv"),sep=",",
                         header=T,row.names = 1)
GSE176078=GSE176078%>%
  filter(!celltype_major %in%c("Normal Epithelial","CAFs",
                               "Cancer Epithelial","Endothelial","PVL"))%>%
  filter(!celltype_minor %in%c("Cycling T-cells","Cycling_Myeloid"))%>%
  unfactor()


table(GSE176078$celltype_minor)

GSE176078$major_celltypes=case_when(
  GSE176078$celltype_subset== "Myeloid_c4_DCs_pDC_IRF7"~ "pDC",
  GSE176078$celltype_minor== "B cells Naive"~ "B",
  GSE176078$celltype_minor== "B cells Memory"~ "B",
  GSE176078$celltype_minor== "DCs"~ "DC",
  GSE176078$celltype_minor== "T cells CD8+"~ "T CD8",
  GSE176078$celltype_minor== "Monocyte"~ "Monocytes",
  
  GSE176078$celltype_minor== "T cells CD4+"~ "T CD4",
  GSE176078$celltype_minor== "Plasmablasts"~ "Plasma",
  GSE176078$celltype_minor== "NKT cells"~ "NKT",
  GSE176078$celltype_minor== "NK cells"~ "NK",
  
  T~GSE176078$celltype_minor
  
)
table(GSE176078$major_celltypes)


table(GSE176078$celltype_subset)

GSE176078$minor_celltypes=case_when(
  GSE176078$celltype_subset== "T_cells_c2_CD4+_T-regs_FOXP3"~ "TCD4 reg",
  GSE176078$celltype_subset== "T_cells_c3_CD4+_Tfh_CXCL13"~ "TCD4 Th", # Tfh
  GSE176078$celltype_subset== "T_cells_c0_CD4+_CCR7"~ "TCD4 naive/memory", #naive/cm
  GSE176078$celltype_subset== "T_cells_c1_CD4+_IL7R"~ "TCD4 memory", #em
  
  GSE176078$celltype_subset== "Myeloid_c11_cDC2_CD1C"~ "cDC2",
  GSE176078$celltype_subset== "Myeloid_c3_cDC1_CLEC9A"~ "cDC1",
  GSE176078$celltype_subset== "Myeloid_c4_DCs_pDC_IRF7"~ "pDC",
  GSE176078$celltype_subset== "Myeloid_c4_DCs_pDC_IRF7"~ "pDC",
  GSE176078$celltype_subset== "B cells Memory"~ "B memory",
  GSE176078$celltype_subset== "B cells Nemory"~ "B naive",
  GSE176078$celltype_subset== "Myeloid_c7_Monocyte_3_FCGR3A"~ "Monocytes CD16",
  GSE176078$celltype_subset== "Myeloid_c8_Monocyte_2_S100A9"~ "Monocytes CD14",
  
  T~GSE176078$major_celltypes
  
)
table(GSE176078$minor_celltypes,GSE176078$major_celltypes)

table(GSE176078$minor_celltypes)

GSE176078$dataset="BRCA_GSE176078"

################## CRC ################
# GSE166555

GSE166555=read.table(paste0(dir,"CRC_GSE166555_metadata.tsv"),sep="\t",
                         header=T,row.names = 1)

head(GSE166555)
table(GSE166555$sample_origin)
GSE166555=GSE166555 %>%
  filter(sample_origin=="Tumor") %>%
  filter(main_cell_type=="Immune") %>%
  filter(!sct_cell_type %in%c("MT-hi","CD8+ IELs","Cycling TA",
                              "CD4+ Activated Fos-hi",
                              "Cycling Monocytes","Cycling B"))%>%
  unfactor()
table(GSE166555$sct_cell_type)

GSE166555$major_celltypes=case_when(
     GSE166555$sct_cell_type== "Follicular"~ "B",
     GSE166555$sct_cell_type== "NKs"~ "NK",
     GSE166555$sct_cell_type== "Tregs"~ "T CD4",
     GSE166555$sct_cell_type== "CD4+ Memory"~ "T CD4",
     GSE166555$sct_cell_type== "CD4+ PD1+"~ "T CD4",
     GSE166555$sct_cell_type== "Macrophages"~ "Macrophage",
     GSE166555$sct_cell_type== "CD69+ Mast"~ "Mast",
     GSE166555$sct_cell_type== "CD69- Mast"~ "Mast",
     GSE166555$sct_cell_type== "CD8+ IL17+"~ "T CD8",
     GSE166555$sct_cell_type== "CD8+ LP"~ "T CD8",
     GSE166555$sct_cell_type== "Inflammatory Monocytes"~ "Monocytes",
     GSE166555$sct_cell_type== "DC1"~ "DC",
     GSE166555$sct_cell_type== "DC2"~ "DC",
     GSE166555$sct_cell_type== "GC"~ "B",
     
       T~GSE166555$sct_cell_type
     
     )
table(GSE166555$major_celltypes)


GSE166555$minor_celltypes=case_when(
  GSE166555$sct_cell_type== "Follicular"~ "B follicular",
  GSE166555$sct_cell_type== "NKs"~ "NK",
  GSE166555$sct_cell_type== "Tregs"~ "TCD4 reg",
  GSE166555$sct_cell_type== "CD4+ Memory"~ "TCD4 memory",
  GSE166555$sct_cell_type== "CD4+ PD1+"~ "TCD4",
  GSE166555$sct_cell_type== "Macrophages"~ "Macrophage",
  GSE166555$sct_cell_type== "CD69+ Mast"~ "Mast",
  GSE166555$sct_cell_type== "CD69- Mast"~ "Mast",
  GSE166555$sct_cell_type== "CD8+ IL17+"~ "TCD8",
  GSE166555$sct_cell_type== "CD8+ LP"~ "TCD8",
  GSE166555$sct_cell_type== "Inflammatory Monocytes"~ "Monocytes",
  GSE166555$sct_cell_type== "DC1"~ "cDC1",
  GSE166555$sct_cell_type== "DC2"~ "cDC2",
  GSE166555$sct_cell_type== "GC"~ "B germinal",
  
  T~GSE166555$sct_cell_type
  
)
table(GSE166555$minor_celltypes)


GSE166555$dataset="CRC_GSE166555"

  table(GSE166555$sct_cell_type,GSE166555$cell_type_imm)

################## Kidney renal cancer ################
# GSE139555
#non-T cells
GSE139555=read.table(paste0(dir,"KIRC_GSE139555_nont_metadata.txt"),sep="\t",
                     header=T,row.names = 1)
head(GSE139555)
table(GSE139555$ident)
table(GSE139555$patient)

GSE139555=GSE139555 %>%
  filter(source=="Tumor") %>%
  filter(patient %in%c("Renal1","Renal2","Renal3"))%>%
  filter(!ident %in%c("ATP8","COX1","CYTB","RNR2","CD9")) %>%
  unfactor()

table(GSE139555$ident)
GSE139555$major_celltypes=GSE139555$ident
GSE139555$major_celltypes=gsub("LST1-COTL1","Monocytes",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("LST1-MS4A7","Monocytes",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("LYZ-FTH1","Monocytes",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("LYZ-MCL1","Monocytes",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("LYZ-TYROBP","Monocytes",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("CD68","Macrophage",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("CD3D","Macrophage",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("NKG7","NK",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("HLA-DQB1","DC",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("NAPSB","DC",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("^B-cells-(.*)","B",GSE139555$major_celltypes)
GSE139555$major_celltypes=gsub("MZB1","Plasma",GSE139555$major_celltypes,fixed = T)
GSE139555$major_celltypes=gsub("IRF7","Mast",GSE139555$major_celltypes,fixed = T)

GSE139555$minor_celltypes=GSE139555$ident
GSE139555$minor_celltypes=gsub("LST1-COTL1","Monocytes CD16",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("LST1-MS4A7","Monocytes CD16",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("LYZ-FTH1","Monocytes CD16-",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("LYZ-MCL1","Monocytes CD16-",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("LYZ-TYROBP","Monocytes CD16-",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("CD68","Macrophage",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("CD3D","Macrophage TCR",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("NKG7","NK",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("HLA-DQB1","DC",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("NAPSB","pDC",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("^B-cells-(.*)","B",GSE139555$minor_celltypes)
GSE139555$minor_celltypes=gsub("MZB1","Plasma",GSE139555$minor_celltypes,fixed = T)
GSE139555$minor_celltypes=gsub("IRF7","Mast",GSE139555$minor_celltypes,fixed = T)

# T cells
GSE139555_T=read.table(paste0(dir,"KIRC_GSE139555_tcell_metadata.txt"),sep="\t",
                     header=T,row.names = 1)
table(GSE139555_T$source)
table(GSE139555_T$patient)
table(GSE139555_T$ident)

GSE139555_T=GSE139555_T %>%
  filter(source=="Tumor") %>%
  filter(patient %in%c("Renal1","Renal2","Renal3"))%>%
  filter(!ident %in%c("3.1-MT","8.4-Chrom", "8.5-Mitosis")) %>%
  unfactor()

GSE139555_T$major_celltypes=case_when(
  grepl("^4(.*)",GSE139555_T$ident)~"T CD4",
  grepl("^8(.*)",GSE139555_T$ident)~"T CD8",
  T~GSE139555_T$ident
  
)
table(GSE139555_T$major_celltypes)
table(GSE139555_T$ident)

GSE139555_T$minor_celltypes=case_when(
  GSE139555_T$ident== "4.1-Trm"~ "TCD4 memory", #rm
  GSE139555_T$ident== "4.2-RPL32"~ "TCD4 Th", 
  GSE139555_T$ident== "4.3-TCF7"~ "TCD4 Th", 
  GSE139555_T$ident== "4.4-FOS"~ "TCD4 Th",
  GSE139555_T$ident== "4.5-IL6ST"~ "TCD4 Th",
  GSE139555_T$ident== "4.6a-Treg"~ "TCD4 reg",
  GSE139555_T$ident== "4.6b-Treg"~ "TCD4 reg",
  GSE139555_T$ident== "8.1-Teff"~ "TCD8 eff",
  GSE139555_T$ident== "8.2-Tem"~ "TCD8 memory", 
  GSE139555_T$ident== "8.3a-Trm"~ "TCD8 memory", 
  GSE139555_T$ident== "8.3b-Trm"~ "TCD8 memory", 
  GSE139555_T$ident== "8.3c-Trm"~ "TCD8 memory",
  GSE139555_T$ident== "8.6-KLRB1"~ "TCD8 cytotoxic",
  T~GSE139555_T$ident
  
)

table(GSE139555_T$minor_celltypes)

GSE139555=rbind(GSE139555,GSE139555_T)

GSE139555$dataset="KIRC_GSE139555"

#remove unclassified cells

GSE139555=GSE139555[!is.na(GSE139555$major_celltypes),]
rownames(GSE139555)=gsub("_","@",rownames(GSE139555))
################## Lung cancer ################

# NSCLC 
# GSE131907
GSE131907=read.table(paste0(dir,"NSCLC_GSE131907_metadata.txt"),sep="\t",
                     header=T,row.names = 1)

table(GSE131907$Cell_type)
table(GSE131907$Sample_Origin)
table(GSE131907$Sample)

GSE131907=GSE131907%>%
  filter(Sample_Origin=="tLung") %>%
  filter(!Cell_type %in%c("Endothelial cells","Epithelial cells",
           "Fibroblasts","Oligodendrocytes","Undetermined"))%>%
  filter(!Cell_subtype %in%c(  "CD8+/CD4+ Mixed Th",
                            "Undetermined"   )) %>%
  unfactor()


GSE131907$major_celltypes=case_when(
  GSE131907$Cell_subtype== "Activated DCs"~ "DC",
  GSE131907$Cell_subtype== "Alveolar Mac"~ "Macrophage",
  GSE131907$Cell_subtype== "CD141+ DCs"~ "DC",
  GSE131907$Cell_subtype== "CD163+CD14+ DCs"~ "DC",
  GSE131907$Cell_subtype== "CD1c+ DCs"~ "DC",
  GSE131907$Cell_subtype== "CD207+CD1a+ LCs"~ "DC",
  GSE131907$Cell_subtype== "CD4+ Th"~ "T CD4",
  GSE131907$Cell_subtype== "CD8 low T"~ "T CD8",
  GSE131907$Cell_subtype== "Cytotoxic CD8+ T"~ "T CD8",
  GSE131907$Cell_subtype== "Exhausted CD8+ T"~ "T CD8",
  GSE131907$Cell_subtype== "Exhausted Tfh"~ "T CD4",
  GSE131907$Cell_subtype== "Follicular B cells"~ "B",
  GSE131907$Cell_subtype== "GC B cells in the DZ"~ "B",
  GSE131907$Cell_subtype== "GC B cells in the LZ"~ "B",
  GSE131907$Cell_subtype== "MALT B cells"~ "B",
  GSE131907$Cell_subtype== "MAST"~ "Mast",
  GSE131907$Cell_subtype== "mo-Mac"~ "Macrophage",
  GSE131907$Cell_subtype== "Naive CD4+ T"~ "T CD4",
  GSE131907$Cell_subtype== "Naive CD8+ T"~ "T CD8",
  GSE131907$Cell_subtype== "pDCs"~ "pDC",
  GSE131907$Cell_subtype== "Plasma cells"~ "Plasma",
  GSE131907$Cell_subtype== "Pleural Mac"~ "Macrophage",
  GSE131907$Cell_subtype== "Treg"~ "T CD4",
  
  T~GSE131907$Cell_subtype
  
)
GSE131907$minor_celltypes=case_when(
  GSE131907$Cell_subtype== "Activated DCs"~ "DC activated",
  GSE131907$Cell_subtype== "Alveolar Mac"~ "Macrophage",
  GSE131907$Cell_subtype== "CD141+ DCs"~ "DC1",
  GSE131907$Cell_subtype== "CD163+CD14+ DCs"~ "DC",
  GSE131907$Cell_subtype== "CD1c+ DCs"~ "DC2",
  GSE131907$Cell_subtype== "CD207+CD1a+ LCs"~ "DC LC",
  GSE131907$Cell_subtype== "CD4+ Th"~ "TCD4 Th",
  GSE131907$Cell_subtype== "CD8 low T"~ "TCD8",
  GSE131907$Cell_subtype== "Cytotoxic CD8+ T"~ "TCD8 cytotoxic",
  GSE131907$Cell_subtype== "Exhausted CD8+ T"~ "TCD8 exhausted",
  GSE131907$Cell_subtype== "Exhausted Tfh"~ "TCD4 Tfh exhausted",
  GSE131907$Cell_subtype== "Follicular B cells"~ "B follicular",
  GSE131907$Cell_subtype== "GC B cells in the DZ"~ "B germinal",
  GSE131907$Cell_subtype== "GC B cells in the LZ"~ "B germinal",
  GSE131907$Cell_subtype== "MALT B cells"~ "B MALT",
  GSE131907$Cell_subtype== "MAST"~ "Mast",
  GSE131907$Cell_subtype== "mo-Mac"~ "Macrophage mo-Mac",
  GSE131907$Cell_subtype== "Naive CD4+ T"~ "TCD4 naive",
  GSE131907$Cell_subtype== "Naive CD8+ T"~ "TCD8 naive",
  GSE131907$Cell_subtype== "pDCs"~ "pDC",
  GSE131907$Cell_subtype== "Plasma cells"~ "Plasma",
  GSE131907$Cell_subtype== "Pleural Mac"~ "Macrophage",
  GSE131907$Cell_subtype== "Treg"~ "TCD4 reg",
  
  T~GSE131907$Cell_subtype
  
)

GSE131907$dataset="NSCLC_GSE131907"

#remove unclassified cells
table(GSE131907[(is.na(GSE131907$major_celltypes)),"Cell_type"])

GSE131907=GSE131907[!is.na(GSE131907$major_celltypes),]

table(GSE131907$major_celltypes)

################## Hepato ################
# GSE140228
GSE140228_drop=read.table(paste0(dir,"LIHC_GSE140228_Droplet_metadata.tsv"),sep="\t",
                          header=T,row.names = 1)
head(GSE140228_drop)

table(GSE140228_drop$Tissue)
GSE140228_drop=GSE140228_drop%>%
  filter(Tissue=="Tumor") %>%
  filter(!celltype_global %in%c("Lymphoid-T-NK-Cycling",
                                "Myeloid-Liver-doublets","ILCs"))%>% 
  filter(!celltype_sub %in%c(  "CD4/CD8-C1-CCR7","CD4/CD8-C2-MKI67" )) %>%
  unfactor()

sort(unique((GSE140228_drop$celltype_global)))
sort(unique((GSE140228_drop$celltype_sub)))

heatmap(table(GSE140228_drop$celltype_sub,GSE140228_drop$celltype_global))
table(GSE140228_drop$celltype_sub,GSE140228_drop$celltype_global)

GSE140228_drop$major_celltypes=case_when(
  # grepl("^CD4/CD8(.*)",GSE140228_drop$celltype_sub)~" T CD4/T CD8",
  grepl("^CD4(.*)",GSE140228_drop$celltype_sub)~"T CD4",
  grepl("^CD8(.*)",GSE140228_drop$celltype_sub)~"T CD8",
  
  grepl("^DC(.*)",GSE140228_drop$celltype_sub)~"DC",
  GSE140228_drop$celltype_sub== "Lymphoid-B"~ "B",
  GSE140228_drop$celltype_sub== "Lymphoid-B-Plasma"~ "Plasma",
  grepl("^Mac(.*)",GSE140228_drop$celltype_sub)~"Macrophage",
  grepl("^Mono(.*)",GSE140228_drop$celltype_sub)~"Monocytes",
  grepl("^NK(.*)",GSE140228_drop$celltype_sub)~"NK",
  grepl("^Mast(.*)",GSE140228_drop$celltype_sub)~"Mast",
  
  T~GSE140228_drop$celltype_sub
  
)
table(GSE140228_drop$celltype_sub,
      GSE140228_drop$major_celltypes)
table(GSE140228_drop$major_celltypes)
sort(unique((GSE140228_drop$celltype_sub)))

GSE140228_drop$minor_celltypes=case_when(
  GSE140228_drop$celltype_sub== "DC-C1-CD1C"~ "DC2",
  GSE140228_drop$celltype_sub== "DC-C2-FCER1A"~ "DC2",
  GSE140228_drop$celltype_sub== "DC-C3-CLEC9A"~ "DC1",
  GSE140228_drop$celltype_sub== "DC-C4-LAMP3"~ "DC",
  GSE140228_drop$celltype_sub== "Mono-C1-CD14"~ "Monocytes CD14",
  GSE140228_drop$celltype_sub== "Mono-C2-FCGR3A"~ "Monocytes CD16",
  GSE140228_drop$celltype_sub== "CD4-C3-ANXA1"~ "TCD4 memory", #cm
  GSE140228_drop$celltype_sub== "CD4-C4-IL7R"~ "TCD4 memory", #cm
  GSE140228_drop$celltype_sub== "CD4-C5-TCF7"~ "TCD4 naive",
  GSE140228_drop$celltype_sub== "CD4-C6-CXCL13"~ "TCD4 exhausted",
  GSE140228_drop$celltype_sub== "CD4-C7-FOXP3"~ "TCD4 reg",
  GSE140228_drop$celltype_sub== "CD8-C3-IL7R"~ "TCD8 memory", #memory
  GSE140228_drop$celltype_sub== "CD8-C4-CX3CR1"~ "TCD8 memory", #em
  GSE140228_drop$celltype_sub== "CD8-C5-SELL"~ "TCD8 memory", #cm
  GSE140228_drop$celltype_sub== "CD8-C6-GZMK"~ "TCD8 memory", #em 
  GSE140228_drop$celltype_sub== "CD8-C7-KLRD1"~ "TCD8 memory",#rm
  GSE140228_drop$celltype_sub== "CD8-C8-PDCD1"~ "TCD8 exhausted",
  GSE140228_drop$celltype_sub== "CD8-C9-SLC4A10"~ "TCD8 MAIT",
  T~GSE140228_drop$major_celltypes
  
)
table(GSE140228_drop$minor_celltypes)
GSE140228_drop$dataset="LIHC_GSE140228_10X"

GSE140228_smart=read.table(paste0(dir,"LIHC_GSE140228_Smartseq2_metadata.tsv"),sep="\t",
                           header=T,row.names = 1)
head(GSE140228_smart)

table(GSE140228_smart$Tissue)
sort(unique(GSE140228_smart$celltype_sub))

GSE140228_smart$Platform="Smartseq2"

GSE140228_smart=GSE140228_smart%>%
  filter(Tissue=="Tumor") %>%
  filter(!celltype_sub %in%c(  "ILCs" ,"CD8-C2-MKI67")) %>%
  
  unfactor()

sort(unique((GSE140228_smart$celltype_sub)))

GSE140228_smart$major_celltypes=case_when(
  grepl("^CD4(.*)",GSE140228_smart$celltype_sub)~"T CD4",
  grepl("^CD8(.*)",GSE140228_smart$celltype_sub)~"T CD8",
  
  grepl("^DC(.*)",GSE140228_smart$celltype_sub)~"DC",
  GSE140228_smart$celltype_sub== "Lymphoid-B"~ "B",
  GSE140228_smart$celltype_sub== "Lymphoid-B-Plasma"~ "Plasma",
  grepl("^Mac(.*)",GSE140228_smart$celltype_sub)~"Macrophage",
  grepl("^Mono(.*)",GSE140228_smart$celltype_sub)~"Monocytes",
  grepl("^NK(.*)",GSE140228_smart$celltype_sub)~"NK",
  grepl("^Mast(.*)",GSE140228_smart$celltype_sub)~"Mast",
  
  T~GSE140228_smart$celltype_sub
  
)
table(
  GSE140228_smart$major_celltypes)

sort(unique((GSE140228_smart$celltype_sub)))

GSE140228_smart$minor_celltypes=case_when(
  GSE140228_smart$celltype_sub== "DC-C1-CD1C"~ "DC2",
  GSE140228_smart$celltype_sub== "DC-C2-FCER1A"~ "DC2",
  GSE140228_smart$celltype_sub== "DC-C3-CLEC9A"~ "DC1",
  GSE140228_smart$celltype_sub== "DC-C4-LAMP3"~ "DC",
  GSE140228_smart$celltype_sub== "Mono-C1-CD14"~ "Monocytes CD14",
  GSE140228_smart$celltype_sub== "Mono-C2-FCGR3A"~ "Monocytes CD16",
  GSE140228_smart$celltype_sub== "CD4-C3-ANXA1"~ "TCD4 memory", #cm
  GSE140228_smart$celltype_sub== "CD4-C4-IL7R"~ "TCD4 memory", #cm
  GSE140228_smart$celltype_sub== "CD4-C5-TCF7"~ "TCD4 naive", 
  GSE140228_smart$celltype_sub== "CD4-C6-CXCL13"~ "TCD4 exhausted",
  GSE140228_smart$celltype_sub== "CD4-C7-FOXP3"~ "TCD4 reg",
  GSE140228_smart$celltype_sub== "CD8-C3-IL7R"~ "TCD8 memory",#cm
  GSE140228_smart$celltype_sub== "CD8-C4-CX3CR1"~ "TCD8 memory", #em
  GSE140228_smart$celltype_sub== "CD8-C5-SELL"~ "TCD8 memory", #cm
  GSE140228_smart$celltype_sub== "CD8-C6-GZMK"~ "TCD8 memory", #em
  GSE140228_smart$celltype_sub== "CD8-C7-KLRD1"~ "TCD8 memory", #rm
  GSE140228_smart$celltype_sub== "CD8-C8-PDCD1"~ "TCD8 exhausted",
  GSE140228_smart$celltype_sub== "CD8-C9-SLC4A10"~ "TCD8 MAIT",
  T~GSE140228_smart$major_celltypes
  
)
table(GSE140228_smart$minor_celltypes)
GSE140228_smart$dataset="LIHC_GSE140228_Smartseq2"

#combine both
inter_colmns=intersect(colnames(GSE140228_drop),
                       colnames(GSE140228_smart))
GSE140228=rbind(GSE140228_drop[,inter_colmns],GSE140228_smart[,inter_colmns])
GSE140228$dataset="LIHC_GSE140228"

################## Melanoma ################

# GSE123139

GSE123139=read.table(paste0(dir,"SKCM_GSE123139_metadata.tsv"),sep="\t",
                     header=T,row.names = 1)
head(GSE123139)

table(GSE123139$Batch.Set.ID)
table(GSE123139$treatment)
table(GSE123139$PatientID,GSE123139$treatment)
table(GSE123139$location)
table(GSE123139$t_nk_mc_group)
table(GSE123139$non_t_nk_mc_group)

sort(unique(GSE123139$Notes))

GSE123139=GSE123139%>%
  filter(treatment=="N") %>%
  filter(!all_mc_group %in%c(  "tumor")) %>%
  filter(!location %in%c(  "LN")) %>%
  filter(!Notes==c("Many empty wells","Mostly empty wells")) %>%
  filter(!non_t_nk_mc_group %in%c(  "osteoclast")) %>%
  filter(!t_nk_mc_group %in%c(  "naive")) %>%
  unfactor()


GSE123139=GSE123139[!(is.na(GSE123139$t_nk_mc_group) &
                        is.na(GSE123139$non_t_nk_mc_group)),]

GSE123139$major_celltypes=case_when(
  GSE123139$non_t_nk_mc_group== "B"~ "B",
  GSE123139$non_t_nk_mc_group== "DC"~ "DC",
  GSE123139$non_t_nk_mc_group== "macrophage"~ "Macrophage",
  GSE123139$non_t_nk_mc_group== "monocyte"~ "Monocytes",
  GSE123139$non_t_nk_mc_group== "non-classic-monocyte"~ "Monocytes",
  GSE123139$non_t_nk_mc_group== "mature-DC"~ "DC",
  GSE123139$non_t_nk_mc_group== "pDC"~ "pDC",
  GSE123139$non_t_nk_mc_group== "plasma"~ "Plasma",
  GSE123139$t_nk_mc_group== "cytotoxic"~ "T CD8",
  GSE123139$t_nk_mc_group== "dysfunctional"~ "T CD8",
  GSE123139$t_nk_mc_group== "NK"~ "NK",
  GSE123139$t_nk_mc_group== "Tfh"~ "T CD4",
  GSE123139$t_nk_mc_group== "Treg"~ "T CD4",
  GSE123139$t_nk_mc_group== "transitional"~ "T CD8",
  GSE123139$t_nk_mc_group== "dysf-cd4"~ "T CD4",
  GSE123139$t_nk_mc_group== "em-cd4"~ "T CD4",
  
  T~"NA"
  
)
table(
  GSE123139$major_celltypes)

GSE123139$minor_celltypes=case_when(
  GSE123139$non_t_nk_mc_group== "B"~ "B",
  GSE123139$non_t_nk_mc_group== "DC"~ "DC",
  GSE123139$non_t_nk_mc_group== "macrophage"~ "Macrophage",
  GSE123139$non_t_nk_mc_group== "monocyte"~ "Monocytes",
  GSE123139$non_t_nk_mc_group== "non-classic-monocyte"~ "Monocytes CD16",
  GSE123139$non_t_nk_mc_group== "mature-DC"~ "DC mature",
  GSE123139$non_t_nk_mc_group== "pDC"~ "pDC",
  GSE123139$non_t_nk_mc_group== "plasma"~ "Plasma",
  GSE123139$t_nk_mc_group== "cytotoxic"~ "TCD8 cytotoxic",
  GSE123139$t_nk_mc_group== "dysfunctional"~ "TCD8 dysfunctional",
  GSE123139$t_nk_mc_group== "NK"~ "NK",
  GSE123139$t_nk_mc_group== "Tfh"~ "TCD4 Th",
  GSE123139$t_nk_mc_group== "Treg"~ "TCD4 reg",
  GSE123139$t_nk_mc_group== "transitional"~ "TCD8 transitional",
  GSE123139$t_nk_mc_group== "dysf-cd4"~ "TCD4 dysfunctional",
  GSE123139$t_nk_mc_group== "em-cd4"~ "TCD4 memory", #em
  
  T~"NA"
  
)
table(GSE123139$minor_celltypes)


GSE123139$dataset="SKCM_GSE123139"


#######
# Put all neta files together as list
# Bar plot for contribution of each dataset to the major cell types
list_data=list(GSE176078, #breast
               GSE166555, #crc
               GSE139555, #kidney renal
               GSE131907, #luad
               GSE140228_drop, #hepato dropin
               GSE140228_smart, #hepato smartseq2
               GSE123139 #melanoma
)

#save datalist

saveRDS(list_data,file = paste0(dir,"curated_metalist.RDS"))

