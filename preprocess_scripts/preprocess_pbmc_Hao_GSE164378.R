library(Seurat)

library(dplyr)


sobj=Seurat::Read10X(data.dir = "C:/Users/X223465/Downloads/Hao_2021")
meta=read.csv("C:/Users/X223465/Downloads/Hao_2021/GSE164378_sc.meta.data_3P.csv",row.names = 1)
hao=CreateSeuratObject(counts = sobj,meta.data = meta)

hao[["percent.mt"]] <- PercentageFeatureSet(hao, pattern = "^MT-")
# VlnPlot(tiss_immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(hao$nFeature_RNA)#min 501, max 5990 
summary(hao$percent.mt) #0-15

# save RDS for referencing purposes
hao=NormalizeData(hao)

remove=hao$celltype.l2 !="CD4 CTL" &
  hao$celltype.l2 !="CD4 Proliferating" &
  hao$celltype.l2 !="gdT" &
  hao$celltype.l2 !="dnT" &
  hao$celltype.l2 !="Treg" &
  hao$celltype.l2 !="MAIT" &
  hao$celltype.l2 !="ILC" &
  hao$celltype.l2 !="Doublet"&
  hao$celltype.l2 !="Eryth" &
  hao$celltype.l2 !="Platelet" &
  hao$celltype.l2 !="HSPC"

hao=hao[,remove]
#
# sort(unique(hao$celltype.l2))

hao$harmonized_celltype=hao$celltype.l2
hao$harmonized_celltype=gsub("B intermediate","B",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("B memory","B",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("B naive","B",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("Plasmablast","B",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("cDC1","DC",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("cDC2","DC",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD14 Mono","Mono",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD16 Mono","Mono",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("NK_CD56bright","NK",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("NK Proliferating","NK",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD4 TCM","T CD4 memory_naive",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD4 TEM","T CD4 memory_naive",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD4 Naive","T CD4 memory_naive",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD8 TCM","T CD8",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD8 TEM","T CD8",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD8 Naive","T CD8",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("CD8 Proliferating","T CD8",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("ASDC","DC",hao$harmonized_celltype)
hao$harmonized_celltype=gsub("pDC","DC",hao$harmonized_celltype)
table(hao$harmonized_celltype)

saveRDS(hao,"Z:/pbmc_hao_ref.RDS")

