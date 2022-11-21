# Get gene association from literature papers and find which genes were assigned to which cell type

library(foreach)
match_genes_to_celltype=function(our_list,lit_list){
  foreach(i=rownames(our_list),.combine=rbind) %do%{
  print(i)
  sources= foreach(j=seq_along(colnames(lit_list)),.combine=c)%do% {
    if (i %in% lit_list[,j,drop=T]){
      colnames(lit_list)[j]
    }else{NA}
  }  
  sources= na.omit(sources)
  if(length(sources)>0){
    paste0(sources,collapse = " ,")
  }else{c("-")}
  
  }
}
library(readxl)
  dirname= "N:/gene_lists_papers"
  gene_sign_file_path="N:./immunsig_scripts_manuscript/immune_signature_discovery_workflow/"
  gene_sig_file_name= "dbscan_cluster_genes_param_0.tsv"
  lit_gene_list= read_excel(path = paste0(dirname,"/arnon_melanoma_2018/NIHMS1508390-supplement-10.xlsx")
                            ,trim_ws = T,sheet = 2,skip = 1)
  gene_list=read.table(file = paste0(gene_sign_file_path,gene_sig_file_name),
                       sep = "\t",stringsAsFactors = F)   



  gene_list$arnon=match_genes_to_celltype(gene_list,lit_gene_list)
  
  lit_gene_list= read_excel(path = paste0(dirname,"/tirosh_melanoma_2016/tirosh_sign.xlsx")
                            ,trim_ws = T,sheet = 1,skip = 1)
  
  gene_list$tirosh=match_genes_to_celltype(gene_list,lit_gene_list)
  
  lit_gene_list= read_excel(path = paste0(dirname,"/li_melanoma_2018/li_sign.xlsx")
                            ,trim_ws = T,sheet = 1,skip = 1)
  
  gene_list$li=match_genes_to_celltype(gene_list,lit_gene_list)
  
  lit_gene_list= read_excel(path = paste0(dirname,"/zilionis_lung_2019/zilionis_sign.xlsx")
                            ,trim_ws = T,sheet = 1,skip = 1)
  
  gene_list$Zilionis=match_genes_to_celltype(gene_list,lit_gene_list)
#add go 
  db_genes_0= read.table(paste0(gene_sign_file_path,"dbscan_cluster_genes_param_0_GO_annot.tsv"),sep = "\t",stringsAsFactors = F)
  
  gene_list$our_annotation=db_genes_0$Go
  
  
  
  write.table(gene_list,paste0(gene_sign_file_path,"dbscan_cluster_genes_lit_annot_updated.tsv"),sep = "\t")
