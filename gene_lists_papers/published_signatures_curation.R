
library(foreach)
library(ComplexHeatmap)
library(bayesbio)
library(readxl)
library(dplyr)

#Script for calculating jaccard index between our signatures and published signatures

lit_dirname= "N:/immunsig_scripts_manuscript/gene_lists_papers/"
dirname="N:./immunsig_scripts_manuscript/immune_signature_discovery_workflow/"



##################################################################
########### Put all published signatures together ################

#Abbas et al. 2005 
abbas= read_excel(path = paste0(lit_dirname,"abbas_2005.xlsx"),sheet = 1)

abbas_filtered=abbas %>%
  select(Cell_type,Gene) %>%
  filter(!Cell_type %in% c("Multiple","Lymphoid","Myeloid"))

abbas_filtered=cbind(abbas_filtered,data.frame(source="Abbas_2005"))

#Becht et al. 2016
becht=read.table(paste0(lit_dirname,"becht.txt"),sep = "\t",
                 stringsAsFactors = F,header=T)
becht_filter=becht%>%
  select(Cell.population,HUGO.symbols)
colnames(becht_filter)=c("Cell_type","Gene")
becht_filter=cbind(becht_filter,data.frame(source="Becht_2016"))

#Charoentong et al. 2016
Charoentong= read_excel(path = paste0(lit_dirname,"/Charoentong.xlsx"),sheet = 1)
Charoentong_f= foreach(i=seq_along(rownames(Charoentong)),.combine = rbind) %do% {
  genes= Charoentong[i,-c(1,2)]
  genes=genes[!is.na(genes)]
  cell_type= Charoentong[i,1]
  data.frame(Cell_type=cell_type,Gene=genes)
}
Charoentong_f=cbind(Charoentong_f,data.frame(source="Charoentong_2017"))
colnames(Charoentong_f)[1]="Cell_type"

#Bindea et al. 2013
bindea= read_excel(path = paste0(lit_dirname,"/bindea.xlsx"),sheet = 1)
bindea_f= foreach(i=seq_along(rownames(bindea)),.combine = rbind) %do% {
  genes= bindea[i,-c(1,2)]
  genes=genes[!is.na(genes)]
  cell_type= bindea[i,1]
  data.frame(Cell_type=cell_type,Gene=genes)
}
bindea_f=cbind(bindea_f,data.frame(source="Bindea_2013"))
colnames(bindea_f)[1]="Cell_type"


# Other published signatures were previously summoned up into published_sign.txt
published=read.table(paste0(lit_dirname,"published_sign.txt"),sep = "\t",
                     stringsAsFactors = F,header=T)
# Newman et al. 2015
cibersort=published%>%
  filter(meta_source=="cibersort") %>%
  select(id_geneset_name,gene_symbol)

cell_types=gsub(pattern = ".*_lm22_(.*)_25822800",x = cibersort$id_geneset_name,replacement = "\\1")
cibersort_f=data.frame(Cell_type=cell_types,Gene=cibersort$gene_symbol,source="Newman_2015")

# Nirmal et al. 2018
imsig=published%>%
  filter(meta_source=="imsig") %>%
  select(id_geneset_name,gene_symbol)
cell_types=gsub(pattern = "imsig_(.*)",x = imsig$id_geneset_name,replacement = "\\1")

imsig_f=data.frame(Cell_type=cell_types,Gene=imsig$gene_symbol,source="Nirmal_2018")

# Angelova et al. 2015
angelova=published%>%
  filter(meta_publication %in% c("https://doi.org/10.1186/s13059-015-0620-6")) %>%
  select(id_geneset_name,gene_symbol)

cell_types=gsub(pattern = "angelova2015[.]",x = angelova$id_geneset_name,replacement = "")
angelova_f=data.frame(Cell_type=cell_types,Gene=angelova$gene_symbol,source="Angelova_2015")

all_sig= rbind(abbas_filtered,angelova_f,bindea_f,cibersort_f,becht_filter,Charoentong_f,imsig_f)
write.table(all_sig,file = paste0(lit_dirname,"all_immune_cell_signatures.txt"),sep = "\t")

##################################################################
########### Curate signatures manually ###########################

# Cell type names are manually curated to have common nomenclature
# all_immune_cell_signatures_curated.txt is the curated version of
# all_immune_cell_signatures.txt

lit_gene_list=read.table( paste0(lit_dirname,"all_immune_cell_signatures_curated.txt"),sep = "\t",header = T)
