library("ggsci")
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(tidyr)
library(cowplot)
library(varhandle)
library(gridExtra)

save_dir="N:./immunsig_scripts_manuscript/benchmarking/"
dir.create(save_dir)

# Kotliarov
tab=read.table(paste0(save_dir,"kotliarov/prediction_scores_pbmc_kotliarov.txt"),header = T,
               sep = "\t")
a=tab$Sensitivity
tab$Sensitivity[1]

wo_space=sapply(a,FUN = function(t){t=unfactor(t)
nm=as.numeric(strsplit(t, " +")[[1]])/100
data.frame(nm)
})
tab[,c(3:ncol(tab))]=do.call(rbind,wo_space)

tab$method_ref=paste0(tab$Method,"_",tab$Reference)
tab$method_ref=as.factor(tab$method_ref)
tab$Dataset="Kotliarov"
m_tab=melt(tab)

m_tab$complete=paste0(m_tab$variable,"_",m_tab$method_ref)
m_tab$Method %>% unique()

m_tab_filtered=m_tab %>%
  filter(Reference== "None" | Reference=="Hao" ) 

figa_kot= subset(m_tab_filtered,Method=="Seurat 200 HVGs"|Method=="Seurat 2k HVGs"|
                         Method== "RF our 199 genes" | Method=="singleR" 
)
figa_kot$Method <- factor(figa_kot$Method, levels = c("Seurat 2k HVGs",  
                                                                  "RF our 199 genes",
                                                                  "singleR" ,
                                                                  "Seurat 200 HVGs"))

figb_kot= subset(m_tab_filtered,Method=="RF our 199 genes"|
                         Method=="RF random 199 genes"|
                         Method=="RF Abbas"|
                         Method=="RF Abbas"|
                         Method== "RF Angelova" | 
                         Method=="RF Charoentong" |
                         Method=="RF Nieto"
)
figb_kot$Method <- factor(figb_kot$Method, levels = c("RF our 199 genes", 
                                                                  "RF Charoentong",
                                                                  "RF Angelova" ,"RF Nieto",
                                                                  "RF Abbas",
                                                                  "RF random 199 genes"))



plot_a_kot=figa_kot %>%
  ggplot( aes(variable, value,fill=Method))+
  geom_bar(position=position_dodge(0.9), colour="black", stat="identity", width=0.9) +
  
  #geom_bar(stat = "identity", aes(fill = Method), position = "dodge",width = 0.2, binwidth=0)+theme_classic()+
  labs(x="Metrics", y = "Metrics score")+ theme_bw(base_size=22)+coord_flip()+
  scale_fill_nejm()+ guides(shape = guide_legend(override.aes = list(size = 0.5)))+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))+facet_wrap("Dataset")+ theme(legend.position="top")+guides(fill=guide_legend(nrow=2,byrow=TRUE))

plot_b_kot=figb_kot %>%
  ggplot( aes(variable, value,fill=Method))+
  geom_bar(position=position_dodge(0.9), colour="black", stat="identity", width=0.9) +
  
  #geom_bar(stat = "identity", aes(fill = Method), position = "dodge",width = 0.2, binwidth=0)+theme_classic()+
  labs(x="Metrics", y = "Metrics score")+ theme_bw(base_size=22)+coord_flip()+
  scale_fill_nejm()+ guides(shape = guide_legend(override.aes = list(size = 0.5)))+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))+facet_wrap("Dataset")+ theme(legend.position="top")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

unique(m_tab$Method)
m_tab_filtered= subset(m_tab,Method=="Seurat 2k HVGs"|
                         Method== "RF our 199 genes" | Method=="singleR" 
)
# svg(paste0(save_dir,"SuppFig5a.svg"),
#     width=20,height = 6)
plot=m_tab_filtered %>%
  ggplot( aes(variable, value,fill=Reference))+
  geom_bar(position=position_dodge(0.9), colour="black", stat="identity", width=0.9) +
  
  #geom_bar(stat = "identity", aes(fill = Method), position = "dodge",width = 0.2, binwidth=0)+theme_classic()+
  labs(x="Metrics", y = "Metrics score")+ theme_bw(base_size=22)+coord_flip()+
  scale_fill_nejm()+facet_nested_wrap("Method")+ guides(shape = guide_legend(override.aes = list(size = 0.5)))+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),legend.position="top",
        panel.spacing = unit(2, "lines"))
# plot
# dev.off()


###sup fig. all together
plot_zheng=readRDS(paste0("N:./immunsig_scripts_manuscript/benchmarking/zheng_sorted2k/",
                  "SuppFig5B_seurat_diffHVGs_zhengsorted_ref_hao.RDS"))
plot_seurat=readRDS(paste0("N:./immunsig_scripts_manuscript/benchmarking/kotliarov/seurat/",
                          "SuppFig5B_seurat_diffHVGs_kotliarov_ref_hao.RDS"))
svg(paste0(save_dir,"SuppFig3.svg"),
    width=13,height = 10)
ggarrange(plot,ggarrange(plot_zheng,plot_seurat,nrow=1,
                               labels = c("B", "C"),
                         font.label = list(size = 22)),nrow=2,labels = c("A", ""),
          font.label = list(size = 22))
dev.off()
#####
#Zheng

tab=read.table(paste0("N:./immunsig_scripts_manuscript/benchmarking/zheng_sorted2k/",
                      "prediction_scores_pbmc_zheng_sorted_2k.txt"),header = T,
               sep = "\t")
a=tab$Sensitivity
tab$Sensitivity[1]

wo_space=sapply(a,FUN = function(t){t=unfactor(t)
nm=as.numeric(strsplit(t, " +")[[1]])/100
data.frame(nm)
})
tab[,c(3:ncol(tab))]=do.call(rbind,wo_space)

tab$method_ref=paste0(tab$Method,"_",tab$Reference)
tab$method_ref=as.factor(tab$method_ref)
tab$Dataset="Zheng"

m_tab=melt(tab)

m_tab$complete=paste0(m_tab$variable,"_",m_tab$method_ref)
m_tab$Method %>% unique()


#m_tab$Method=unfactor(m_tab$Method)
m_tab$Method <- factor(m_tab$Method, levels = c(  "Seurat 2k HVGs" ,
                                                  "RF our 182 genes",
                                                  "RF Charoentong",
                                                  "singleR",
                                                  "Seurat 180 HVGs" ,
                                                  "RF Angelova",
                                                  "RF Nieto",
                                                  "RF Abbas",
                                                  "RF random 182 genes"
                                                  
))


figa_zheng= m_tab[grepl("Seurat",m_tab$Method) |
                    grepl("our",m_tab$Method)|
                    grepl("singleR",m_tab$Method)  ,]
plot_a_zheng=figa_zheng %>%
  ggplot( aes(variable, value,fill=Method))+
  geom_bar(position=position_dodge(0.9), colour="black", stat="identity", width=0.9) +
  
  #geom_bar(stat = "identity", aes(fill = Method), position = "dodge",width = 0.2, binwidth=0)+theme_classic()+
  labs(x="", y = "Metrics score")+ theme_bw(base_size=22)+coord_flip()+
  scale_fill_nejm()+ guides(shape = guide_legend(override.aes = list(size = 0.5)))+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))+facet_wrap("Dataset")+ theme(legend.position="top")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

figb_zheng= m_tab[!(grepl("Seurat",m_tab$Method)|
                      grepl("singleR",m_tab$Method)),]

plot_b_zheng=figb_zheng %>%
  ggplot( aes(variable, value,fill=Method))+
  geom_bar(position=position_dodge(0.9), colour="black", stat="identity", width=0.9) +
  
  #geom_bar(stat = "identity", aes(fill = Method), position = "dodge",width = 0.2, binwidth=0)+theme_classic()+
  labs(x="", y = "Metrics score")+ theme_bw(base_size=22)+coord_flip()+
  scale_fill_nejm()+ guides(shape = guide_legend(override.aes = list(size = 0.5)))+ 
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))+facet_wrap("Dataset")+ theme(legend.position="top")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


svg(paste0(save_dir,"Fig3.svg"),
    width=12,height = 7)
grid.arrange(plot_a_kot,plot_a_zheng,nrow=1)

dev.off()

svg(paste0(save_dir,"Fig4.svg"),
    width=20,height = 10)
grid.arrange(plot_b_kot,plot_b_zheng,nrow=1)

dev.off()


