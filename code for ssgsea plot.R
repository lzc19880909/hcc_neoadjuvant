library(ggplot2)
library(GSVA)
library(dplyr)
library(Hmisc)
library(ggpubr)
library(tidyverse)
library(reshape2)

###cellmarker_data###
cellmarker <- read.csv("pancancer_cellmarker.csv",header = T)
colnames(cellmarker)[2] <- "celltype" 
cellmarker <- split(as.matrix(cellmarker$Metagene) ,cellmarker$celltype)

###expression_data###
exp <- read.table("expression_data.txt",header = T)
exp_1 <- as.matrix(exp)
gsva_pancancer <- gsva(exp_1 , cellmarker , method = "ssgsea")

###add_group_info###
group.info <- read.table("group.info.txt",header = T)
idx <- str_replace_all(colnames(gsva_pancancer),"\\.","-")
group <- group.info$group[match(idx,group.info$sample)]

###add_cell_type###
typeidx <- read.table("ssgsea.type.txt",header = T,sep = "\t")
gsva_pancancer <- gsva_pancancer[match(typeidx$Cell.type,rownames(gsva_pancancer)),]
data_gsva = melt(gsva_pancancer)
colnames(data_gsva) = c("Celltype", "sample", "Composition")
data_gsva$sample <- str_replace_all(data_gsva$sample,"\\.","-")
gsva_df <- merge(data_gsva,group.info,by = "sample")
gsva_df$group <- factor(gsva_df$group,levels = c("pre-treatment","post-treatment"),ordered = T)

###plot_for_ssgsea###
ggplot(gsva_df, aes(x = Celltype, y = Composition))+ 
  labs(y="Enrichment Score", x = "", title = "",size = 14)+
  geom_rect(aes(xmin="Activated CD4 T cell",xmax="Type 2 T helper cell",ymin=-Inf,ymax=Inf),fill="grey95", alpha=0.9)+
  stat_summary(aes(group = group),fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", color = "black", position = position_dodge(0.8), width = 0.4, size = 0.5)+
  geom_boxplot(aes(fill = group), position = position_dodge(0.8), width = 0.8, outlier.alpha = 0, coef =0)+
  scale_fill_manual(values = c("#748ab3","#e8777b")) +
  stat_summary(aes(shape = group),fun.y = "mean", fun.args = list(mult = 1), geom = "point",position = position_dodge(0.8),size = 2,color = "black",fill = "white")+
  geom_jitter(aes(color = group),size = 0.8,shape = 20,position = position_jitterdodge(jitter.width = 0.6,dodge.width = 0.8)) +
  scale_color_manual(values = c("#1d4a9b","#e5171a")) +
  scale_shape_manual(values = c(21,21)) +
  ylim(0,1.2)+
  theme_bw() + 
  guides(fill=guide_legend(title ="Stage"),shape = "none",color = "none")+
  theme(plot.title = element_text(size = 12,color="black",hjust = -0.05), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title= element_text(size = 12),
        text = element_text(family = "serif")) + 
  stat_compare_means(aes(group =  group),
                     label = "p.signif",
                     hide.ns = T,size = 5)

ggsave("ssGSEA.pdf",height = 7,width = 10)