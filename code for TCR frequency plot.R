library(immunarch)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggthemes)
###path_sampleid_file###
metadata <- read.table("metadata.txt",header = T)
###TCR_clonetype_number_data###
data <- read.table("clonetypenumber.txt",header = T,sep = "\t")
meta <- cbind(metadata,data[match(metadata$sample_id,data$Sample),1:2])
meta <- arrange(meta,meta$patient,meta$stage)
meta$id <- substr(meta$sample_id,1,9)
meta$ID <- str_c("P",meta$patient,"_",meta$stage)
immdata$data <- lapply(1:nrow(meta), function(x){
  data <- read.table(file = meta$file_name[x],header = T)
  colnames(data)[c(1,5,6,7,3,4)] <- c("Clones", "V.name", "D.name", "J.name", "CDR3.nt", "CDR3.aa")
  data
})
attributes(immdata$data)$names <- meta$ID
immdata$meta <- meta[2:6]
immdata$meta$Sample <- immdata$meta$ID
###overlap_clonetypes_analysis###
over <- lapply(1:8, function(x){
  com <- intersect(unique(immdata$data[[names(immdata$data)[post][x]]]$CDR3.aa),unique(immdata$data[[names(immdata$data)[pre][x]]]$CDR3.aa))
  uni <- setdiff(unique(immdata$data[[names(immdata$data)[post][x]]]$CDR3.aa),unique(immdata$data[[names(immdata$data)[pre][x]]]$CDR3.aa))
  c(length(com),length(uni))
})
over <- as.data.frame(Reduce(rbind,over))
colnames(over) <- c("common","unique")
over$patient <- str_replace(names(immdata$data)[post],"_Post-treatment","")
over$patient <- str_replace(over$patient,"P","")
over <- arrange(over,desc(unique))
over.melt <- melt(over)
over.melt$variable <- factor(over.melt$variable,levels = c("unique","common"),labels = c("Newly","Overlap"),ordered = T)
over.melt$patient <- factor(over.melt$patient,levels = over$patient, ordered = T)
###overlap_clonetype_number_plot###
p1 <- ggplot(over.melt) +
  geom_bar(aes(x = patient, y = value, fill = variable,color = variable), position = "stack", stat = "identity", data = over.melt)+
  scale_fill_manual(values = c("#e3b6ad","#b05046"))+
  scale_color_manual(values = c("#e3b6ad","#b05046"))+
  theme_few()+
  theme(legend.position = c(0.8,0.89),
    legend.direction = "horizontal",
    axis.text.y=element_text(size=12, family = "serif"),  
    axis.title.y=element_text(size=14,family = "serif"),
    axis.text.x=element_text(size=12,family = "serif"))+
  labs(x = "", y = "Number of unique CDR3(aa) clones", title = , subtitle = "",legend = "", fill = "",color="")
###overlapped_clonetypes_frequency###
data <- lapply(1:8, function(x){
  com <- intersect(unique(immdata$data[[names(immdata$data)[post][x]]]$CDR3.aa),unique(immdata$data[[names(immdata$data)[pre][x]]]$CDR3.aa))
  uni <- setdiff(unique(immdata$data[[names(immdata$data)[post][x]]]$CDR3.aa),unique(immdata$data[[names(immdata$data)[pre][x]]]$CDR3.aa))
  overlap <- as.data.frame(cbind("overlap",com,names(immdata$data)[post][x]))
  colnames(overlap) <- c("Type","clone","ID")
  data <- immdata$data[[names(immdata$data)[post][x]]][immdata$data[[names(immdata$data)[post][x]]]$CDR3.aa %in% overlap$clone,]
  data
})
overlap <- list()
overlap$data <- data
attributes(overlap$data)$names <- meta$ID[post]
overlap$meta <- meta[post,2:6]
overlap$meta$Sample <- overlap$meta$ID
freq_over <- as.data.frame(repClonality(overlap$data, "homeo"))
freq2_over <- freq_over
freq2_over$`Large (1e-04 < X <= 0.01)` <- freq2_over$`Medium (1e-04 < X <= 0.001)`+freq2_over$`Large (0.001 < X <= 0.01)`
freq2_over <- cbind(freq2_over,overlap$meta[,2:3])
freq2_over.melt <- melt(freq2_over, id.vars = c("patient", "stage"), measure.vars = c("Rare (0 < X <= 1e-05)","Small (1e-05 < X <= 1e-04)","Large (1e-04 < X <= 0.01)","Hyperexpanded (0.01 < X <= 1)"))
freq2_over.melt$patient <- factor(freq2_over.melt$patient)
freq2_over.melt$stage <- factor(freq2_over.melt$stage, levels = c("Pre-treatment","Post-treatment"),ordered = T)
freq2_over.melt$type <- "Overlap"
###newly_clonetypes_frequency###
data2 <- lapply(1:8, function(x){
  uni <- setdiff(unique(immdata$data[[names(immdata$data)[post][x]]]$CDR3.aa),unique(immdata$data[[names(immdata$data)[pre][x]]]$CDR3.aa))
  newly <- as.data.frame(cbind("newly",uni,names(immdata$data)[post][x]))
  colnames(newly) <- c("Type","clone","ID")
  data <- immdata$data[[names(immdata$data)[post][x]]][immdata$data[[names(immdata$data)[post][x]]]$CDR3.aa %in% newly$clone,]
  data
})
newly <- list()
newly$data <- data2
attributes(newly$data)$names <- meta$ID[post]
newly$meta <- meta[post,2:6]
newly$meta$Sample <- newly$meta$ID
freq_new <- as.data.frame(repClonality(newly$data, "homeo"))
freq2_new <- freq_new
freq2_new$`Large (1e-04 < X <= 0.01)` <- freq2_new$`Medium (1e-04 < X <= 0.001)`+freq2_new$`Large (0.001 < X <= 0.01)`
freq2_new <- cbind(freq2_new,newly$meta[,2:3])
freq2_new.melt <- melt(freq2_new, id.vars = c("patient", "stage"), measure.vars = c("Rare (0 < X <= 1e-05)","Small (1e-05 < X <= 1e-04)","Large (1e-04 < X <= 0.01)","Hyperexpanded (0.01 < X <= 1)"))
freq2_new.melt$patient <- factor(freq2_new.melt$patient)
freq2_new.melt$stage <- factor(freq2_new.melt$stage, levels = c("Pre-treatment","Post-treatment"),ordered = T)
freq2_new.melt$type <- "Newly"
###overlap_newly_frequency_plot##
freq_all <- rbind(freq2_new.melt,freq2_over.melt)
freq_all$patient <- factor(freq_all$patient,levels = over$patient, ordered = T)
p2 <- ggplot(freq_all,aes(patient,value,fill = variable))+
  geom_bar(stat = "identity",position = "fill",color = "black",width = 0.8,size = 0.25)+
  scale_fill_manual(values = c("#808080","#e7e7e7","#C7C3DB","#4A4A90"))+
  facet_grid(rows = vars(type))+ theme_tufte() +
  scale_y_continuous(labels = scales::percent) + 
  guides(fill=guide_legend(title = "Clonetype group")) +
  theme(legend.position = "bottom", 
        strip.text = element_text(face = "bold",size = 14,family = "serif"),
        axis.text.y = element_text(size = 12,family = "serif"),
        axis.text.x = element_text(size = 12,family = "serif"),
        axis.title.x = element_text(size = 14,color = "black",family = "serif"),
        axis.title.y = element_text(size = 14,color = "black",family = "serif")) +
  labs(x = 'Patient', y = 'Relative abundance, perc')
###p1+p2###
p1/p2+plot_layout(heights = c(1,1.5))+plot_annotation(tag_levels = "A")&theme(plot.tag = element_text(size = 14,family = "serif",face = "bold",vjust = 1.5,hjust = -0.5))
ggsave("overlap_newly_frequency.pdf",width = 9,height = 10)
