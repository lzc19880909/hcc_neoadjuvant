library(tidyverse)
library(circlize)
library(RColorBrewer)
library(data.table)
library(cols4all)
library(ggthemes)
library(reshape2)
library(rstatix)
library(ggpubr)
###expression_data###
data <- read.table("expression_data.txt",header = T,sep = "\t")
###add_group_info###
group.info <- read.table("group.info.txt",header = T)
colnames(data) <- str_replace_all(colnames(data),"\\.","-")
df <- data[,2:39]
group <- group.info$group[match(colnames(df),group.info$sample)]
patient <- group.info$patient[match(colnames(df),group.info$sample)]
ID <- str_c("P",patient,"_",group)
colnames(df) <- ID
info <- as.data.frame(cbind(group,patient,ID))
info$group <- factor(info$group,levels = c("pre-treatment","post-treatment"),ordered = T)
info$patient <- as.numeric(info$patient)
info <- arrange(info,group,patient)
idx <- match(info$ID,colnames(df))
df <- df[,idx]
df$gene <- data$Symbol
###add_immunegenes_info###
immune <-read.table("immune.genelist.txt",header = F,sep = "\t") 
colnames(immune) <- c("cell","gene")
tdf <- df[df$gene %in% immune$gene,]
rownames(tdf) <- tdf$gene
tdf <- tdf[,1:38]
exp <- tdf
exp <- exp[match(immune$gene,rownames(exp)),]
###add_phenotype_data###
phe <-fread("pheno.txt",header = T)
pdata <- phe[,c(1,16,17)]
colnames(pdata) <- c("patient","CLD","mRECIST")
#CLD:"changes in longest diameter"
pheno <- merge(info,pdata)
pheno <- pheno[,c(3,1,2,5,4)]
colnames(pheno) <- c("ID","Patient","Stage","mRECIST","CLD")
pheno <- arrange(pheno,Stage,Patient)
###calculate immune score###
df_sum <- cbind(exp,immune)
data.sum <- lapply(1:38, function(x){
  data2 <- df_sum[,c(x,39)]
  colnames(data2) <- c("patient","cell")
  stat <- data2 %>% group_by(cell) %>% summarise(immunescore = sum(patient))
  stat$ID <- colnames(df_sum)[x]
  stat
})
score <- Reduce(rbind,data.sum)
score$immunescore <- log2(score$immunescore+1)
score_merge <- merge(score,pheno,by = "ID")
score_merge <- score_merge[score_merge$cell %in% c("Tumor reactive T cells","Cytotoxic cells"),]
score_merge$cell <- factor(score_merge$cell,levels = c("Tumor reactive T cells","Cytotoxic cells"),ordered = T)
###compared by treatment###
colors <- brewer.pal(8,"Set1")
ggplot(score_merge,aes(Stage, immunescore, fill = Stage))+
  geom_line(aes(group = Patient),
            size = 0.25,color="grey55")+#图层在下，就不会显示到圆心的连线
  geom_point(shape = 21,
             size = 2,
             stroke = 0.6,color="black")+
  scale_x_discrete(expand = c(-1.05, 0)) + # 坐标轴起始
  scale_fill_manual(values = c(colors[2],colors[1]))+
  facet_wrap(~cell,nrow = 2, scales = "free")+
  stat_compare_means(method = "t.test",paired = TRUE,comparisons=list(c("pre-treatment", "post-treatment")),label = "p.format",size =4)+
  geom_rangeframe() + # 坐标轴分离
  theme_tufte() +
  theme(legend.position = c(0.9,0.1),
        strip.text = element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14, color = "black",face = "bold")) +
  labs(x = ' ', y = 'Log2(Immune score)')
ggsave("immunescore.treatment.pdf",width = 4,height = 6)

###compared by mRECIST###
score_merge$mRECIST <- factor(score_merge$mRECIST,levels = c("SD","PR","CR"),ordered = T)
ggplot(score_merge[score_merge$Stage == "post-treatment",],aes(mRECIST, immunescore, fill = mRECIST))+
  geom_boxplot()+
  geom_point(shape = 21,
             size = 2,
             stroke = 0.6,color="black")+
  scale_fill_manual(values = c("#61739f","#c8d8ed","#efe8e5"))+
  facet_wrap(~cell,nrow = 2, scales = "free")+
  stat_compare_means(method = "wilcox.test", comparisons=list(c("CR", "PR"),c("CR", "SD"),c("SD", "PR")),size = 3)+
  geom_rangeframe() + 
  theme_tufte() +
  theme(legend.position = c(0.9,0.1),  
        strip.text = element_text(face = "bold",size = 10),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14, color = "black", face = "bold")) +
  labs(x = ' ', y = 'Log2(Immune score)')
ggsave("immunescore.mRECIST.pdf",width = 6,height = 6)
