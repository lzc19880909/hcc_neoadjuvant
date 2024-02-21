library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(data.table)
library(cols4all)
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
exp <- apply(tdf, 1, scale)
rownames(exp) <- colnames(tdf)
exp <- t(exp)
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
colors <- brewer.pal(8,"Set1")
col_fun = colorRamp2(c(-2, 0, 1.5),c(colors[2],"white","#c13b2f"))
pheno$Stage <- factor(pheno$Stage,levels = c("pre-treatment","post-treatment"),labels = c("Pre-treatment","Post-treatment"),ordered = T)
pheno$mRECIST <- factor(pheno$mRECIST,levels = c("SD","PR","CR"),ordered = T)
pheno$Patient <- factor(as.character(pheno$Patient),levels <- as.character(c(2:20)),ordered = T)
pheno <- arrange(pheno,Stage,mRECIST,Patient)
exp <- exp[,match(pheno$ID,colnames(exp))]
col_mRECIST <- c("#61739f","#c8d8ed","#efe8e5")
attributes(col_mRECIST)$names <- levels(pheno$mRECIST)
col_stage <- c( "#426a9a","#c13b2f")
attributes(col_stage)$names <- levels(pheno$Stage)
color2 <- brewer.pal(8,"Pastel1")
col_fun2 = colorRamp2(c(0,20,40), c("#aac2dc","white","#e5a4a4"))
col_cld <- rep("darkgrey",38)
col_cld[pheno$CLD==-1] <- "#d82d1e"
col_patient <- col_fun2(seq(2, 20))
attributes(col_patient)$names <- levels(pheno$Patient)
###heatmap_plot###
ha <- HeatmapAnnotation(Stage = anno_block(gp = gpar(fill = col_stage,col = col_stage),
                                           labels = c("Pre-treatment", "Post-treatment"), 
                                           labels_gp = gpar(col = "white", fontsize = 12,font = 2)),
                        mRECIST = pheno$mRECIST, 
                        CLD = anno_points(pheno$CLD,gp = gpar(col = col_cld)),
                        Patient = anno_simple(1:38,pch = as.character(pheno$Patient),col = col_fun2,pt_gp = gpar(col = "black"), pt_size = unit(3, "mm")),
                        col = list(mRECIST = col_mRECIST, Patient = col_patient, Stage = col_stage), 
                        show_annotation_name = TRUE,
                        annotation_name_side = "left",
                        annotation_name_gp = gpar(fontsize = 10,font = 2),
                        annotation_height = c(1.25, 1, 2.5,1),
                        gap = unit(c(0.2,0.1,0.1,0.1), "cm"))

ht_opt(heatmap_border = TRUE)
ht <- Heatmap(exp, cluster_rows = F, cluster_columns = F,col = col_fun,row_split = immune$cell,column_split = pheno$Stage,column_title = NULL,row_names_gp = gpar(fontsize=10,font=3),
              top_annotation = ha,show_column_names = F,show_row_names = T,row_title_rot= 0,row_title_gp = gpar(fontsize =10, font = 2),
              heatmap_legend_param = list(title = "Expression"))
draw(ht,merge_legend = TRUE)
pdf(file="heatmap.pdf",width = 10,height = 7)
draw(ht,merge_legend = TRUE)
dev.off()