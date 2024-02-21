library(tidyverse)

assessment <- read.csv('~/plot/assessment.csv')

df_long <- assessment %>% 
  pivot_longer(3:5,names_to = "type",values_to = "value") %>% 
  mutate(type=factor(type, 
                     levels = c("RECIST","mRECIST","pathologic.response")),
         value=ifelse(value>=0,0-value,value))  %>%
  mutate(group1=paste0("Patient ",as.character(patient.ID))) %>%
  mutate(value=ifelse(group1=='Patient 1' & type=='pathologic.response',0,value)) %>%
  filter(group1!='Patient 1')



p <- ggplot(df_long,aes(x = reorder(patient.ID, order),y=value,fill=type)) + 
    geom_bar(stat='identity',position = position_dodge(width = 0.5), width = 0.6,na.rm=TRUE) +
    theme_classic()+
    theme(legend.position = "top",
          legend.justification = "left",
          legend.key.size = unit(0.5, "cm"),
          axis.text.x = element_text(angle=30,hjust = 1,vjust = 1),
          legend.title = element_blank(),
          legend.direction = 'vertical',
          axis.text = element_text(size = 12.5),
          axis.title = element_text(size = 12.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,'cm'))+
    labs(x="Patient",y="Response (%)")+
    scale_fill_manual(values = c("#C6E2FF","#8F8FBD", "#2E8B57"), labels=c("Percentage change in tumor size (RECIST v1.1) from baseline","Percentage change in tumor size (mRECIST) from baseline","Pathological response"))+
    scale_y_continuous(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2), labels = function(x){x*100}) +
    geom_hline(yintercept = -0.7,lty="dashed",linewidth=0.2) +
    geom_hline(yintercept = -0.3,lty="dashed",linewidth=0.2) +
    geom_hline(yintercept = 0,lty="dashed", linewidth=0.1) +
    annotate(geom = "text",x=0.95,y=-0.775,label="-70%") +
    annotate(geom = "text",x=0.95,y=-0.375,label="-30%") + 
    scale_x_discrete(breaks=c('2-SII','2-SIV','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18', '19', '20'),
                     labels=c('','','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18', '19', '20')) +
    coord_cartesian(clip = 'off',ylim=c(-1,0.2),xlim=c(1,20)) +
    annotate('text',x=1.5,y=-1.12,label='2',size=5,angle=30,color='#696969', family='sans') +
    annotate('segment',x=0.95,xend=2.05,y=-1.085,yend=-1.085,color='black',cex=.5)+
    annotate('text',x=0.5,y=-1.175,label="Segment[II]",parse=T,size=4,angle=30,color='#696969', family='sans') +
    annotate('text',x=1.5,y=-1.175,label="Segment[IV]",parse=T,size=4,angle=30,color='#696969', family='sans') 
p
ggsave("Waterfall plot of tumour regression.pdf",p,width=10,height=6)

  

