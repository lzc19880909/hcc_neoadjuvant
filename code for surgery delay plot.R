library(tidyverse)

surgery <- read.csv('~/plot/surgery.csv')

surgery <- surgery %>% 
  mutate(group1=as.character(order)) %>% 
  arrange(desc(day)) %>% 
  mutate(Patient = fct_reorder(factor(group1), day, .na_rm = TRUE),
         color = ifelse(day > 0,"high","low"),
         add = ifelse(Patient==5|Patient==15|Patient==19,day+1,NA),
         shape = as.factor(ifelse(!is.na(add),"Delay due to COVID-19",NA))) %>%
  group_by(group1)


p <- ggplot(surgery, aes(x=reorder(group1, order), day, shape = shape))+
      geom_point( aes(x=reorder(group1, order),add,shape = shape)) +
      scale_shape_manual(values = 8, na.translate = FALSE) +
      geom_bar(aes(fill = color), stat = "identity",width = 0.8, position = position_dodge(width = 0.5), na.rm=TRUE)+
      theme_classic()+
      theme(axis.text.x = element_text(angle=60,hjust = 1,vjust = 1), 
            legend.position = c(0.2,0.9),
            # legend.direction = "horizontal",
            legend.title = element_blank(),
            line = element_line(colour="grey",linewidth=1),
            text = element_text(family = "serif"),
            axis.text = element_text(size = 12.5),
            axis.title = element_text(size = 12.5))+
      labs(x="Patient",y="Days between treatment start and surgery")+
      scale_y_continuous(limits = c(-10, 55),breaks = c(-10, -5, 0, 5, 10, 15, 20, 25,30,35,40,45,50,55),labels = function(x){x+50})+
      geom_hline(yintercept = 0,lty="dashed", linewidth=0.1)+
      geom_hline(yintercept = 42,lty="dashed", linewidth=0.2,color="black")+
      scale_fill_manual(values = c("#4A708B","#6A5ACD"),
                        labels =c("Surgery > 50 days after treatment start","Surgery < 50 days after treatment start"),
                        na.translate = FALSE)+
      guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2))+
      annotate(geom = "text",x=1.5,y=40,label="+ 6 Weeks", family = "serif" ,size = 4)+
      annotate(geom = "text",x=19,y=37.5,label="Significant \nsurgery delay", family = "serif" ,size = 4)


p
ggsave("Surgery delay.pdf",p ,width=10,height=6)

