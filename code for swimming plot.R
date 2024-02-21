library(tidyverse)

date <- read.csv('~/plot/date.csv')

recist <- read.csv('~/plot/recist.csv')

post <- read.csv('~/plot/post.csv')

recist <- recist %>%
          mutate(type=paste0(response," by ",cata),
                 Patient = factor(No.)) %>%
          group_by(Patient,cata) %>%
          top_n(1,date)

post <- post %>%
        mutate(Patient = factor(No.)) %>%
        group_by(Patient)

date <- date %>%
        mutate(firstdate=ifelse(date.of.the.1st.PD.1<=date.of.the.1st.SBRT,date.of.the.1st.PD.1,date.of.the.1st.SBRT),
               firstpdl1=ifelse((date.of.the.1st.PD.1-firstdate)/6*2 > 2, (date.of.the.1st.PD.1-firstdate-6)/30.4375+2, (date.of.the.1st.PD.1-firstdate)/6*2),
               firstsbrt=ifelse((date.of.the.1st.SBRT-firstdate)/6*2 > 2, (date.of.the.1st.SBRT-firstdate-6)/30.4375+2, (date.of.the.1st.SBRT-firstdate)/6*2),
               secondsbrt=ifelse((date.of.the.2nd.SBRT-firstdate)/6*2 > 2, (date.of.the.2nd.SBRT-firstdate-6)/30.4375+2, (date.of.the.2nd.SBRT-firstdate)/6*2),
               thirdsbrt=ifelse((date.of.the.3rd.SBRT-firstdate)/6*2 > 2, (date.of.the.3rd.SBRT-firstdate-6)/30.4375+2, (date.of.the.3rd.SBRT-firstdate)/6*2),
               secondpdl1=ifelse((date.of.the.2nd.PD.1-firstdate)/6*2 > 2, (date.of.the.2nd.PD.1-firstdate-6)/30.4375+2, (date.of.the.2nd.PD.1-firstdate)/6*2),
               surgery=(date.of.surgery-firstdate-6)/30.4375+2,
               recurrent=(date.of.the.diagnosis.of.recurrence-firstdate-6)/30.4375+2,
               retreatment=(date.of.the.treatment.of.recurrence-firstdate-6)/30.4375+2,
               duration=(survival.date-firstdate-6)/30.4375+2,
               treaton =ifelse(ongoing==1,duration+0.15,''),
               rfa=(RFA-firstdate-6)/30.4375+2,
               # Patient = fct_reorder(factor(order), duration),
               Patient=order,
               line = 'Treatment ongoing') %>%
        group_by(Patient)

data_recist <- merge(date,recist,by = "Patient") %>%
                mutate(value=(date-firstdate-6)/30.4375+2) %>%
                dplyr::select(Patient,type,value) %>%
                filter(value >= 2)

data_post <- merge(date,post,by = "Patient") %>%
              mutate(value=(date-firstdate-6)/30.4375+2) %>%
              dplyr::select(Patient,type,value) %>%
              filter(value >= 2)

df_long <- date %>% 
  pivot_longer(c(15:22,24,25),names_to = "type",values_to = "value") %>% 
  mutate(type=factor(type)) %>%
  dplyr::select(Patient,type,value) %>%
  bind_rows(data_recist) %>%
  bind_rows(data_post) %>%
  mutate(type_1=case_when(type=="recurrent"~"Recurrent diagnosis",
                          type=="retreatment"~"Recurrent treatment",
                          type=="secondpdl1"~"neoadjuvant PD-1",
                          type=="secondsbrt"~"SBRT",
                          type=="surgery"~"Surgery",
                          type=="thirdsbrt"~"SBRT",
                          type=="treaton"~"Treatment ongoing",
                          type=="firstpdl1"~"neoadjuvant PD-1",
                          type=="firstsbrt"~"SBRT",
                          type=="rfa"~"RFA",
                          type=="post surgery treatment"~"adjuvant PD-1",
                          type=="follow-up"~"Follow-up",
                          TRUE~type)) %>%
  filter(type_1!='Follow-up' & type_1!='Treatment ongoing') %>%
  mutate(Patient1=case_when(type=="firstpdl1"~21-as.numeric(Patient)+0.2,
                            type_1=="SBRT"~21-as.numeric(Patient)-0.2,
                            grepl("RECIST v1.1",type_1)~21-as.numeric(Patient)+0.2,
                            grepl("mRECIST",type_1)~21-as.numeric(Patient)-0.2,
                            TRUE~21-as.numeric(Patient)))

date$Patient<-factor(date$Patient,levels = c('20','19','18','17','16','15','14','13','12','11','10','9','8','7','6','5','4','3','2','1'))
  
p <- ggplot(data = date, aes(x=duration, y=Patient))+
  geom_bar(stat = "identity", fill= "#C6E2FF", width = 0.5, show.legend = FALSE)+
  theme_classic()

p

df_long$type_1 <- factor(df_long$type_1, levels=c('neoadjuvant PD-1', 'adjuvant PD-1', 'SBRT', 'Surgery', 'RFA', 'Recurrent diagnosis', 'Recurrent treatment', 'CR by mRECIST', 'PR by mRECIST', 'SD by mRECIST', 'CR by RECIST v1.1', 'PR by RECIST v1.1', 'SD by RECIST v1.1'))

p1<-p+geom_point(data = df_long, aes(value, Patient1, shape=type_1, color=type_1), size=2.5,na.rm = TRUE)+
  scale_colour_manual(values = c("#700699","#DA3EA7","#1A900B","#CD8599","#8470FF","#8B658B","#8470FF","#5F9EA0","#8B0000","#CD8500","#156759","#8B0000","#000080"))+
  #scale_fill_manual(labels =c("CR by mRECIST","PR by mRECIST","PR by RECIST","Recurrent diagnosis","Recurrent treatment",
                             #"SD by mRECIST","SD by RECIST","2nd PD-1 treatment","2nd SBRT","Surgery","3rd SBRT","Treatment ongoing"))+
  scale_shape_manual(values= c(21,21,2,3,8,25,25,23,23,23,15,15,62))+
  geom_segment(data = date %>% dplyr::filter(ongoing == 1),
               aes(x=duration+0.1, xend=duration+0.3, y=Patient, yend=Patient, linetype = line),
               linewidth = 0.75, color="#8470FF",
               arrow = arrow(length = unit(0.15,"cm"))) +
  guides(linetype = guide_legend('test', order = 2), color = guide_legend('a', order = 1), shape = guide_legend('a', order = 1))+
  labs(x="Time since treatment start",y='Patient')+
  geom_vline(xintercept = 2,lty="dashed", linewidth=0.2,color="black")+
  annotate(geom = "text",x=1,y=20.75,label="Week 1", family = "serif" ,size = 3.5)+
  scale_x_continuous(limits = c(-0.25, 24),expand=c(0,0),breaks = c(0,0.333,0.667,1,1.333,1.667,2,3.9,5.9,7.9,9.9,11.9,13.9,15.9,17.9,19.9,21.9,23.9),
                     labels=c('D1','D2','D3','D4','D5','D6','D7','M2','M4','M6','M8','M10','M12','M14','M16','M18','M20','M22'))+
  scale_y_discrete(expand = c(0,1))+
  theme_classic()+
  theme(legend.key.size = unit(12, "point"),
        legend.title = element_blank(),
        legend.position = c(0.875,0.35),
        line = element_line(colour="#000000",linewidth=0.5),
        text = element_text(family = "serif",size=11),
        axis.text.x = element_text(angle=60,hjust = 1,vjust = 1),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12.5))
  # annotate("segment", x = 18.95, xend = 19.15, y = 3.55, yend = 3.55,
  #          linewidth = 0.3, arrow = arrow(length = unit(0.1,"cm")))
p1


ggsave("Swimming plot of treatment time.pdf",p1,width=10,height=6.57)

