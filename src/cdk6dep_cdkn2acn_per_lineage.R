INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(dplyr)
library(readr)
library(ggplot2)

CDK6dep_CDKN2ACN<-read.csv(file.path(INPUT_PATH, "CDK6dep_CDKN2ACN.csv"))

if (file.exists(file.path(OUTPUT_PATH, "CDK6_CDKN2A_LinA.csv"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_CDKN2A_LinA.csv"))}
if (file.exists(file.path(OUTPUT_PATH, "CDK6_CDKN2A_LinB.csv"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_CDKN2A_LinB.csv"))}
if (file.exists(file.path(OUTPUT_PATH, "CDK6_CDKN2A_LinC_q.csv"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_CDKN2A_LinC_q.csv"))}
                           
CDK6_CDKN2A_LinA<-CDK6dep_CDKN2ACN %>%
  group_by(Lineage) %>%
  summarize(
    coef = cor.test(CDK6, CDKN2A)$estimate,
    p_value = cor.test(CDK6, CDKN2A)$p.value
  )
write_csv(CDK6_CDKN2A_LinA, file.path(OUTPUT_PATH,"CDK6_CDKN2A_LinA.csv"))

CDK6_CDKN2A_LinB<-CDK6dep_CDKN2ACN %>%
  group_by(Subtype2) %>%
  summarize(
    coef = cor.test(CDK6, CDKN2A)$estimate,
    p_value = cor.test(CDK6, CDKN2A)$p.value
  )
write_csv(CDK6_CDKN2A_LinB, file.path(OUTPUT_PATH,"CDK6_CDKN2A_LinB.csv"))

CDK6_CDKN2A_LinC<-read.csv(file.path(INPUT_PATH,"CDK6_CDKN2A_LinC.csv"))

CDK6_CDKN2A_LinC_q<-p.adjust(CDK6_CDKN2A_LinC$p_value,method="BH")
write.csv(CDK6_CDKN2A_LinC_q,file.path(OUTPUT_PATH,"CDK6_CDKN2A_LinC_q.csv"))

ggplot(CDK6_CDKN2A_LinC,aes(x=coef,y=log10P))+geom_point(stroke=2, aes(shape=class,fill=color,colour=color,alpha=significance,size=significance))+scale_shape_manual(values=c(21,1))+scale_color_manual(values=c("#011892","#ff7e79", "#ff9200","#005492","#0095ff","#018e00","#ff2600","#521b92","#931100","#918f00","#999999"))+scale_fill_manual(values=c("#011892","#ff7e79", "#ff9200","#005492","#0095ff","#018e00","#ff2600","#521b92","#931100","#918f00","#999999"))+scale_alpha_manual(values=c(0.4,0.8))+scale_size_manual(values=c(5,10))+theme(aspect.ratio=1)
