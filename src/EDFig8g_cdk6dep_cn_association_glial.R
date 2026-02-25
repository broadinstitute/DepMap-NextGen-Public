INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

CDK6vsCN_g<-read.csv(file.path(INPUT_PATH, "CDK6vsCN_glial.csv"))

if (file.exists(file.path(OUTPUT_PATH, "CDK6_CNassociation_glial.txt"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_CNassociation_glial.txt"))}
if (file.exists(file.path(OUTPUT_PATH, "CDK6_CNassociation_gq.csv"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_CNassociation_gq.csv"))}
                     
for (j in 7:ncol(CDK6vsCN_g)){
  corr1<-cor.test(CDK6vsCN_g[,j],CDK6vsCN_g[,6])
  corp1<-corr1$p.value
  core1<-corr1$estimate
  cat(j,corp1,core1,"\n",file=file.path(OUTPUT_PATH, "CDK6_CNassociation_glial.txt"),append=TRUE)
}

CDK6_CNassociation_g<-read.csv(file.path(INPUT_PATH, "CDK6_CNassociation_glial.csv"))

CDK6_CNassociation_gq<-p.adjust(CDK6_CNassociation_g$Corr_p,method="BH")
write.csv(CDK6_CNassociation_gq,file.path(OUTPUT_PATH, "CDK6_CNassociation_gq.csv"))


ggplot(CDK6_CNassociation_g,aes(x=Corr,y=log10P))+geom_point(aes(alpha = Significance, size = Significance,colour=EnrichedLoc))+scale_alpha_manual (values=c(0.3,0.6))+scale_size_manual(values=c(2.5,4.5))+scale_colour_manual(values=c("#931751", "#005492","#aaaaaa"))+theme(aspect.ratio=1)
