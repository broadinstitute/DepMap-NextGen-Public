INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

CDK6vsCN<-read.csv(file.path(INPUT_PATH,"CDK6vsCN.csv"))

if (file.exists(file.path(OUTPUT_PATH, "CDK6_CNassociation.txt"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_CNassociation.txt"))}

for (j in 7:ncol(CDK6vsCN)){
  corr1<-cor.test(CDK6vsCN[,j],CDK6vsCN[,6])
  corp1<-corr1$p.value
  core1<-corr1$estimate
  cat(j,corp1,core1,"\n",file=file.path(OUTPUT_PATH, "CDK6_CNassociation.txt"), append=TRUE)
}

CDK6_CNassociation<-read.csv(file.path(INPUT_PATH,"CDK6_CNassociation.csv"))

ggplot(CDK6_CNassociation,aes(x=Corr,y=log10P))+geom_point(aes(size=Significance,alpha=Significance,colour=EnrichedLoc))+scale_size_manual(values=c(2,4))+scale_alpha_manual(values=c(0.2,0.4))+scale_colour_manual(values=c("#931751", "#005492","#aaaaaa"))+theme(aspect.ratio=1)
