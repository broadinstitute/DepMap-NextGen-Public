INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

RepMut_Org<-read.csv(file.path(INPUT_PATH, "RepMut_Org.csv"))
WntTop5_Org<-read.csv(file.path(INPUT_PATH, "WntTop5_Org.csv"))

if (file.exists(file.path(OUTPUT_PATH, "WntTop5_MutCorr.txt"))) {file.remove(file.path(OUTPUT_PATH, "WntTop5_MutCorr.txt"))}

for(i in 6:10) {
  for (j in 4:ncol(RepMut_Org)){
    mwut<-wilcox.test(WntTop5_Org[,i]~RepMut_Org[,j])
    mwup<-mwut$p.value
    ttest<-t.test(WntTop5_Org[,i]~RepMut_Org[,j])
    ttese<-ttest$estimate
    ttesp<-ttest$p.value
    cat(i,j,ttese,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "WntTop5_MutCorr.txt"),append=TRUE)
  }
}

WntTop5_MutCorr<-read.csv(file.path(INPUT_PATH, "WntTop5_MutCorr.csv"))

ggplot(WntTop5_MutCorr,aes(x=MUT.WT,y=log10P,colour=class))+geom_point(alpha=0.6, size=3)+theme(aspect.ratio=1)+scale_color_manual(values=c("#bbbbbb","#931751","#005492"))
