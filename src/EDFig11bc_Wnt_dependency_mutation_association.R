INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

PDACc_Chronos_wnt3a<- read.csv(file.path(INPUT_PATH, "PDACc_Chronos_wnt3a.csv"))
if (file.exists(file.path(INPUT_PATH, "PDACc_associated_dep_wnt3a.txt"))) {file.remove(file.path(INPUT_PATH, "PDACc_associated_dep_wnt3a.txt"))}
if (file.exists(file.path(INPUT_PATH, "PDACc_associated_dep_wnt3a_q.csv"))) {file.remove(file.path(INPUT_PATH, "PDACc_associated_dep_wnt3a_q.csv"))}
if (file.exists(file.path(INPUT_PATH, "WntTop5_MutCorr.txt"))) {file.remove(file.path(INPUT_PATH, "WntTop5_MutCorr.txt"))}



for (i in 8:ncol(PDACc_Chronos_wnt3a)) {
    corr2<-cor.test(PDACc_Chronos_wnt3a[,i],PDACc_Chronos_wnt3a[,7])
    corp2<-corr2$p.value
    corc2<-corr2$estimate
    cat(i,corc2,corp2,"\n",file=file.path(INPUT_PATH, "PDACc_associated_dep_wnt3a.txt"), append=TRUE)
}

PDACc_associated_dep_wnt3a<-read.csv(file.path(INPUT_PATH, "PDACc_associated_dep_wnt3a.csv"))

PDACc_associated_dep_wnt3a_q<-p.adjust(PDACc_associated_dep_wnt3a$Corr_p,method="BH")
write.csv(PDACc_associated_dep_wnt3a_q,file.path(INPUT_PATH, "PDACc_associated_dep_wnt3a_q.csv"))

RepMut_Org<-read.csv(file.path(INPUT_PATH, "RepMut_Org.csv"))
WntTop5_Org<-read.csv(file.path(INPUT_PATH, "WntTop5_Org.csv"))


for(i in 4:8) {
  for (j in 4:ncol(RepMut_Org)){
    mwut<-wilcox.test(WntTop5_Org[,i]~RepMut_Org[,j])
    mwup<-mwut$p.value
    ttest<-t.test(WntTop5_Org[,i]~RepMut_Org[,j])
    ttese<-ttest$estimate
    ttesp<-ttest$p.value
    cat(i,j,ttese,ttesp,mwup,"\n",file=file.path(INPUT_PATH, "WntTop5_MutCorr.txt"),append=TRUE)
  }
}

WntTop5_MutCorr<-read.csv(file.path(INPUT_PATH, "WntTop5_MutCorr.csv"))
ggplot(WntTop5_MutCorr,aes(x=DeltaChronos,y=log10P,colour=class))+geom_point(alpha=0.6, size=3)+theme(aspect.ratio=1)+scale_color_manual(values=c("#bbbbbb","#931751","#005492"))

