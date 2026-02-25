INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

ChronosCNS_AN<-read.csv(file.path(INPUT_PATH, "ChronosCNS_AN.csv"))
ChronosCNS_OZ<-read.csv(file.path(INPUT_PATH, "ChronosCNS_OZ.csv"))
ChronosOrg_AN<-read.csv(file.path(INPUT_PATH, "ChronosOrg_AN.csv"))
ChronosOrg_OZ<-read.csv(file.path(INPUT_PATH, "ChronosOrg_OZ.csv"))

ChronosCNS<-cbind(ChronosCNS_AN, ChronosCNS_OZ)
ChronosOrg<-cbind(ChronosOrg_AN, ChronosOrg_OZ)

if (file.exists(file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.txt"))) {file.remove(file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.txt"))}
if (file.exists(file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.txt"))) {file.remove(file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.txt"))}
if (file.exists(file.path(OUTPUT_PATH, "CNS2Dvs3D_q.csv"))) {file.remove(file.path(OUTPUT_PATH, "CNS2Dvs3D_q.csv"))}
if (file.exists(file.path(OUTPUT_PATH, "Org2Dvs3D_q.csv"))) {file.remove(file.path(OUTPUT_PATH, "Org2Dvs3D_q.csv"))}

for(i in 5:ncol(ChronosCNS)) {
  mwut<-wilcox.test(ChronosCNS[,i]~ChronosCNS[,4])
  mwup<-mwut$p.value
  ttest<-t.test(ChronosCNS[,i]~ChronosCNS[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  cat(i,tteste,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.txt"), append=TRUE)
}

for(i in 5:ncol(ChronosOrg)) {
  mwut<-wilcox.test(ChronosOrg[,i]~ChronosOrg[,4])
  mwup<-mwut$p.value
  ttest<-t.test(ChronosOrg[,i]~ChronosOrg[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  cat(i,tteste,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.txt"), append=TRUE)
}


CNSOrg_2Dvs3D_p<-read.csv(file.path(INPUT_PATH, "CNSOrg_2Dvs3D_p.csv"))

CNS2Dvs3D_q<-p.adjust(CNSOrg_2Dvs3D_p$CNS,method="BH")
Org2Dvs3D_q<-p.adjust(CNSOrg_2Dvs3D_p$Org,method="BH")
write.csv(CNS2Dvs3D_q, file.path(OUTPUT_PATH, "CNS2Dvs3D_q.csv"))
write.csv(Org2Dvs3D_q, file.path(OUTPUT_PATH, "Org2Dvs3D_q.csv"))
