INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

if (file.exists(file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.txt"))) {file.remove(file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.txt"))}

if (file.exists(file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.txt"))) {file.remove(file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.txt"))}


Chronos2DN_AN<-read.csv(file.path(INPUT_PATH, "Chronos2DN_AN.csv"))
Chronos2DN_OZ<-read.csv(file.path(INPUT_PATH, "Chronos2DN_OZ.csv"))
Chronos3DN_AN<-read.csv(file.path(INPUT_PATH, "Chronos3DN_AN.csv"))
Chronos3DN_OZ<-read.csv(file.path(INPUT_PATH, "Chronos3DN_OZ.csv"))
Chronos2DO_AN<-read.csv(file.path(INPUT_PATH, "Chronos2DO_AN.csv"))
Chronos2DO_OZ<-read.csv(file.path(INPUT_PATH, "Chronos2DO_OZ.csv"))
Chronos3DO_AN<-read.csv(file.path(INPUT_PATH, "Chronos3DO_AN.csv"))
Chronos3DO_OZ<-read.csv(file.path(INPUT_PATH, "Chronos3DO_OZ.csv"))

Chronos2DN<-cbind(Chronos2DN_AN,Chronos2DN_OZ)
Chronos3DN<-cbind(Chronos3DN_AN,Chronos3DN_OZ)
Chronos2DO<-cbind(Chronos2DO_AN,Chronos2DO_OZ)
Chronos3DO<-cbind(Chronos3DO_AN,Chronos3DO_OZ)

for(i in 4:ncol(Chronos2DO)) {
  mwut<-wilcox.test(Chronos2DO[,i],Chronos3DO[,i])
  mwup<-mwut$p.value
  ttest<-t.test(Chronos2DO[,i],Chronos3DO[,i])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  cat(i,tteste,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.txt"), append=TRUE)
}

for(i in 4:ncol(Chronos2DN)) {
  mwut<-wilcox.test(Chronos2DN[,i],Chronos3DN[,i])
  mwup<-mwut$p.value
  ttest<-t.test(Chronos2DN[,i],Chronos3DN[,i])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  cat(i,tteste,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.txt"), append=TRUE)
}
