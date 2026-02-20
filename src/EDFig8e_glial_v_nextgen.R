INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

if (file.exists(file.path(OUTPUT_PATH, "GlialvsMesenchymal.txt"))) {file.remove(file.path(OUTPUT_PATH, "GlialvsMesenchymal.txt"))}

GvsN1<-read.csv(file.path(INPUT_PATH, "GlialvsNextGen_1.csv"))
GvsN2<-read.csv(file.path(INPUT_PATH, "GlialvsNextGen_2.csv"))
GvsN<-cbind(GvsN1,GvsN2)

for(i in 7:ncol(GvsN)) {
  ttest1<-t.test(GvsN[,i]~GvsN[,6]) 
  ttese1<-ttest1$estimate
  ttesp1<-ttest1$p.value
  MWUt1<-wilcox.test(GvsN[,i]~GvsN[,6])
  MWUp1<-MWUt1$p.value
  cat(i,ttese1,ttesp1,MWUp1,"\n",file=file.path(OUTPUT_PATH, "GlialvsNextGen.txt"),append=TRUE)
}
