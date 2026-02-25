INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

GvsN1<-read.csv(file.path(INPUT_PATH, "GlialvsNextGen_1.csv"))
GvsN2<-read.csv(file.path(INPUT_PATH, "GlialvsNextGen_2.csv"))
GvsN<-cbind(GvsN1,GvsN2)

if (file.exists(file.path(OUTPUT_PATH, "GlialvsNextGen.txt"))) {file.remove(file.path(OUTPUT_PATH, "GlialvsNextGen.txt"))}

for(i in 6:ncol(GvsN)) {
  ttest1<-t.test(GvsN[,i]~GvsN[,5]) 
  ttese1<-ttest1$estimate
  ttesp1<-ttest1$p.value
  MWUt1<-wilcox.test(GvsN[,i]~GvsN[,5])
  MWUp1<-MWUt1$p.value
  cat(i,ttese1,ttesp1,MWUp1,"\n",file=file.path(OUTPUT_PATH, "GlialvsNextGen.txt"),append=TRUE)
}

