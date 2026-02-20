INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

if (file.exists(file.path(OUTPUT_PATH, "GlialvsMesenchymal.txt"))) {file.remove(file.path(OUTPUT_PATH, "GlialvsMesenchymal.txt"))}

GvsM1<- read.csv(file.path(INPUT_PATH, "GlialvsMesenchymal_1.csv"))
GvsM2<- read.csv(file.path(INPUT_PATH, "GlialvsMesenchymal_2.csv"))
GvsM<- cbind(GvsM1,GvsM2)

for(i in 7:ncol(GvsM)) {
  ttest1<-t.test(GvsM[,i]~GvsM[,6]) 
  ttese1<-ttest1$estimate
  ttesp1<-ttest1$p.value
  MWUt1<-wilcox.test(GvsM[,i]~GvsM[,6])
  MWUp1<-MWUt1$p.value

  cat(i,ttese1,ttesp1,MWUp1,"\n",file=file.path(OUTPUT_PATH, "GlialvsMesenchymal.txt"),append=TRUE)
}

