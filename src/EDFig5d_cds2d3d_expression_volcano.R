INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

CNS2D3D_1<-read.csv(file.path(INPUT_PATH,"CNS2D3D_1.csv"))
CNS2D3D_2<-read.csv(file.path(INPUT_PATH,"CNS2D3D_2.csv"))
CNS2D3D_t<-cbind(CNS2D3D_1,CNS2D3D_2)

if (file.exists(file.path(INPUT_PATH, "CNS2D3D_t.txt"))) {file.remove(file.path(INPUT_PATH, "CNS2D3D_t.txt"))}


for (i in 3:ncol(CNS2D3D_t)) {
  mwut<-wilcox.test(CNS2D3D_t[,i]~CNS2D3D_t[,2])
  mwup<-mwut$p.value
  ttest<-t.test(CNS2D3D_t[,i]~CNS2D3D_t[,2])
  ttesp<-ttest$p.value
  ttese<-ttest$estimate
  cat(i,ttese,ttesp,mwup,"\n",file=file.path(OUTPUT_PATH, "CNS2D3D_t.txt"), append=TRUE)
}
