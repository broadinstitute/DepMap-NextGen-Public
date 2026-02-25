INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

PDACc_Chronos_wnt3a<-read.csv(file.path(INPUT_PATH, "PDACc_Chronos_wnt3a.csv"))

if (file.exists(file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a.txt"))) {file.remove(file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a.txt"))}
if (file.exists(file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a_q.csv"))) {file.remove(file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a_q.csv"))}

for (i in 8:ncol(PDACc_Chronos_wnt3a)) {
  corr2<-cor.test(PDACc_Chronos_wnt3a[,i],PDACc_Chronos_wnt3a[,7])
  corp2<-corr2$p.value
  corc2<-corr2$estimate
  cat(i,corc2,corp2,"\n",file=file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a.txt"), append=TRUE)
}

PDACc_associated_dep_wnt3a<-read.csv(file.path(INPUT_PATH, "PDACc_associated_dep_wnt3a.csv"))

PDACc_associated_dep_wnt3a_q<-p.adjust(PDACc_associated_dep_wnt3a$Corr_p,method="BH")
write.csv(PDACc_associated_dep_wnt3a_q, file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a_q.csv"))
