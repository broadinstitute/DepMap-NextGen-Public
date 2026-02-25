INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

CDK6_MutAssoc<-read.csv(file.path(INPUT_PATH, "CDK6_GBM_MutAssociation.csv"))

if (file.exists(file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result.txt"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result.txt"))}
if (file.exists(file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result_q.csv"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result_q.csv"))}

for (j in 6:ncol(CDK6_MutAssoc)){
  ttest<-t.test(CDK6_MutAssoc[,5]~CDK6_MutAssoc[,j])
  ttese<-ttest$estimate
  ttesp<-ttest$p.value
  MWUt<-wilcox.test(CDK6_MutAssoc[,5]~CDK6_MutAssoc[,j])
  MWUp<-MWUt$p.value
  cat(j,ttese,ttesp,MWUp,"\n",file=file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result.txt"),append=TRUE)
}


CDK6_MutAssoc_Result<-read.csv(file.path(INPUT_PATH, "CDK6_MutAssoc_Result.csv"))

CDK6_MutAssoc_Result_q<-p.adjust(CDK6_MutAssoc_Result$MWU_p,method="BH")
write.csv(CDK6_MutAssoc_Result_q, file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result_q.csv"))
