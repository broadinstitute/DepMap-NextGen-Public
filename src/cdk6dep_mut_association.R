INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

CDK6_MutAssoc<-read.csv(file.path(INPUT_PATH, "CDK6_GBM_MutAssociation.csv"))

if (file.exists(file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result.csv"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result.csv"))}

mean_wts <- c()
mean_muts <- c()
tps <- c()
mwups <- c()
for (j in 6:ncol(CDK6_MutAssoc)){
  ttest<-t.test(CDK6_MutAssoc[,5]~CDK6_MutAssoc[,j])
  ttese<-ttest$estimate
  ttesp<-ttest$p.value
  MWUt<-wilcox.test(CDK6_MutAssoc[,5]~CDK6_MutAssoc[,j])
  MWUp<-MWUt$p.value
  mean_wts <- c(mean_wts, ttese['mean in group FALSE'])
  mean_muts <- c(mean_muts, ttese['mean in group TRUE'])
  tps <- c(tps, ttesp)
  mwups <- c(mwups, MWUp)
}
ttest_df <- data.frame(list(Mean_WT=mean_wts,Mean_Mut=mean_muts))
ttest_df['Mut-WT'] <- ttest_df$Mean_Mut - ttest_df$Mean_WT
ttest_df$t_p <- tps
ttest_df$MWU_p <- mwups
ttest_df$log10P <- -log10(ttest_df$MWU_p)
ttest_df$MWU_q <- p.adjust(ttest_df$MWU_p, method='BH')
rownames(ttest_df) <- colnames(CDK6_MutAssoc)[6:ncol(CDK6_MutAssoc)]
write.csv(ttest_df, file.path(OUTPUT_PATH, "CDK6_MutAssoc_Result.csv"), quote=FALSE)
