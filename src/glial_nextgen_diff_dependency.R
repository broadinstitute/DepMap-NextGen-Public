INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

GvsN1<-read.csv(file.path(INPUT_PATH, "GlialvsNextGen_1.csv"), check.names = FALSE)
GvsN2<-read.csv(file.path(INPUT_PATH, "GlialvsNextGen_2.csv"), check.names = FALSE)
GvsN<-cbind(GvsN1,GvsN2)
ess_annot_df <- read.csv(file.path(INPUT_PATH, 'CommonEssential_25q3p.csv'), col.names = c('symbol', 'CommonEssential'))

if (file.exists(file.path(OUTPUT_PATH, "GlialvsNextGen.csv"))) {file.remove(file.path(OUTPUT_PATH, "GlialvsNextGen.csv"))}

mean_glial <- c()
mean_ng <- c()
tps <- c()
mwups <- c()
for(i in 6:ncol(GvsN)) {
  ttest1<-t.test(GvsN[,i]~GvsN[,5]) 
  ttese1<-ttest1$estimate
  ttesp1<-ttest1$p.value
  MWUt1<-wilcox.test(GvsN[,i]~GvsN[,5])
  MWUp1<-MWUt1$p.value
  mean_glial <- c(mean_glial, ttese1['mean in group Glial'])
  mean_ng <- c(mean_ng, ttese1['mean in group non-glialNextGen'])
  tps <- c(tps, ttesp1)
  mwups <- c(mwups, MWUp1)
}
ttest_df <- data.frame(list(mean_Glial=mean_glial,mean_NextGen=mean_ng))
ttest_df['Glial-NextGen'] <- ttest_df$mean_Glial - ttest_df$mean_NextGen
ttest_df$t_p <- tps
ttest_df$MWU_p <- mwups
ttest_df$log10P <- -log10(ttest_df$MWU_p)
rownames(ttest_df) <- colnames(GvsN)[6:ncol(GvsN)]

ttest_df$CommonEssential <- rownames(ttest_df) %in% ess_annot_df$symbol
ttest_df$MWU_q <- p.adjust(ttest_df$MWU_p, method='BH')
ttest_df$Significance <- ifelse((ttest_df$`Glial-NextGen` < -0.3) & (ttest_df$MWU_q < 0.05), 'significant', 'nonsignificant')

ttest_df <- ttest_df[!is.na(ttest_df$t_p) & !ttest_df$CommonEssential, ]
ttest_df <- ttest_df[order(ttest_df$`Glial-NextGen`), ]
write.csv(ttest_df, file.path(OUTPUT_PATH, "GlialvsNextGen.csv"))

