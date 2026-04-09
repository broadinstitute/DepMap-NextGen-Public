INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

GvsM1<-read.csv(file.path(INPUT_PATH,"GlialvsMesenchymal_1.csv"), check.names=FALSE)
GvsM2<-read.csv(file.path(INPUT_PATH,"GlialvsMesenchymal_2.csv"), check.names=FALSE)
GvsM<-cbind(GvsM1,GvsM2)
ess_annot_df <- read.csv(file.path(INPUT_PATH, 'CommonEssential_25q3p.csv'), col.names = c('symbol', 'CommonEssential'))

if (file.exists(file.path(OUTPUT_PATH, "GlialvsMesenchymal.csv"))) {file.remove(file.path(OUTPUT_PATH, "GlialvsMesenchymal.csv"))}

mean_glial <- c()
mean_mes <- c()
tps <- c()
mwups <- c()
for(i in 6:ncol(GvsM)) {
  ttest1<-t.test(GvsM[,i]~GvsM[,5]) 
  ttese1<-ttest1$estimate
  ttesp1<-ttest1$p.value
  MWUt1<-wilcox.test(GvsM[,i]~GvsM[,5])
  MWUp1<-MWUt1$p.value
  mean_glial <- c(mean_glial, ttese1['mean in group Glial'])
  mean_mes <- c(mean_mes, ttese1['mean in group Mesenchymal'])
  tps <- c(tps, ttesp1)
  mwups <- c(mwups, MWUp1)
}
ttest_df <- data.frame(list(mean_Glial=mean_glial,mean_Mesenchymal=mean_mes))
ttest_df['Glial-Mesenchymal'] <- ttest_df$mean_Glial - ttest_df$mean_Mesenchymal
ttest_df$t_p <- tps
ttest_df$MWU_p <- mwups
ttest_df$log10P <- -log10(ttest_df$MWU_p)
rownames(ttest_df) <- colnames(GvsM)[6:ncol(GvsM)]

ttest_df$CommonEssential <- rownames(ttest_df) %in% ess_annot_df$symbol
ttest_df$MWU_q <- p.adjust(ttest_df$MWU_p, method='BH')
ttest_df$class <- ifelse((ttest_df$`Glial-Mesenchymal` < -0.3) & (ttest_df$MWU_q < 0.05), 'q-significant', ifelse((ttest_df$`Glial-Mesenchymal` < -0.3) & (ttest_df$MWU_p < 0.05), 'p-significant', 'nonsignificant'))

ttest_df <- ttest_df[!is.na(ttest_df$t_p) & !ttest_df$CommonEssential, ]
ttest_df <- ttest_df[order(ttest_df$`Glial-Mesenchymal`), ]
write.csv(ttest_df, file.path(OUTPUT_PATH, "GlialvsMesenchymal.csv"))


