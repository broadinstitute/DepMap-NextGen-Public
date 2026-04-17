INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

ChronosCNS_AN<-read.csv(file.path(INPUT_PATH, "ChronosCNS_AN.csv"), check.names = FALSE)
ChronosCNS_OZ<-read.csv(file.path(INPUT_PATH, "ChronosCNS_OZ.csv"), check.names = FALSE)
ChronosOrg_AN<-read.csv(file.path(INPUT_PATH, "ChronosOrg_AN.csv"), check.names = FALSE)
ChronosOrg_OZ<-read.csv(file.path(INPUT_PATH, "ChronosOrg_OZ.csv"), check.names = FALSE)

ChronosCNS<-cbind(ChronosCNS_AN, ChronosCNS_OZ)
ChronosOrg<-cbind(ChronosOrg_AN, ChronosOrg_OZ)

if (file.exists(file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.csv"))) {file.remove(file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.csv"))}
if (file.exists(file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.csv"))) {file.remove(file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.csv"))}

mean_2dn <- c()
mean_3dn <- c()
tps <- c()
mwups <- c()
for(i in 5:ncol(ChronosCNS)) {
  mwut<-wilcox.test(ChronosCNS[,i]~ChronosCNS[,4])
  mwup<-mwut$p.value
  ttest<-t.test(ChronosCNS[,i]~ChronosCNS[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  mean_2dn <- c(mean_2dn, tteste['mean in group 2DN'])
  mean_3dn <- c(mean_3dn, tteste['mean in group 3DN'])
  tps <- c(tps, ttesp)
  mwups <- c(mwups, mwup)
}
ttest_df <- data.frame(list(mean_2DN=mean_2dn,mean_3DN=mean_3dn))
ttest_df['3DN-2DN'] <- ttest_df$mean_3DN - ttest_df$mean_2DN
ttest_df$t_p <- tps
ttest_df$MWU_p <- mwups
ttest_df$log10P <- -log10(ttest_df$MWU_p)
ttest_df$MWU_q <- p.adjust(ttest_df$MWU_p, method='BH')
rownames(ttest_df) <- colnames(ChronosCNS)[5:ncol(ChronosCNS)]
ttest_df <- ttest_df[!is.na(ttest_df$t_p), ]
write.csv(ttest_df, file.path(OUTPUT_PATH, "Chronos 2DN vs 3DN.csv"))


mean_2do <- c()
mean_3do <- c()
tps <- c()
mwups <- c()
for(i in 5:ncol(ChronosOrg)) {
  mwut<-wilcox.test(ChronosOrg[,i]~ChronosOrg[,4])
  mwup<-mwut$p.value
  ttest<-t.test(ChronosOrg[,i]~ChronosOrg[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  mean_2do <- c(mean_2do, tteste['mean in group 2DO'])
  mean_3do <- c(mean_3do, tteste['mean in group 3DO'])
  tps <- c(tps, ttesp)
  mwups <- c(mwups, mwup)
}
ttest_df <- data.frame(list(mean_2DO=mean_2do,mean_3DO=mean_3do))
ttest_df['3DO-2DO'] <- ttest_df$mean_3DO - ttest_df$mean_2DO
ttest_df$t_p <- tps
ttest_df$MWU_p <- mwups
ttest_df$log10P <- -log10(ttest_df$MWU_p)
ttest_df$MWU_q <- p.adjust(ttest_df$MWU_p, method='BH')
rownames(ttest_df) <- colnames(ChronosOrg)[5:ncol(ChronosOrg)]
ttest_df <- ttest_df[!is.na(ttest_df$t_p), ]
write.csv(ttest_df, file.path(OUTPUT_PATH, "Chronos 2DO vs 3DO.csv"))
