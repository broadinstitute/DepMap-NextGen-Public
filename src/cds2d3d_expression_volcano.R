INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

CNS2D3D_1<-read.csv(file.path(INPUT_PATH,"CNS2D3D_1.csv"), check.names = FALSE)
CNS2D3D_2<-read.csv(file.path(INPUT_PATH,"CNS2D3D_2.csv"), check.names = FALSE)
CNS2D3D_t<-cbind(CNS2D3D_1,CNS2D3D_2)

if (file.exists(file.path(INPUT_PATH, "CNS2D3D_t.csv"))) {file.remove(file.path(INPUT_PATH, "CNS2D3D_t.csv"))}

mean_2d <- c()
mean_org <- c()
tps <- c()
mwups <- c()
for (i in 3:ncol(CNS2D3D_t)) {
  mwut<-wilcox.test(CNS2D3D_t[,i]~CNS2D3D_t[,2])
  mwup<-mwut$p.value
  ttest<-t.test(CNS2D3D_t[,i]~CNS2D3D_t[,2])
  ttesp<-ttest$p.value
  ttese<-ttest$estimate
  mean_2d <- c(mean_2d, ttese['mean in group 2D'])
  mean_org <- c(mean_org, ttese['mean in group 3D'])
  tps <- c(tps, ttesp)
  mwups <- c(mwups, mwup)
}
ttest_df <- data.frame(list(`2D_mean`=mean_2d,Org_mean=mean_org), check.names = FALSE)
ttest_df['2D-Org'] <- ttest_df$`2D_mean` - ttest_df$Org_mean
ttest_df$t_p <- tps
ttest_df$MWUp <- mwups
ttest_df$log10P <- -log10(ttest_df$MWUp)
ttest_df$MWUq <- p.adjust(ttest_df$MWUp, method='BH')
rownames(ttest_df) <- unlist(lapply(colnames(CNS2D3D_t)[3:ncol(CNS2D3D_t)], function(x){strsplit(x, ' ')[[1]][1]}))
ttest_df <- ttest_df[!is.na(ttest_df$t_p), ]
ttest_df <- ttest_df[order(ttest_df$MWUp), ]

write.csv(ttest_df, file.path(OUTPUT_PATH, "CNS2D3D_t.csv"), quote=FALSE)

