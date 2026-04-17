INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

CDK6vsCN<-read.csv(file.path(INPUT_PATH,"CDK6vsCN.csv"))
chr_loc_df <- read.csv(file.path(INPUT_PATH,"ChrLocus.csv"), row.names=1)

if (file.exists(file.path(OUTPUT_PATH, "CDK6_CNassociation.csv"))) {file.remove(file.path(OUTPUT_PATH, "CDK6_CNassociation.csv"))}

corrs <- c()
corps <- c()
for (j in 7:ncol(CDK6vsCN)){
  corr1<-cor.test(CDK6vsCN[,j],CDK6vsCN[,6])
  corp1<-corr1$p.value
  core1<-corr1$estimate
  corrs <- c(corrs, core1)
  corps <- c(corps, corp1)
}
corr_df <- data.frame(list(Corr=corrs,Corr_p=corps))
rownames(corr_df) <- colnames(CDK6vsCN)[7:ncol(CDK6vsCN)]
corr_df$log10P <- -log10(corr_df$Corr_p)
corr_df$Corr_q <- p.adjust(corr_df$Corr_p, method='BH')
corr_df$Significance <- (corr_df$Corr_q < 0.05) & (abs(corr_df$Corr) > 0.3)
corr_df <- merge(corr_df, chr_loc_df, by='row.names')
corr_df <- data.frame(corr_df[, -1], row.names = corr_df[, "Row.names"])
corr_df <- corr_df[order(corr_df$Corr_q), ]
write.csv(corr_df, file.path(OUTPUT_PATH, "CDK6_CNassociation.csv"), quote = FALSE)

CDK6_CNassociation<-read.csv(file.path(OUTPUT_PATH,"CDK6_CNassociation.csv"))

ggplot(CDK6_CNassociation,aes(x=Corr,y=log10P))+geom_point(aes(size=Significance,alpha=Significance,colour=EnrichedLoc))+scale_size_manual(values=c(2,4))+scale_alpha_manual(values=c(0.2,0.4))+scale_colour_manual(values=c("#931751", "#005492","#aaaaaa"))+theme(aspect.ratio=1)
