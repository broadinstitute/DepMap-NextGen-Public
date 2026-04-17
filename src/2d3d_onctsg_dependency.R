INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

mONCTSG_Gene<-read.csv(file.path(INPUT_PATH, "mONCTSG_Gene.csv")) 
mONCTSG_MutFree<-read.csv(file.path(INPUT_PATH, "mONCTSG_MutFree.csv")) 
mONCTSG_Pathway<-read.csv(file.path(INPUT_PATH, "mONCTSG_Pathway.csv")) 
onc_tsg_gene_annot <- read.csv(file.path(INPUT_PATH, "OncTSG_annotations.csv"), row.names = 1)

if (file.exists(file.path(OUTPUT_PATH, "mONCTSG_Gene_Result2.csv"))) {file.remove(file.path(OUTPUT_PATH, "mONCTSG_Gene_Result2.csv"))}
if (file.exists(file.path(OUTPUT_PATH, "mONCTSG_Pathway_Result.csv"))) {file.remove(file.path(OUTPUT_PATH, "mONCTSG_Pathway_Result.csv"))}

mean_2ds <- c()
mean_3do <- c()
tps <- c()
mwups <- c()
for(i in 5:ncol(mONCTSG_Gene)) {
  mwut<-wilcox.test(mONCTSG_Gene[,i]~mONCTSG_Gene[,4])
  mwup<-mwut$p.value
  ttest<-t.test(mONCTSG_Gene[,i]~mONCTSG_Gene[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  mean_2ds <- c(mean_2ds, tteste['mean in group 2DS'])
  mean_3do <- c(mean_3do, tteste['mean in group 3DO'])
  tps <- c(tps, ttesp)
  mwups <- c(mwups, mwup)
}
ttest_df <- data.frame(list(mean_2DS=mean_2ds,mean_3DO=mean_3do))
ttest_df$t_p <- tps
ttest_df$MWU_p <- mwups
ttest_df$log10P <- -log10(ttest_df$MWU_p)
ttest_df$MWU_q <- p.adjust(ttest_df$MWU_p, method='BH')
rownames(ttest_df) <- colnames(mONCTSG_Gene)[5:ncol(mONCTSG_Gene)]
ttest_df <- ttest_df[ttest_df$MWU_q < 0.05, ]
write.csv(ttest_df, file.path(OUTPUT_PATH, "mONCTSG_Gene_Result2.csv"), quote=FALSE)

mutfree_p <- c()
for(i in 5:ncol(mONCTSG_MutFree)) {
  mwut<-wilcox.test(mONCTSG_MutFree[,i]~mONCTSG_MutFree[,4])
  mwup<-mwut$p.value
  mutfree_p <- c(mutfree_p, mwup)
}
mutfree_df <- data.frame(list(mutfree_p=mutfree_p, mutfree_q=p.adjust(mutfree_p, method='BH')))
mutfree_df$mf_confirmation <- mutfree_df$mutfree_q < 0.05
rownames(mutfree_df) <- colnames(mONCTSG_MutFree)[5:ncol(mONCTSG_MutFree)]

ttest_df <- merge(ttest_df, mutfree_df, by='row.names')
ttest_df <- data.frame(ttest_df[, -1], row.names = ttest_df[, "Row.names"])
ttest_df <- merge(ttest_df, onc_tsg_gene_annot, by='row.names')
ttest_df <- data.frame(ttest_df[, -1], row.names = ttest_df[, "Row.names"])
ttest_df$significance <- ttest_df$MWU_q < 0.05
ttest_df$class <- ifelse(ttest_df$significance, 'significant', 'insignificant')
ttest_df <- ttest_df[c('oncTSG', 'mean_2DS', 'mean_3DO', 't_p', 'MWU_p', 'log10P', 'MWU_q', 'mutfree_q', 'significance', 'mf_confirmation', 'class')]
ttest_df <- ttest_df[order(ttest_df$MWU_q), ]
write.csv(ttest_df, file.path(OUTPUT_PATH, "mONCTSG_Gene_Result2.csv"), quote=FALSE)


mean_2ds <- c()
mean_3do <- c()
tps <- c()
mwups <- c()
for(i in 5:ncol(mONCTSG_Pathway)) {
  mwut<-wilcox.test(mONCTSG_Pathway[,i]~mONCTSG_Pathway[,4])
  mwup<-mwut$p.value
  ttest<-t.test(mONCTSG_Pathway[,i]~mONCTSG_Pathway[,4])
  tteste<-ttest$estimate
  ttesp<-ttest$p.value
  mean_2ds <- c(mean_2ds, tteste['mean in group 2DS'])
  mean_3do <- c(mean_3do, tteste['mean in group 3DO'])
  tps <- c(tps, ttesp)
  mwups <- c(mwups, mwup)
}
ttest_df <- data.frame(list(mean_2DS=mean_2ds,mean_3DO=mean_3do))
ttest_df$t_p <- tps
ttest_df$MWU_p <- mwups
ttest_df$log10P <- -log10(ttest_df$MWU_p)
ttest_df$MWU_q <- p.adjust(ttest_df$MWU_p, method='BH')
ttest_df$significance <- ttest_df$MWU_q < 0.05
rownames(ttest_df) <- colnames(mONCTSG_Pathway)[5:ncol(mONCTSG_Pathway)]
ttest_df <- ttest_df[-which(rownames(ttest_df) == 'Grand.Average'), ]
ttest_df$oncTSG <- ifelse(grepl('ONC', rownames(ttest_df)), 'ONC', 'TSG')
rownames(ttest_df) <- gsub('TGF\\.beta', 'TGF-beta', rownames(ttest_df))
rownames(ttest_df) <- gsub('RTK\\.RAS', 'RTK/RAS', rownames(ttest_df))
rownames(ttest_df) <- gsub('\\.', ' ', rownames(ttest_df))
ttest_df <- ttest_df[c('oncTSG', 'mean_2DS', 'mean_3DO', 't_p', 'MWU_p', 'log10P', 'MWU_q', 'significance')]
ttest_df <- ttest_df[order(ttest_df$MWU_q), ]
write.csv(ttest_df, file.path(OUTPUT_PATH, "mONCTSG_Pathway_Result.csv"), quote=FALSE)


mOT_gene2<-read.csv(file.path(OUTPUT_PATH, "mONCTSG_Gene_Result2.csv"))
mOT_pathway<-read.csv(file.path(OUTPUT_PATH, "mONCTSG_Pathway_Result.csv"))

ggplot(mOT_gene2,aes(x=mean_2DS,y=mean_3DO))+geom_point(shape=21, alpha=0.5, aes(fill=oncTSG,colour=oncTSG,size=log10P))+scale_color_manual(values=c("#931751","#005492"))+scale_fill_manual(values=c("#931751","#005492"))+scale_alpha_manual(values=c(0.6,0.3,0.3))+scale_size(range=c(4,16))+theme(aspect.ratio=1)+xlim(-2.2,0.4)+ylim(-2.2,0.4)

ggplot(mOT_pathway,aes(x=mean_2DS,y=mean_3DO))+geom_point(alpha=0.6, aes(shape=significance, colour=oncTSG,fill=oncTSG,size=log10P))+scale_shape_manual(values=c(1,21))+scale_color_manual(values=c("#931751","#005492"))+scale_fill_manual(values=c("#931751","#005492"))+scale_size(range=c(3,30))+theme(aspect.ratio=1)+xlim(-0.6,0.2)+ylim(-0.6,0.2)

