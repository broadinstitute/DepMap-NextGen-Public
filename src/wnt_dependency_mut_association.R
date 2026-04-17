INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)
library(ggplot2)

RepMut_Org<-read.csv(file.path(INPUT_PATH, "RepMut_Org.csv"))
WntTop5_Org<-read.csv(file.path(INPUT_PATH, "WntTop5_Org.csv"))
oncokb_annot_df <- read.csv(file.path(INPUT_PATH, "OncoKB_mutant_class.csv"))[c('Mut', 'oncoKB', 'class')]

if (file.exists(file.path(OUTPUT_PATH, "WntTop5_MutCorr.csv"))) {file.remove(file.path(OUTPUT_PATH, "WntTop5_MutCorr.csv"))}

wnt_genes <- c()
mut_strings <- c()
mut_genes <- c()
mut_types <- c()
mean_wts <- c()
mean_muts <- c()
tps <- c()
mwups <- c()
for(i in 6:10) {
  for (j in 4:ncol(RepMut_Org)){
    wnt_genes <- c(wnt_genes, colnames(WntTop5_Org)[i])
    mut_string <- colnames(RepMut_Org)[j]
    mut_string <- sub('HLA\\.A', 'HLA-A', mut_string)
    mut_string <- sub('HLA\\.B', 'HLA-B', mut_string)
    mut_strings <- c(mut_strings, mut_string)
    mut_genes <- c(mut_genes, strsplit(mut_string, "_")[[1]][1]) 
    mut_types <- c(mut_types, strsplit(mut_string, "_")[[1]][2]) 
    
    mwut<-wilcox.test(WntTop5_Org[,i]~RepMut_Org[,j])
    mwup<-mwut$p.value
    ttest<-t.test(WntTop5_Org[,i]~RepMut_Org[,j])
    ttese<-ttest$estimate
    ttesp<-ttest$p.value
    mean_wts <- c(mean_wts, ttese['mean in group FALSE'])
    mean_muts <- c(mean_muts, ttese['mean in group TRUE'])
    tps <- c(tps, ttesp)
    mwups <- c(mwups, mwup)
  }
}
ttest_df <- data.frame(list(WNT=wnt_genes,MutGene=mut_genes,MutType=mut_types,Mut=mut_strings,
                            mean_WT=mean_wts,mean_MUT=mean_muts))
ttest_df['MUT-WT'] <- ttest_df$mean_MUT - ttest_df$mean_WT
ttest_df$t_p <- tps
ttest_df$MWU_p <- mwups
ttest_df$log10P <- -log10(ttest_df$MWU_p)
ttest_df$MWU_q <- p.adjust(ttest_df$MWU_p, method='BH')
ttest_df <- merge(ttest_df, oncokb_annot_df, by='Mut')
ttest_df <- ttest_df[c('WNT', 'MutGene', 'MutType', 'Mut', 'mean_WT', 'mean_MUT', 'MUT-WT', 't_p', 'MWU_p', 'log10P', 'MWU_q', 'oncoKB', 'class')]
ttest_df <- ttest_df[order(ttest_df$MWU_p, ttest_df$class), ]
write.csv(ttest_df, file.path(OUTPUT_PATH, "WntTop5_MutCorr.csv"), quote=FALSE, row.names=FALSE)

WntTop5_MutCorr<-read.csv(file.path(OUTPUT_PATH, "WntTop5_MutCorr.csv"))

ggplot(WntTop5_MutCorr,aes(x=MUT.WT,y=log10P,colour=class))+geom_point(alpha=0.6, size=3)+theme(aspect.ratio=1)+scale_color_manual(values=c("#bbbbbb","#931751","#005492"))
