INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(stats)

PDACc_Chronos_wnt3a<-read.csv(file.path(INPUT_PATH, "PDACc_Chronos_wnt3a.csv"), check.names = FALSE)
all_genesets <- read.csv(file.path(INPUT_PATH, "geneset_table_updated.csv"))
wnt3a_genes <- all_genesets[all_genesets$Geneset == "GOBP_WNT_SIGNALING_PATHWAY", "GeneSymbol"]

if (file.exists(file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a.csv"))) {file.remove(file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a.csv"))}

corrs <- c()
corps <- c()
for (i in 8:ncol(PDACc_Chronos_wnt3a)) {
  corr2<-cor.test(PDACc_Chronos_wnt3a[,i],PDACc_Chronos_wnt3a[,7])
  corp2<-corr2$p.value
  corc2<-corr2$estimate
  corrs <- c(corrs, corc2)
  corps <- c(corps, corp2)
}
corr_df <- data.frame(list(Corr=corrs,Corr_p=corps))
rownames(corr_df) <- colnames(PDACc_Chronos_wnt3a)[8:ncol(PDACc_Chronos_wnt3a)]
corr_df$log10P <- -log10(corr_df$Corr_p)
corr_df$MsigDB_WNT <- rownames(corr_df) %in% wnt3a_genes
corr_df$Corr_q <- p.adjust(corr_df$Corr_p, method='BH')
corr_df <- corr_df[order(corr_df$Corr), ]
write.csv(corr_df, file.path(OUTPUT_PATH, "PDACc_associated_dep_wnt3a.csv"))
