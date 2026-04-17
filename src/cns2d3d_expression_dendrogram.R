INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(dplyr)

CNS2D3D<-read.csv(file.path(INPUT_PATH,"CNS_NextGen2D3D_HV.csv"))

if (file.exists(file.path(OUTPUT_PATH, "CNS2D3Dcor.csv"))) {file.remove(file.path(OUTPUT_PATH, "CNS2D3Dcor.csv"))}
if (file.exists(file.path(OUTPUT_PATH, "CNS2D3Dcor_list.csv"))) {file.remove(file.path(OUTPUT_PATH, "CNS2D3Dcor_list.csv"))}

CNS2D3D.cor <- dplyr::select(CNS2D3D, -Gene, -Variance) %>%
  cor(use="pairwise.complete.obs")
CNS2D3D.dist <- as.dist(1 - CNS2D3D.cor)
CNS2D3D.tree <- hclust(CNS2D3D.dist, method="complete")
plot(CNS2D3D.tree)
write.csv(CNS2D3D.cor,file.path(OUTPUT_PATH,"CNS2D3Dcor.csv"))

adherent_samples <- rownames(CNS2D3D.cor)[rownames(CNS2D3D.cor) %>% grepl('\\.2D', .)]
nextgen_samples <- rownames(CNS2D3D.cor)[rownames(CNS2D3D.cor) %>% grepl('\\.3D', .)]
adh_adh <- CNS2D3D.cor[adherent_samples, adherent_samples]
adh_adh_values <- adh_adh[lower.tri(adh_adh)]
adh_adh_df <- data.frame(list(PairClass=rep('Adh:Adh', times=length(adh_adh_values)), Corr=adh_adh_values))
adh_adh_df <- adh_adh_df[order(adh_adh_df$Corr), ]

sph_sph <- CNS2D3D.cor[nextgen_samples, nextgen_samples]
sph_sph_values <- sph_sph[lower.tri(sph_sph)]
sph_sph_df <- data.frame(list(PairClass=rep('Sph:Sph', times=length(sph_sph_values)), Corr=sph_sph_values))
sph_sph_df <- sph_sph_df[order(sph_sph_df$Corr), ]

adh_sph_values <- CNS2D3D.cor[adherent_samples, nextgen_samples]
cross_vals <- c()
same_vals <- c()
for (i in 1:nrow(adh_sph_values)){
  for (j in 1:ncol(adh_sph_values)){
    if (strsplit(adherent_samples[i], '\\.')[[1]][1] == strsplit(nextgen_samples[j], '\\.')[[1]][1]) {
      same_vals <- c(same_vals, adh_sph_values[i, j])
    } else {
      cross_vals <- c(cross_vals, adh_sph_values[i, j])
    }
  }
}
cross_df <- data.frame(list(PairClass=rep('Adh:Sph', times=length(cross_vals)), Corr=cross_vals))
cross_df <- cross_df[order(cross_df$Corr), ]
same_df <- data.frame(list(PairClass=rep('Same model', times=length(same_vals)), Corr=same_vals))
same_df <- same_df[order(same_df$Corr), ]
longform_df <- do.call(rbind, list(adh_adh_df, sph_sph_df, cross_df, same_df))
write.csv(longform_df, file.path(OUTPUT_PATH,"CNS2D3Dcor_list.csv"), row.names = FALSE)

library(gplots)
library(scales)
library(RColorBrewer)

CNS2D3D.dend <- as.dendrogram(CNS2D3D.tree)
color.scheme <- alpha(brewer.pal(9,"YlOrBr"),0.7)
heatmap.2(CNS2D3D.cor, 
          Rowv = CNS2D3D.dend, 
          Colv = CNS2D3D.dend, 
          dendrogram = "row", 
          trace = "none", density.info = "none", 
          col = color.scheme,
          xlab = "NA", ylab = "NA")
