INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(dplyr)
library(viridis)

DMHan<-read.csv(file.path(INPUT_PATH,"DepMapHan.csv"))

if (file.exists(file.path(OUTPUT_PATH, "DMHan.csv"))) {file.remove(file.path(OUTPUT_PATH, "DMHan.csv"))}

DMHan.cor <- dplyr::select(DMHan,-Gene) %>%
  cor(use="pairwise.complete.obs")
DMHan.dist <- as.dist(1 - DMHan.cor)
DMHan.tree <- hclust(DMHan.dist, method="complete")
plot(DMHan.tree)
write.csv(DMHan.cor, file.path(OUTPUT_PATH,"DMHan.csv"))

library(gplots)
library(scales)

DMHan.dend <- as.dendrogram(DMHan.tree)
color.scheme <- alpha(viridis(50),0.7)
heatmap.2(DMHan.cor, 
          Rowv = DMHan.dend, 
          Colv = DMHan.dend, 
          dendrogram = "row", 
          trace = "none", density.info = "none", 
          breaks = seq(0.4, 1, length.out = 51),
          col = color.scheme,
          xlab = "NA", ylab = "NA")
