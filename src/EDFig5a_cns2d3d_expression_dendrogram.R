INPUT_PATH = file.path("..", "data")
OUTPUT_PATH = file.path("..", "processed")

library(dplyr)

CNS2D3D<-read.csv(file.path(INPUT_PATH,"CNS_NextGen2D3D_HV.csv"))

if (file.exists(file.path(OUTPUT_PATH, "CNS2D3Dcor.csv"))) {file.remove(file.path(OUTPUT_PATH, "CNS2D3Dcor.csv"))}

CNS2D3D.cor <- dplyr::select(CNS2D3D, -Gene, -Variance) %>%
  cor(use="pairwise.complete.obs")
CNS2D3D.dist <- as.dist(1 - CNS2D3D.cor)
CNS2D3D.tree <- hclust(CNS2D3D.dist, method="complete")
plot(CNS2D3D.tree)
write.csv(CNS2D3D.cor,file.path(OUTPUT_PATH,"CNS2D3Dcor.csv"))

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
