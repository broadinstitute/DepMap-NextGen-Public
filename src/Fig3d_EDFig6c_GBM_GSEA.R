### Gene Set Enrichment Analysis (GSEA) for GBM clusters
# Purpose: This script performs GSEA using metaprograms and generates enrichment plots.

# Define input, output, and figure directories
INPUT_PATH <- file.path("..", "data")
OUTPUT_PATH <- file.path("..", "processed")
FIGURE_PATH <- file.path("..", "figures")

# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(scales)
library(readxl)
library(dplyr)

# Load volcano data for ranking genes
volcano_data <- read.csv(file.path(OUTPUT_PATH, "Cell_State_exp_volcano_LogFC.csv"))
gene_data <- volcano_data[order(volcano_data$LogFC, decreasing = TRUE), ]

# Load metaprogram gene sets
Mprogram <- read_xlsx(file.path(INPUT_PATH, "Meta_Programs.xlsx"))
colnames(Mprogram) <- gsub("^MP\\d+\\s+", "", colnames(Mprogram))

# Convert MProgram into a gene list format for GSEA
metaprograms <- lapply(Mprogram, function(column) {
  na.omit(as.character(column))
})

# Convert list to TERM2GENE format
TERM2GENE <- stack(metaprograms)
names(TERM2GENE) <- c("gene", "term")
TERM2GENE <- TERM2GENE[, c("term", "gene")]

# Prepare ranked gene list for GSEA
gene_list <- gene_data$LogFC
names(gene_list) <- gene_data$Gene

# Add small noise to break ties
gene_list <- gene_list + rnorm(length(gene_list), mean = 0, sd = 1e-5)
# Sort genes in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!is.na(gene_list)]

# Run GSEA
gsea_result <- GSEA(
  geneList = gene_list,
  TERM2GENE = TERM2GENE,
  pvalueCutoff = 0.001,
  minGSSize = 15,
  maxGSSize = 100,
  eps = 1e-10
)

# Filter top 10 enriched pathways
explore_res <- gsea_result@result
top10_gsea <- explore_res %>% arrange(p.adjust) %>% head(10)

# Generate ridge plot for top pathways
pdf(file=file.path(FIGURE_PATH, "Top10_GSEA_Ridgeplot.pdf"), width = 6, height = 5)
ridgeplot(gsea_result) + labs(x = "Enrichment Distribution") +
  scale_fill_gradientn(colors = c("#124985", "#717788", "#a5abbd")) +
  theme(legend.position = "right")
dev.off()

# Generate enrichment plots for top sets
pdf(file=file.path(FIGURE_PATH, "Astro_Enrich.pdf"), width = 6, height = 5)
gseaplot2(gsea_result, geneSetID = 1, title = gsea_result$Description[1])
dev.off()

pdf(file=file.path(FIGURE_PATH, "EMT_Enrich.pdf"), width = 6, height = 5)
gseaplot2(gsea_result, geneSetID = 4, title = gsea_result$Description[4])
dev.off()

