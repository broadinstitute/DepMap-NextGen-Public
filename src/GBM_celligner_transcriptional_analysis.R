### Celligner-Based GBM Transcriptional Analysis
# Purpose: This script generates multiple figures for GBM transcriptional profiles, including UMAP, stacked bar, SOX2, OLIG1, and volcano plots.

# Define input, output, and figure directories
INPUT_PATH <- file.path("..", "data")
OUTPUT_PATH <- file.path("..", "processed")
FIGURE_PATH <- file.path("..", "figures")

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Load datasets
celligner <- read.csv(
  file.path(INPUT_PATH,
            "celligner_coordinates_w_hcmi.csv"),
  stringsAsFactors = FALSE
)

model_info <- read.csv(
  file.path(INPUT_PATH,
            "model_metadata.csv"),
  stringsAsFactors = FALSE
)

# Map annotations
celligner_gbm <- celligner %>%
  dplyr::filter(
    lineage == "CNS/Brain",
    subtype %in% c("Glioblastoma", "Glioblastoma Multiforme", "Gliosarcoma", "Diffuse Glioma")
  )

# Map NextGen / Traditional using IsNextGen in model metadata
celligner_gbm$modelType <- model_info$IsNextGen[match(celligner_gbm$ModelID, model_info$ModelID)]
celligner_gbm$modelType <- gsub("True",  "NextGen",     celligner_gbm$modelType)
celligner_gbm$modelType <- gsub("False", "Traditional", celligner_gbm$modelType)

# Treat Diffuse Glioma tumors as "Tumor"
celligner_gbm$modelType <- ifelse(
  celligner_gbm$subtype == "Diffuse Glioma",
  "Tumor",
  celligner_gbm$modelType
)

# Define TCGA vs DepMap shapes
celligner_gbm$shape <- ifelse(
  grepl("tcga", celligner_gbm$type, ignore.case = TRUE),
  "triangle",
  "circle"
)

# Colors
celligner_gbm <- celligner_gbm %>%
  mutate(custom_color = dplyr::case_when(
    modelType == "Traditional" ~ "#2b6999",  # GBM Traditional line
    modelType == "NextGen"     ~ "#F58900",  # GBM NextGen line
    modelType == "Tumor"       ~ "#93ae93",  # GBM Tumor
    TRUE                       ~ "gray80"    # Others
  ))

# Generate and save UMAP Plot
pdf(file = file.path(FIGURE_PATH, "Celligner_organoid.pdf"),
    width = 6.1, height = 4.4)

umap_plot <- ggplot(celligner_gbm, aes(x = umap1, y = umap2, color = custom_color, shape = shape)) +
  # Background "other" points
  geom_point(
    data = subset(celligner_gbm, custom_color == "gray80"),
    aes(x = umap1, y = umap2),
    size = 1, alpha = 0.6, color = "gray80"
  ) +
  # Tumors (green triangles)
  geom_jitter(
    data = subset(celligner_gbm, custom_color == "#93ae93"),
    aes(x = umap1, y = umap2, color = custom_color),
    size = 2, alpha = 0.6, shape = 17
  ) +
  # NextGen (orange circles + outline)
  geom_jitter(
    data = subset(celligner_gbm, custom_color == "#F58900"),
    aes(x = umap1, y = umap2, color = custom_color),
    size = 2, alpha = 0.6, shape = 16
  ) +
  geom_jitter(
    data = subset(celligner_gbm, custom_color == "#F58900"),
    aes(x = umap1, y = umap2),
    size = 2, alpha = 0.6, shape = 1, stroke = 0.3, color = "#F58900"
  ) +
  # Traditional (blue circles + outline)
  geom_jitter(
    data = subset(celligner_gbm, custom_color == "#2b6999"),
    aes(x = umap1, y = umap2, color = custom_color),
    size = 2, alpha = 0.6, shape = 16
  ) +
  geom_jitter(
    data = subset(celligner_gbm, custom_color == "#2b6999"),
    aes(x = umap1, y = umap2),
    size = 2, alpha = 0.6, shape = 1, stroke = 0.3, color = "#2b6999"
  ) +
  scale_shape_manual(
    values = c("circle" = 16, "triangle" = 17),
    labels = c("circle" = "DepMap", "triangle" = "TCGA")
  ) +
  scale_color_manual(
    values = c("#2b6999", "#93ae93", "#F58900", "gray80"),
    labels = c("GBM Traditional Line", "GBM Tumor", "GBM NextGen Line", "Others")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16)),
    shape = guide_legend(override.aes = list(color = "black"))
  ) +
  labs(x = "UMAP 1", y = "UMAP 2", color = NULL, shape = NULL) +
  ggpubr::theme_pubr() +
  theme(
    text          = element_text(size = 12, color = "black"),
    legend.text   = element_text(size = 12, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    legend.position = "right"
  )

print(umap_plot)
dev.off()

# DepMap GBM models only (for downstream plots)
filter_GBM <- celligner_gbm %>%
  dplyr::filter(
    type == "DepMap Model",
    subtype %in% c("Glioblastoma", "Glioblastoma Multiforme")
  )

# Define coarse cluster based on UMAP2
filter_GBM$cluster <- ifelse(filter_GBM$umap2 > 9.5, "Glial", "Mesenchymal")

# Stacked bar frequencies of model growth patterns/clusters
p_growth_labels <- ggplot(filter_GBM, aes(x = modelType, fill = GrowthPattern)) +
  geom_bar(position = "fill", stat = "count") +
  facet_wrap(~ cluster) +
  labs(
    x   = "Model Type",
    y   = "Proportion",
    fill = "Growth Pattern"
  ) +
  ggpubr::theme_pubr() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5),
    legend.position = "right"
  )+
  geom_text(
    stat = "count",
    aes(label = ..count..),
    position = position_fill(vjust = 0.5),
    size = 4, color = "white"
  ) +
  scale_fill_manual(
    values = c(
      "Spheroid" = "#FF8C00",
      "Adherent" = "#36648B",
      "Mixed"    = "gray70"
    )
  )

ggsave(
  file.path(FIGURE_PATH, "freq_model_cluster.pdf"),
  p_growth_labels,
  width = 4, height = 3.9, units = "in", dpi = 300
)

# SOX2 / OLIG1 rank plots
exp_nextgen_path <- file.path(INPUT_PATH, "next_gen_expression.csv")
exp_trad_path <- file.path(INPUT_PATH, "traditional_expression.csv") 

exp_nextgen <- read.csv(
  exp_nextgen_path,
  stringsAsFactors = FALSE
)

exp_trad <- read.csv(
  exp_trad_path,
  stringsAsFactors = FALSE
)

exp_all <- rbind(exp_nextgen, exp_trad)

colnames(exp_all) <- gsub("\\..*$", "", colnames(exp_all))

# SOX2 and OLIG1
exp_depmap <- exp_all[, c("ModelID", "SOX2", "OLIG1")]
rownames(exp_depmap) <- exp_depmap$ModelID
exp_depmap$ModelID <- NULL
exp_depmap$type <- filter_GBM$modelType[match(rownames(exp_depmap), filter_GBM$ModelID)]
exp_depmap$type <- ifelse(is.na(exp_depmap$type), "Other", exp_depmap$type)

# Rank by SOX2 and OLIG1
exp_depmap <- exp_depmap[order(exp_depmap$SOX2, decreasing = TRUE), ]
exp_depmap$SOX2_rank <- seq_len(nrow(exp_depmap))
exp_depmap <- exp_depmap[order(exp_depmap$OLIG1, decreasing = TRUE), ]
exp_depmap$OLIG1_rank <- seq_len(nrow(exp_depmap))

# OLIG1 plot
p_olig1 <- ggplot(exp_depmap, aes(x = OLIG1_rank, y = OLIG1)) +
  geom_point(data = subset(exp_depmap, type == "Other"),
             color = "grey", size = 2, alpha = 0.6) +
  geom_point(data = subset(exp_depmap, type == "Traditional"),
             color = "#36648B", size = 2, alpha = 0.8) +
  geom_point(data = subset(exp_depmap, type == "NextGen"),
             color = "#FF8C00", size = 2, alpha = 0.8) +
  labs(title = "OLIG1 Expression", x = "Rank", y = "Expression log(TPM)") +
  ggpubr::theme_pubr()

ggsave(
  file.path(FIGURE_PATH, "OLIG1.pdf"),
  p_olig1,
  width = 4.7, height = 3.9, units = "in", dpi = 300
)

# SOX2 plot
p_sox2 <- ggplot(exp_depmap, aes(x = SOX2_rank, y = SOX2)) +
  geom_point(data = subset(exp_depmap, type == "Other"),
             color = "grey", size = 2, alpha = 0.6) +
  geom_point(data = subset(exp_depmap, type == "Traditional"),
             color = "#36648B", size = 2, alpha = 0.8) +
  geom_point(data = subset(exp_depmap, type == "NextGen"),
             color = "#FF8C00", size = 2, alpha = 0.8) +
  labs(title = "SOX2 Expression", x = "Rank", y = "Expression log(TPM)") +
  ggpubr::theme_pubr()

ggsave(
  file.path(FIGURE_PATH, "SOX2.pdf"),
  p_sox2,
  width = 4.7, height = 3.9, units = "in", dpi = 300
)

# Glial vs Mesenchymal Gene Expression
# Glial / Mesenchymal IDs
Glial_IDs <- filter_GBM %>%
  dplyr::filter(cluster == "Glial", modelType != "Traditional") %>%
  dplyr::pull(ModelID)

Mes_IDs <- filter_GBM %>%
  dplyr::filter(cluster == "Mesenchymal") %>%
  dplyr::pull(ModelID)

exp_volcano <- exp_all  # already cleaned colnames, deduped above

# Subset to glial/mes IDs
exp_filtered <- exp_volcano[exp_volcano$ModelID %in% c(Glial_IDs, Mes_IDs), ]

# Add cluster label
exp_filtered$cluster <- ifelse(exp_filtered$ModelID %in% Glial_IDs, "Glial", "Mesenchymal")

# Identify non-gene meta columns to drop
meta_cols <- intersect(colnames(exp_filtered),
                       c("DepMapID", "ModelID", "lineage", "subtype",
                         "type", "modelType", "GrowthPattern", "cluster"))

expression_data <- exp_filtered[, setdiff(colnames(exp_filtered), meta_cols)]
cluster_labels  <- exp_filtered$cluster

# Initialize vectors
logFC    <- numeric(ncol(expression_data))
p_values <- numeric(ncol(expression_data))

# Compute logFC and p-values
for (i in seq_len(ncol(expression_data))) {
  glial_values <- expression_data[cluster_labels == "Glial", i]
  mes_values   <- expression_data[cluster_labels == "Mesenchymal", i]
  
  if (length(glial_values[!is.na(glial_values)]) > 1 &&
      length(mes_values[!is.na(mes_values)])   > 1) {
    
    mean_glial <- mean(glial_values, na.rm = TRUE)
    mean_mes   <- mean(mes_values,   na.rm = TRUE)
    logFC[i]   <- log2(mean_glial + 1) - log2(mean_mes + 1)
    
    t_test     <- t.test(glial_values, mes_values)
    p_values[i] <- t_test$p.value
  } else {
    logFC[i]    <- NA
    p_values[i] <- NA
  }
}

# Build volcano data frame
valid_idx <- !is.na(logFC) & !is.na(p_values)

volcano_data <- data.frame(
  Gene      = colnames(expression_data)[valid_idx],
  LogFC     = logFC[valid_idx],
  # first pass FDR to get logP, then re-do FDR cleanly below
  LogPValue = -log10(p.adjust(p_values[valid_idx], method = "fdr"))
)

# Back-transform to raw p-values and re-do FDR
volcano_data$raw_p_value <- 10^(-volcano_data$LogPValue)
volcano_data$q_value     <- p.adjust(volcano_data$raw_p_value, method = "fdr")

# Significance categories
volcano_data$Significance <- ifelse(
  volcano_data$q_value < 0.001 & volcano_data$LogFC >  1, "Significant Positive",
  ifelse(volcano_data$q_value < 0.001 & volcano_data$LogFC < -1, "Significant Negative",
         "Not Significant")
)

# Top 10 up / down
top_right_genes <- volcano_data %>%
  dplyr::filter(q_value < 0.001, LogFC > 1) %>%
  dplyr::arrange(q_value) %>%
  head(10)

top_left_genes <- volcano_data %>%
  dplyr::filter(q_value < 0.001, LogFC < -1) %>%
  dplyr::arrange(q_value) %>%
  head(10)

top_genes <- rbind(top_right_genes, top_left_genes)

# Plot
p_volcano <- ggplot(volcano_data, aes(x = LogFC, y = -log10(q_value), label = Gene)) +
  geom_point(aes(color = Significance), alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Significant Positive" = "#727400",
      "Significant Negative" = "#2f4858",
      "Not Significant"      = "grey"
    )
  ) +
  ggrepel::geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 3, max.overlaps = Inf, box.padding = 0.3, point.padding = 0.3
  ) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(Gene %in% c("SOX2", "OLIG1"), as.character(Gene), "")),
    size = 3, fontface = "bold", max.overlaps = Inf,
    box.padding = 0.5, point.padding = 0.5
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-log10(q-value)"
  ) +
  ggpubr::theme_pubr() +
  theme(
    text            = element_text(size = 12),
    legend.position = "none"
  )

ggsave(
  file.path(FIGURE_PATH, "Cell_State_Volcano_LogFC.pdf"),
  p_volcano,
  width = 4.9, height = 4.1, units = "in", dpi = 300
)

# Save data
write.csv(
  volcano_data,
  file.path(OUTPUT_PATH, "Cell_State_exp_volcano_LogFC.csv"),
  row.names = FALSE
)

