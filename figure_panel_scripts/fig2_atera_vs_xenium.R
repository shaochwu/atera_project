#!/usr/bin/env Rscript
# Fig 2 -- Atera vs Xenium (5K / 280) platform comparison (panels a-f, h-j).
library(Seurat)
library("glmGamPoi")
plan("multisession", workers = 10)
library(ggplot2)

options(future.globals.maxSize = Inf)
options("future.globals.maxSize" = Inf)
Sys.setenv("R_FUTURE_GLOBALS_MAXSIZE" = "Inf")

library(future)
library(magrittr)
library(dplyr)
plan("sequential")
library(tidyr)
library(patchwork)
library(scales)
library(ggtext)
getOption("future.globals.maxSize")



library(arrow)

atera <- readRDS("/path/to/data/atera_Cell_type_and_state.rds")
head(atera@meta.data)

atera$cell_types <- recode(atera$cell_types,
                           "DC" = "Myeloid",
                           "DCIS" = "Cancer Cells",
                           "IBC" = "Cancer Cells",
                           "Enothelial" = "Endothelial",
                           "Pericyte" = "Pericytes")

xenium_5k <- readRDS("/path/to/data/Xenium5K_Cell_type_and_state.rds")
head(xenium_5k@meta.data)
xenium_5k$cell_types <- recode(xenium_5k$cell_types,
                               "DC" = "Myeloid",
                               "DCIS" = "Cancer Cells",
                               "IBC" = "Cancer Cells",
                               "Enothelial" = "Endothelial",
                               "Pericyte" = "Pericytes")


xenium_280 <- readRDS("/path/to/data/xenium_280_s1_bottom.rds")
xenium_280$cell_types <- recode(xenium_280$cell_types,
                                "DC" = "Myeloid",
                                "DCIS" = "Cancer Cells",
                                "IBC" = "Cancer Cells",
                                "Enothelial" = "Endothelial",
                                "Pericyte" = "Pericytes")

xenium_5k <- subset(xenium_5k, subset = cell_types != "Low quality")
atera     <- subset(atera,     subset = cell_types != "Low quality")



######################## fig 2a and b ####################################################################

# No of transcript count per cell 
count_density_df <- bind_rows(
  data.frame(log2count = xenium_280$log2count, cell_types = xenium_280$cell_types, dataset = "Xenium 280"),
  data.frame(log2count = xenium_5k$log2count,  cell_types = xenium_5k$cell_types,  dataset = "Xenium 5k"),
  data.frame(log2count = atera$log2counts,     cell_types = atera$cell_types,      dataset = "Atera")
)

# 2. Shared y-axis scale
count_x_min <- floor(min(count_density_df$log2count, na.rm = TRUE))
count_x_max <- ceiling(max(count_density_df$log2count, na.rm = TRUE))

# 3. Set order: Xenium 280 -> Xenium 5k -> Atera BC
count_density_df$dataset <- factor(count_density_df$dataset,
                                   levels = c("Xenium 280", "Xenium 5k", "Atera"))


count_density_plot <- ggplot(count_density_df, aes(x = log2count, color = dataset, fill = dataset)) +
  geom_density(alpha = 0.25, linewidth = 1.4) +
  scale_fill_manual(values = c("Xenium 280" = "#C85A54", "Xenium 5k" = "#E4A84A", "Atera" = "#2A7F62")) +
  scale_color_manual(values = c("Xenium 280" = "#C85A54", "Xenium 5k" = "#E4A84A", "Atera" = "#2A7F62")) +
  labs(x = "log2(nCount + 1)", y = "Density") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face = "bold", size = 18),
        legend.position = "top",
        axis.line = element_line(linewidth = 1),
        axis.text.x = element_text(face = "bold", size = 18),
        axis.text.y = element_text(face = "bold", size = 18),
        legend.text = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 18)
  ) + labs(x = "log2(nCount+1)")

count_density_plot



##### No. of genes detected per cell
gene_density_df <- bind_rows(
  data.frame(n_genes = xenium_280$nFeature_Xenium, cell_types = xenium_280$cell_types, dataset = "Xenium 280"),
  data.frame(n_genes = xenium_5k$nFeature_Xenium,  cell_types = xenium_5k$cell_types,  dataset = "Xenium 5k"),
  data.frame(n_genes = atera$nFeature_RNA,         cell_types = atera$cell_types,      dataset = "Atera")
)

gene_density_df <- gene_density_df %>% rename(Platform = "dataset")
head(gene_density_df)
gene_y_min <- floor(min(gene_density_df$n_genes, na.rm = TRUE))
gene_y_max <- ceiling(max(gene_density_df$n_genes, na.rm = TRUE))
gene_density_df <- gene_density_df %>%
  mutate(Platform = factor(Platform,
                           levels = c("Xenium 280", "Xenium 5k", "Atera")),
         cell_types = factor(cell_types, levels = rev(c("Basal", "Normal Luminal", "Cancer Cells",
                                                        "Myeloid", "Tcells", "PlasmaCells", "Bcells", "FBs", "Pericytes", "Endothelial"
         ))))

gene_density_plot <- ggplot(gene_density_df, aes(x = log2(n_genes), color = Platform, fill = Platform)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  # facet_wrap(~cell_types, scales = "fixed",  ncol = 5) +
  scale_fill_manual(values = c("Xenium 280" = "#C85A54", "Xenium 5k" = "#E4A84A", "Atera" = "#2A7F62")) +
  scale_color_manual(values = c("Xenium 280" = "#C85A54", "Xenium 5k" = "#E4A84A", "Atera" = "#2A7F62")) +
  labs(x = "log2(nCount + 1)", y = "Density") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face = "bold", size = 18),
        legend.position = "top",
        axis.line = element_line(linewidth = 1),
        axis.text.x = element_text(face = "bold", size = 18),
        axis.text.y = element_text(face = "bold", size = 18),
        legend.text = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 18)
  ) + labs(x = "log2(Number of detected genes per cell)")

gene_density_plot


######################## fig 2C ####################################################################

# Function to build transcript count by gene

make_transcript_bin_df <- function(obj, assay_name, platform_name) {
  
  counts <- GetAssayData(obj, assay = assay_name, layer = "counts")
  max_tx <- max(counts)
  
  cell_summary <- apply(counts, 2, function(x) {
    sapply(1:max_tx, function(i) sum(x == i))
  })
  rownames(cell_summary) <- 1:max_tx
  
  df <- as.data.frame(t(cell_summary))
  df$cell <- rownames(df)
  
  df_long <- df %>%
    pivot_longer(cols      = -cell,
                 names_to  = "transcript_count",
                 values_to = "n_genes") %>%
    mutate(
      transcript_count = as.numeric(transcript_count),
      platform         = platform_name,
      count_bin        = cut(
        transcript_count,
        breaks = c(0, 1, 2, 5, 10, 20, 30, 40, 50, Inf),
        labels = c("1", "2", "3–5", "6–10", "11–20", "21–30",
                   "31–40", "41–50",  ">50"),
        right  = TRUE
      )
    )
  
  return(df_long)
}


# distribution for each platform

transcript_bin_atera <- make_transcript_bin_df(atera,
                                               assay_name    = "RNA",
                                               platform_name = "Atera")

transcript_bin_xenium_5k <- make_transcript_bin_df(xenium_5k,
                                                   assay_name    = "Xenium",
                                                   platform_name = "Xenium 5K")

transcript_bin_xenium_280 <- make_transcript_bin_df(xenium_280,
                                                    assay_name    = "RNA",
                                                    platform_name = "Xenium 280")



# Combine

transcript_bin_df <- bind_rows(transcript_bin_atera, transcript_bin_xenium_5k, transcript_bin_xenium_280) %>%
  mutate(platform = factor(platform,
                           levels = c("Xenium 5K", "Xenium 280", "Atera")))

transcript_bin_df <- bind_rows(transcript_bin_atera, transcript_bin_xenium_5k, transcript_bin_xenium_280) %>%
  mutate(platform = factor(platform,
                           levels = c("Xenium 5K", "Xenium 280", "Atera")))

# checking 
transcript_bin_df %>%
  group_by(platform, count_bin) %>%
  summarise(n = n(), .groups = "drop")


transcript_bin_plot_df <- transcript_bin_df %>%
  mutate(
    log2_genes = log2(n_genes),
    count_bin = factor(count_bin,
                       levels = c("1", "2", "3–5", "6–10", "11–20", "21–30",
                                  "31–40", "41–50",  ">50"))
  )

transcript_bin_plot_df <- transcript_bin_plot_df %>%
  filter(n_genes > 0)

transcript_bin_plot_df <- transcript_bin_plot_df %>% mutate(
  dataset = factor(platform, levels = c("Xenium 280", "Xenium 5k", "Atera")))


transcript_bin_sina_plot <- ggplot(transcript_bin_plot_df, aes(x = count_bin, y = log2_genes, color = platform)) +
  geom_sina(aes(group = interaction(count_bin, platform)),
            alpha = 0.05, size = 0.04, scale = "width", maxwidth = 0.35,
            position = position_dodge(width = 0.6), # Crucial for grouping
            seed = 42) +
# 2. Boxplot
  geom_boxplot(aes(group = interaction(count_bin, platform)),
               width = 0.2, position = position_dodge(width = 0.6), outlier.shape = NA, alpha = 0.5,
               color = "black") +
  scale_color_manual(values = c("Xenium 280" = "#C85A54", "Xenium 5K" = "#E4A84A", "Atera" = "#2A7F62")) +
  theme_classic() + theme_classic() + theme(
    strip.text      = element_text(face = "bold", size = 10),
    axis.title      = element_text(face = "bold", size = 10),
    legend.position = "top",
    axis.text.x     = element_text(face = "bold", size = 10),
    axis.text.y     = element_text(face = "bold", size = 10),
    legend.text     = element_text(face = "bold", size = 10),
    legend.title    = element_text(face = "bold", size = 10),
    axis.line       = element_line(linewidth = 1)
  ) +
  labs(
    x = "Transcript count (binned)",
    y = "log2(Number of genes per cell)"
  )

transcript_bin_sina_plot

######################## Fig 2d: Spatial image plot ####################################################################

xenium_5k_coords <- GetTissueCoordinates(xenium_5k)
head(xenium_5k_coords)
colnames(xenium_5k_coords) <- c("x_centroid", "y_centroid", "cell")
xenium_5k <- AddMetaData(
  xenium_5k,
  metadata = xenium_5k_coords[, c("x_centroid", "y_centroid")])
head(xenium_5k@meta.data[, c("x_centroid", "y_centroid")])

xenium_5k$log2_transcript_counts <- log2(xenium_5k$nCount_Xenium + 1)
xenium_280$log2_transcript_counts <- log2(xenium_280$transcript_counts + 1)
atera$log2_transcript_counts <- log2(atera$transcript_counts + 1)

spatial_all_vals <- c(
  xenium_5k$log2_transcript_counts,
  xenium_280$log2_transcript_counts,
  atera$log2_transcript_counts
)
spatial_lims <- quantile(spatial_all_vals, c(0.05, 0.95), na.rm = TRUE)
spatial_mid_val <- mean(spatial_lims) # midpoint for the diverging scale
spatial_bwr_scale <- scale_color_gradient2(
  low      = "blue",
  mid      = "white",
  high     = "red",
  midpoint = spatial_mid_val,
  limits   = spatial_lims,
  oob      = scales::squish,
  name     = "log2 transcript\ncounts"
)

plot_spatial <- function(obj, title) {
  ggplot(
    obj@meta.data,
    aes(x = x_centroid, y = y_centroid, color = log2_transcript_counts)
  ) +
    geom_point(size = 0.15) +
    bwr_scale +
    coord_fixed() +
    scale_y_reverse() +
    theme_classic() +
    theme(
      panel.background  = element_rect(fill = "white"),
      plot.background   = element_rect(fill = "white"),
      axis.line         = element_line(color = "black"), 
      axis.text = element_text(color = "black", face = "bold"), # visible
      axis.ticks        = element_line(color = "black"),   # visible
      axis.title        = element_text(color = "black", face = "bold"),   # visible
      plot.title        = element_text(color = "black", hjust = 0.5, face = "bold", size = 12),
      legend.background = element_rect(fill = "white"),
      legend.text       = element_text(color = "black"),
      legend.title      = element_text(color = "black"),
    ) +
    labs(title = title, x = "X centroid", y = "Y centroid")
}

p1 <- ggplot(
  xenium_5k@meta.data,
  aes(x = y_centroid, y = x_centroid, color = log2_transcript_counts)
) +
  geom_point(size = 0.15) +
  bwr_scale +
  coord_fixed() +
  scale_y_reverse() +
  theme_classic() +theme( plot.title= element_text(color = "black", hjust = 0.5, face = "bold", size = 12),
                          axis.title        = element_text(color = "black", face = "bold"), 
                          axis.text = element_text(color = "black", face = "bold"))+
  labs(title = "Xenium 5K", x = "X centroid", y = "Y centroid")

p2 <- plot_spatial(xenium_280,        "Xenium 280")
p3 <- plot_spatial(atera,          "Atera")

# Equal widths for all three panels
spatial_log2_count_fig <-   p2 +p1+ p3 +
  plot_layout(guides = "collect", widths = c(1, 1, 1)) &
  theme(legend.position = "right")


spatial_log2_count_fig


######################## Fig 1e and f ####################################################################


celltype_order <- c("Basal", "Normal Luminal", "Cancer Cells", "Myeloid", "Tcells", "PlasmaCells", "FBs", "Bcells", "Pericytes", "Endothelial")

celltype_cols <- c("Basal" = "#F08A80", "Normal Luminal" = "#4DB3E6", "Cancer Cells" = "#B7B51A", "Myeloid" = "#41C0C5", "Tcells" = "#E56AB3", "PlasmaCells" = "#D86AE3", "Bcells" = "#D89B1D",
                   "FBs" = "#43C08B", "Pericytes" = "#8F88F0", "Endothelial" = "#57C53B")

celltype_labels <- paste0("<span style='color:", celltype_cols, "'>", names(celltype_cols), "</span>")
names(celltype_labels) <- names(celltype_cols)


######### AVG transcript per cell types and per cells
transcript_ridge_df <- bind_rows(
  data.frame(transcript_count = xenium_280$nCount_Xenium, cell_types = xenium_280$cell_types, dataset = "Xenium 280"),
  data.frame(transcript_count = xenium_5k$nCount_Xenium,  cell_types = xenium_5k$cell_types,  dataset = "Xenium 5k"),
  data.frame(transcript_count = atera$nCount_RNA,         cell_types = atera$cell_types,      dataset = "Atera")
)
transcript_ridge_df <- transcript_ridge_df %>% filter(!is.na(cell_types), cell_types != "Low quality")
transcript_ridge_df <- transcript_ridge_df %>% mutate(dataset = factor(dataset, levels = c("Xenium 280", "Xenium 5k", "Atera")))
transcript_ridge_df <- transcript_ridge_df %>%
  mutate(dataset = factor(dataset,
                          levels = c("Xenium 280", "Xenium 5k", "Atera")),
         cell_types = factor(cell_types, levels = rev(c("Basal", "Normal Luminal", "Cancer Cells",
                                                        "Myeloid", "Tcells", "PlasmaCells", "Bcells", "FBs", "Pericytes", "Endothelial"
         ))))

transcript_ridge_plot <- ggplot(transcript_ridge_df,
                                aes(x = log2(transcript_count), y = cell_types, fill = dataset)) + geom_density_ridges(alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("Xenium 280" = "#C85A54", "Xenium 5k" = "#E4A84A", "Atera" = "#2A7F62")) +
  scale_y_discrete(labels = celltype_labels) +
  theme_classic(base_size = 12) +
  labs(x = "log2(Trancript count)", y = "Cell Types",
       fill = "dataset") +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = ggtext::element_markdown(face = "bold", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    legend.position = "top",
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12)
  )

transcript_ridge_plot


##. Gene count by cell types

gene_ridge_df <- gene_density_df %>% filter(!is.na(cell_types)) %>%
  mutate(Platform = factor(
    Platform, levels = c("Xenium 280", "Xenium 5k", "Atera")),
    cell_types = factor(cell_types,
                        levels = rev(celltype_order)))

gene_ridge_plot <- ggplot(gene_ridge_df,
                          aes(x = log2(n_genes), y = cell_types, fill = Platform)) + geom_density_ridges(alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("Xenium 280" = "#C85A54", "Xenium 5k" = "#E4A84A", "Atera" = "#2A7F62")) +
  scale_y_discrete(labels = celltype_labels) +
  theme_classic(base_size = 12) +
  labs(x = "log2(Number of detected genes per cell)", y = "Cell Types",
       fill = "Platform") +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = ggtext::element_markdown(face = "bold", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    legend.position = "top",
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12)
  )

gene_ridge_plot



######################## Fig 2e ####################################################################
ct_col <- "cell_types"
# Compute cell-type proportions in each object
get_proportions <- function(obj, col) {
  tab <- table(obj@meta.data[[col]])
  prop <- as.data.frame(tab / sum(tab))
  colnames(prop) <- c("cell_type", "proportion")
  prop
}

prop_atera      <- get_proportions(atera, ct_col)
prop_xenium_280 <- get_proportions(xenium_280, ct_col)

# Keep only shared cell types
prop_corr_df <- inner_join(prop_atera, prop_xenium_280, by = "cell_type",
                           suffix = c("_Atera", "_Xenium280"))

# Spearman correlation
prop_corr_r <- cor(prop_corr_df$proportion_Atera, prop_corr_df$proportion_Xenium280, method = "spearman")
prop_corr_p <- cor.test(prop_corr_df$proportion_Atera, prop_corr_df$proportion_Xenium280,
                        method = "spearman")$p.value

celltype_cols <- c("Basal" = "#F08A80", "Normal Luminal" = "#4DB3E6", "Cancer Cells" = "#B7B51A",
                   "Myeloid" = "#41C0C5", "Tcells" = "#E56AB3", "PlasmaCells" = "#D86AE3", "Bcells" = "#D89B1D",
                   "FBs" = "#43C08B", "Pericytes" = "#8F88F0", "Endothelial" = "#57C53B")



celltype_labels <- paste0(
  "<span style='color:", celltype_cols, "'>",
  names(celltype_cols),
  "</span>"
)
names(celltype_labels) <- names(celltype_cols)

prop_corr_plot <- ggplot(prop_corr_df, aes(x = proportion_Atera, y = proportion_Xenium280)) +
  geom_point(aes(colour = cell_type), size = 8, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_smooth(method = "glm", se = FALSE, colour = "black", linewidth = 0.8) +
  annotate("text",
           x = min(prop_corr_df$proportion_Atera) * 1.5,
           y = max(prop_corr_df$proportion_Xenium280) * 0.9,
           label = sprintf("R = %.2f, P = %.3f", prop_corr_r, prop_corr_p), hjust = 0, size = 6,
           fontface = "italic") +
  scale_colour_manual(values = celltype_cols, labels = celltype_labels) +
  guides(colour = guide_legend(override.aes = list(size = 0, alpha = 0))) +
  labs(
    x = "Cell-type proportion (Atera)",
    y = "Cell-type proportion (Xenium 280)",
    colour = NULL) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", face = "bold", size = 22),
    # axis.ticks = element_line(color = "black", size = 14),
    axis.line = element_line(color = "black", linewidth = 1.5),
    axis.title = element_text(color = "black", face = "bold", size = 16),
    legend.text = ggtext::element_markdown(face = "bold", size = 24),
    legend.title = element_blank(),
    plot.title = element_text(color = "black", hjust = 0.5, face = "bold", size = 12),
    legend.box.margin = margin(0, 0, 0, -90)
  )

prop_corr_plot


######################## Fig 2h-j ####################################################################

#. Marker gene transcript-count distributions across platforms

make_gene_bin_celltype_df <- function(obj,
                                      assay_name,
                                      platform_name,
                                      celltype_col = "cell_types") {
  
  counts <- GetAssayData(obj, assay = assay_name, layer = "counts")
  
  celltypes <- obj@meta.data[[celltype_col]]
  names(celltypes) <- rownames(obj@meta.data)
  
  x <- summary(counts)
  
  df <- tibble(
    gene = rownames(counts)[x$i],
    cell = colnames(counts)[x$j],
    transcript_count = x$x,
    celltype = celltypes[colnames(counts)[x$j]],
    platform = platform_name
  ) %>%
    mutate(
      count_bin = cut(
        transcript_count,
        breaks = c(0, 1, 2, 5, 10, 20, 30, 40, 50, Inf),
        labels = c("1", "2", "3–5", "6–10", "11–20", "21–30",
                   "31–40", "41–50",  ">50"),
        right = TRUE
      )
    )
  
  df
}

gene_bin_atera <- make_gene_bin_celltype_df(
  atera,
  assay_name = "RNA",
  platform_name = "Atera",
  celltype_col = "cell_types"
)

gene_bin_xenium_5k <- make_gene_bin_celltype_df(
  xenium_5k,
  assay_name = "Xenium",
  platform_name = "Xenium 5K",
  celltype_col = "cell_types"
)

gene_bin_xenium_280 <- make_gene_bin_celltype_df(
  xenium_280,
  assay_name = "RNA",
  platform_name = "Xenium 280",
  celltype_col = "cell_types"
)

gene_bin_df <- bind_rows(gene_bin_atera, gene_bin_xenium_5k, gene_bin_xenium_280)


# Common gene markers across all the three platforms
gene_order <- c(
  'LAMA3', 'LAMB3', 'LAMC2',   ## Basal
  "GATA3", "FOXA1", "ANKRD30A", "KIT", ## Normal Luminal
  "EPCAM", "VEGFA", "VTCN1", "MDM2", ## Cancer Cells
  "CD68", "CD163", "ITGAX", "CHIT1", ## Myeloid
  "TRAC", "CD3E", "IL7R", "CXCR4", "TIGIT", ## T cells
  "TENT5C", "TNFRSF13C", ## Plasma cells
  "MS4A1", "BANK1", "CD27", # Bcels
  "COL4A1", "SFRP1", "SFRP4",  ## Fibroblasts
  "RGS5", "PDGFRB", "MYH11", "MYLK", ## Pericytes
  "PECAM1", "FLT1", "PLVAP", "CLEC14A", "AQP1" ## Endothelial
)

celltype_order <- c(
  "Basal",
  "Normal Luminal",
  "Cancer Cells",
  "Myeloid",
  "Tcells",
  "PlasmaCells",
  "Bcells",
  "FBs",
  "Pericytes",
  "Endothelial"
)

celltype_cols <- c("Basal" = "#F08A80", "Normal Luminal" = "#4DB3E6", "Cancer Cells" = "#B7B51A",
                   "Myeloid" = "#41C0C5", "Tcells" = "#E56AB3", "PlasmaCells" = "#D86AE3", "Bcells" = "#D89B1D",
                   "FBs" = "#43C08B", "Pericytes" = "#8F88F0", "Endothelial" = "#57C53B")


marker_all_df <- gene_bin_df %>% filter(platform %in% c("Atera", "Xenium 280", "Xenium 5K"),
                                        gene %in% gene_order) %>%
  mutate(platform = factor(platform, levels = c("Xenium 280", "Xenium 5K", "Atera")),
         celltype = factor(celltype, levels = celltype_order),
         gene = factor(gene, levels = gene_order))



marker_all_plot <- ggplot(marker_all_df,
                          aes(x = platform, y = log2(transcript_count + 1), fill = platform)) +
  geom_violin(scale = "width",
              color = "black", linewidth = 0.2) +
  facet_grid(gene ~ celltype, scales = "free_y", switch = "y") +
  scale_y_continuous(position = "right") +
  theme_classic(base_size = 13) +
  theme(
    strip.placement = "outside",
    strip.background.y = element_blank(),
    strip.background.x = element_rect(
      fill = "grey90",
      colour = "black",
      linewidth = 0.3
    ),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    strip.text.x = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y.right = element_text(size = 12, face = "bold"),
    legend.position = "none",
    
    # space between cell-type panels
    panel.spacing.x = unit(1, "cm")
  ) + scale_fill_manual(
    values = c(
      "Xenium 280" = "#C85A54",
      "Xenium 5K"  = "#E4A84A",
      "Atera"      = "#2A7F62"
    )
  ) +
  labs(
    x = "Platform",
    y = "log2(transcript count per cell)"
  )

marker_all_plot



# Xenium 280
gene_order_xenium_280 <- c(
  "KRT14", "KRT5", "MYLK", "COL17A1", "ACTA2",       ## Basal
  "KIT", "VTCN1", "KRT23", "KRT6B", "PTN", "AGR3",
  "FOXA1", "MLPH", "SCUBE2", "ESR1",                 ## Normal Luminal
  "CD68", "CD74", "CD163", "LYZ", "ITGAX",           ## Myeloid
  "TRAC", "CD3E", "CXCR4", "IL2RG", "IL7R", "CD27",  ## T cells
  "TENT5C", "CYTIP",                                 ## Plasma cells
  "PECAM1", "BIRC3", "MS4A1", "BANK1", "TNFRSF13C",
  "ARHGDIB", "GPR183",                               ## B cells
  "MMP2", "LAMB1", "LAMA2", "PDGFRB", "C1R",         ## Fibroblasts
  "RGS5", "NDUFA4L2", "COL4A1", "IGFBP5", "HEYL",    ## Pericytes
  "VWF", "FLT1", "PLVAP", "AQP1", "CLEC14A"          ## Endothelial
)

#  Atera
gene_order_atera <- c(
  "TNS4", "MYLK", "COL17A1", "TRIM29", "TP63",       ## Basal
  "GABRP", "CCL28", "PROM1", "KIT", "EHF",
  "ELAPOR1", "AGR3", "CA12", "ADCY1", "SCUBE2",      ## Normal Luminal
  "CSF1R", "CYBB", "C1QA", "C1QC", "CD14",           ## Myeloid
  "TRAC", "CD3E", "TMC8", "TRBC2", "ACAP1",          ## T cells
  "MZB1", "SLAMF7", "DERL3", "POU2AF1", "JCHAIN",    ## Plasma cells
  "MS4A1", "BANK1", "CD37", "CD22",                  ## B cells
  "LTB", "LUM", "GEM", "ABCA6", "THBS2", "COL6A3",   ## Fibroblasts
  "EPAS1", "CAV1", "NOTCH3", "NEURL1B", "MYO1B",     ## Pericytes
  "VWF", "SHANK3", "DDHD", "AQP1", "PLVAP"           ## Endothelial
)

# Xenium 5K
gene_order_xenium_5k <- c(
  "DST", "THBS1", "MYL9", "MYH11", "MYLK",           ## Basal
  "CA12", "TSPAN13", "GATA3", "HSPB8", "DDR1",       ## Cancer Cells / Luminal
  "NAMPT", "MS4A6A", "C1QC", "CTSL",                 ## Myeloid
  "IL7R", "TRAC", "TRBC1", "HNRNPH1", "EEF1G",       ## T cells
  "PIM2", "TENT5C", "XBP1",                          ## Plasma cells
  "MS4A1", "CD52",                                   ## B cells
  "POSTN", "DCN", "LUM", "CCL2", "BGN", "C11orf96",  ## Fibroblasts
  "COL4A1", "EPAS1", "COL4A2",                       ## Pericytes
  "AQP1", "ADAMTS1", "FLT1", "PLVAP"                 ## Endothelial
)


marker_atera_df <- gene_bin_df %>% filter(platform == "Atera",
                                          gene %in% gene_order_atera) %>%
  mutate(celltype = factor(celltype, levels = celltype_order),
         gene = factor(gene, levels = (gene_order_atera)))

marker_atera_plot <- ggplot(marker_atera_df, aes(x = "", y = log2(transcript_count + 1), fill = celltype
)) +
  geom_violin(scale = "width", color = "black", linewidth = 0.2) +
  facet_grid(gene ~ celltype, scales = "free_y", switch = "y", labeller = labeller(celltype = celltype_labels)
  ) +
  scale_y_continuous(position = "right") + scale_fill_manual(values = celltype_cols) +
  theme_classic(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.background.y = element_blank(),
    strip.background.x = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    strip.text.x = ggtext::element_markdown(
      hjust = 0.5, face = "bold", size = 12, angle = 90),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y.right = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.spacing.x = unit(1, "cm")) +
  labs(x = NULL,
       y = "log2(transcript count per cell + 1)")


marker_atera_plot

marker_xenium_280_df <- gene_bin_df %>% filter(platform == "Xenium 280", gene %in% gene_order_xenium_280) %>%
  mutate(celltype = factor(celltype, levels = celltype_order),
         gene = factor(gene, levels = (gene_order_xenium_280)))


marker_xenium_280_plot <- ggplot(marker_xenium_280_df, aes(x = "", y = log2(transcript_count + 1), fill = celltype
)) +
  geom_violin(scale = "width", color = "black", linewidth = 0.2) +
  facet_grid(gene ~ celltype, scales = "free_y", switch = "y", labeller = labeller(celltype = celltype_labels)
  ) +
  scale_y_continuous(position = "right") + scale_fill_manual(values = celltype_cols) +
  theme_classic(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.background.y = element_blank(),
    strip.background.x = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    strip.text.x = ggtext::element_markdown(
      angle = 90, hjust = 0, vjust = 0.5, face = "bold", size = 12),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y.right = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.spacing.x = unit(1, "cm")) +
  labs(x = NULL,
       y = "log2(transcript count per cell + 1)")


marker_xenium_280_plot


marker_xenium_5k_df <- gene_bin_df %>% filter(platform == "Xenium 5K", gene %in% gene_order_xenium_5k) %>%
  mutate(celltype = factor(celltype, levels = celltype_order),
         gene = factor(gene, levels = (gene_order_xenium_5k)))

celltype_cols <- c("Basal" = "#F08A80", "Normal Luminal" = "#4DB3E6", "Cancer Cells" = "#B7B51A",
                   "Myeloid" = "#41C0C5", "Tcells" = "#E56AB3", "PlasmaCells" = "#D86AE3", "Bcells" = "#D89B1D",
                   "FBs" = "#43C08B", "Pericytes" = "#8F88F0", "Endothelial" = "#57C53B")

celltype_labels <- paste0("<span style='color:", celltype_cols, "'>",
                          names(celltype_cols), "</span>")

names(celltype_labels) <- names(celltype_cols)

marker_xenium_5k_df <- gene_bin_df %>% filter(
  platform == "Xenium 5K", gene %in% gene_order_xenium_5k) %>%
  mutate(celltype = factor(celltype, levels = celltype_order),
         gene = factor(gene, levels = gene_order_xenium_5k))

marker_xenium_5k_plot <- ggplot(marker_xenium_5k_df, aes(x = "", y = log2(transcript_count + 1), fill = celltype
)) +
  geom_violin(scale = "width", color = "black", linewidth = 0.2) +
  facet_grid(gene ~ celltype, scales = "free_y", switch = "y", labeller = labeller(celltype = celltype_labels)
  ) +
  scale_y_continuous(position = "right") + scale_fill_manual(values = celltype_cols) +
  theme_classic(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.background.y = element_blank(),
    strip.background.x = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    strip.text.x = ggtext::element_markdown(
      angle = 90, hjust = 0, vjust = 0.5, face = "bold", size = 12),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y.right = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.spacing.x = unit(1, "cm")) +
  labs(x = NULL,
       y = "log2(transcript count per cell + 1)")

marker_xenium_5k_plot