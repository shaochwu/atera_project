#!/usr/bin/env Rscript
# Combined 3-platform UMAP figure (stacked, shared cell-type legend).
suppressMessages({library(ggplot2); library(ggrastr); library(patchwork)})

CACHE <- "/path/to/platform_comparison/results/manuscript_cache"
OUT   <- "/path/to/atera_demo/results/infercnv_run/manuscript_panels/figure1_panelF_umaps"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

plat_labels <- c(atera = "Atera", Xenium280 = "Xenium 280", Xenium5K = "Xenium 5K")
plat_cols   <- c(Atera = "#4FA381", "Xenium 280" = "#D95F52", "Xenium 5K" = "#E5A552")
order_plat  <- c("Xenium280","Xenium5K","atera")

harmonize_ct <- function(ct) {
  ct <- as.character(ct)
  map <- c(
    "DCIS"="Cancer Cells","IBC"="Cancer Cells","Cancer Cells"="Cancer Cells","Cancer cells"="Cancer Cells",
    "Basal"="Basal",
    "Normal Luminal"="Normal Luminal","Normal luminal"="Normal Luminal","Luminal"="Normal Luminal",
    "Myeloid"="Myeloid","DC"="Myeloid","Macrophage"="Myeloid","Macrophages"="Myeloid",
    "Tcells"="Tcells","T cells"="Tcells","Tcell"="Tcells",
    "Bcells"="Bcells","B cells"="Bcells","Bcell"="Bcells",
    "PlasmaCells"="PlasmaCells","Plasma Cells"="PlasmaCells","Plasma cells"="PlasmaCells","Plasmablast"="PlasmaCells",
    "FBs"="FBs","Fibroblasts"="FBs","Fibroblast"="FBs","FB"="FBs",
    "Pericyte"="Pericytes","Pericytes"="Pericytes",
    "Endothelial"="Endothelial","Enothelial"="Endothelial","Endothelium"="Endothelial")
  out <- unname(map[ct]); out[is.na(out)] <- "Low quality"; out
}
ct_cols <- c(
  "Basal" = "#F08A80","Normal Luminal" = "#4DB3E6","Cancer Cells" = "#B7B51A",
  "Myeloid" = "#41C0C5","Tcells" = "#E56AB3","PlasmaCells" = "#D86AE3",
  "Bcells" = "#D89B1D","FBs" = "#43C08B","Pericytes" = "#8F88F0","Endothelial" = "#57C53B",
  "Low quality" = "#9E9E9E")
ct_order <- c("Cancer Cells","Basal","Normal Luminal","FBs","Pericytes","Endothelial",
              "Myeloid","Tcells","Bcells","PlasmaCells","Low quality")  # Low quality last

save_pdf <- function(plot, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), plot, width = w, height = h,
         device = "pdf", useDingbats = FALSE, dpi = 300)
  ggsave(file.path(OUT, paste0(name, ".png")), plot, width = w, height = h, dpi = 200)
}

umap_plot <- function(p, show_legend) {
  df <- readRDS(file.path(CACHE, paste0("umap_", p, ".rds")))
  df$cell_type <- harmonize_ct(df$cell_type)
  df <- df[df$cell_type %in% ct_order, ]
  df$cell_type <- factor(df$cell_type, levels = ct_order)   # all types -> shared legend
  df <- df[order(df$cell_type != "Low quality"), ]          # draw Low quality underneath
  lab <- plat_labels[[p]]
  ptsize <- if (nrow(df) > 1e5) 0.2 else 0.4
  ggplot(df, aes(UMAP_1, UMAP_2, color = cell_type)) +
    rasterise(geom_point(size = ptsize, stroke = 0), dpi = 300, dev = "ragg") +
    scale_color_manual(values = ct_cols, drop = FALSE, limits = ct_order, name = "Cell type") +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(title = lab, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(color = plat_cols[[lab]], face = "bold", hjust = 0.5),
          axis.text = element_text(color = "black"),
          legend.key.size = unit(0.4, "cm"),
          legend.position = if (show_legend) "right" else "none")
}

# legend on atera panel: it has all types (incl. Low quality) so every rasterized
# key glyph has data and renders a colored swatch
pl <- lapply(order_plat, function(p) umap_plot(p, show_legend = (p == "atera")))
comb <- pl[[1]] / pl[[2]] / pl[[3]] + plot_layout(ncol = 1)
TARGET_H <- 12                                   # shared height across all 3 composites
nat_w <- 5.2; nat_h <- 13
save_pdf(comb, "Fig1F_umap_all3_vertical", nat_w * TARGET_H / nat_h, TARGET_H)
cat("Combined vertical UMAP written to", OUT, "\n")
