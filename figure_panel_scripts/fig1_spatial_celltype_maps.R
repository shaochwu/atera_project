#!/usr/bin/env Rscript
# Cell-state / cell-type spatial maps for all 3 platforms (Fig 1D style).
suppressMessages({library(ggplot2); library(ggrastr); library(dplyr); library(RColorBrewer); library(patchwork)})

CACHE <- "/path/to/platform_comparison/results/manuscript_cache"
OUT   <- "/path/to/atera_demo/results/infercnv_run/manuscript_panels/spatial_cellstates_all3"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

plat_cols <- c(Atera = "#4FA381", "Xenium 280" = "#D95F52", "Xenium 5K" = "#E5A552")

# ---- load per-platform data.frames (cell_id, x, y, Cell_states, cell_types) ----
atera <- read.csv("/path/to/atera_demo/cell_annotations.csv",
                  stringsAsFactors = FALSE)
x280  <- readRDS(file.path(CACHE, "cellstate_Xenium280.rds"))
x5k   <- readRDS(file.path(CACHE, "cellstate_Xenium5K.rds"))

dats <- list(
  atera     = list(df = atera, lab = "Atera",      flip = FALSE, transpose = FALSE),
  Xenium280 = list(df = x280,  lab = "Xenium 280", flip = TRUE,  transpose = FALSE),
  Xenium5K  = list(df = x5k,   lab = "Xenium 5K",  flip = FALSE, transpose = TRUE)
)

# ---- global lineage-grouped color map (shared state -> shared color) ----
assign_group <- function(s) {
  s <- as.character(s)
  g <- rep("Other", length(s))
  g[grepl("Basal", s)] <- "Basal"                                   # before Tumor
  g[grepl("Luminal DCIS|DCIS[0-9]|IBC|Proliferating Tumor|Lumen Border|Ductal Luminal|Tumor", s)] <- "Tumor"
  g[grepl("Normal Luminal", s)] <- "NormalLuminal"
  g[grepl("FBs|Fibro|Smooth Muscle|Myofibro|Pericyte", s)] <- "Stroma"
  g[grepl("Endothelial|Enothelial|Lymphatic", s)] <- "Endothelial"
  g[grepl("TRM|TAM|MONO|DC|Myeloid|Macro|Mast|Neutro|Granulo", s)] <- "Myeloid"
  g[grepl("CD4|CD8|Treg|Tcell|T cell|_Tcells|NK|TCF7|ISG_T|prolif_T", s)] <- "Tcell"
  g[grepl("Bcell|B cell|PCs|Plasma", s)] <- "Bcell"
  g[grepl("Other|Low quality|Low counts|Unknown|Unassigned|Unannot", s, ignore.case = TRUE)] <- "Other"
  g
}
ramps <- list(
  Tumor         = c("#fcbba1","#99000d"),
  Basal         = c("#bcbddc","#3f007d"),
  NormalLuminal = c("#f768a1","#ae017e"),
  Stroma        = c("#fdd0a2","#7f2704"),
  Endothelial   = c("#66c2a4","#01665e"),
  Myeloid       = c("#c7e9c0","#00441b"),
  Tcell         = c("#c6dbef","#08306b"),
  Bcell         = c("#fee391","#cc9a06"),
  Other         = c("grey75","grey55")
)
group_order <- c("Tumor","Basal","NormalLuminal","Stroma","Endothelial",
                 "Myeloid","Tcell","Bcell","Other")
build_global_state_cols <- function(all_states) {
  all_states <- sort(unique(all_states[!is.na(all_states)]))
  grp <- assign_group(all_states)
  cols <- c(); ordered <- c()
  for (g in group_order) {
    st <- sort(all_states[grp == g])
    if (!length(st)) next
    pal <- if (length(st) == 1) ramps[[g]][2] else colorRampPalette(ramps[[g]])(length(st))
    cols <- c(cols, setNames(pal, st)); ordered <- c(ordered, st)
  }
  list(cols = cols, order = ordered)
}
all_states <- c(atera$Cell_states, x280$Cell_states, x5k$Cell_states)
gs <- build_global_state_cols(all_states)
state_cols <- gs$cols; state_order <- gs$order

# cell_types: harmonize each platform's labels to canonical names that match the
# user-specified palette keys (DCIS/IBC -> Cancer Cells; Enothelial -> Endothelial;
# Pericyte -> Pericytes; T cells -> Tcells; etc.). Low-quality/unannotated -> Other.
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
  out <- unname(map[ct]); out[is.na(out)] <- "Low quality"  # Others/Low counts/Unannot
  out
}
# user-specified cell_types palette (canonical keys) + grey "Low quality"
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

# ---- plotting coords per platform (matching Fig 1D), then a SHARED x/y window ----
for (p in names(dats)) {
  o <- dats[[p]]; df <- o$df
  if (o$transpose) { df$px <- df$y; df$py <- df$x } else { df$px <- df$x; df$py <- df$y }
  if (o$flip) { df$px <- -df$px; df$py <- -df$py }
  dats[[p]]$df <- df
}
# identical axis window for all platforms (same as Fig 1D): common x/y span
# (max across platforms), each tissue shifted into a common 0-based frame.
fin <- lapply(dats, function(o) { d <- o$df; d[is.finite(d$px) & is.finite(d$py), ] })
Sx  <- max(sapply(fin, function(d) diff(range(d$px))))
Sy  <- max(sapply(fin, function(d) diff(range(d$py))))
for (p in names(dats)) {
  d <- dats[[p]]$df; f <- fin[[p]]
  d$px <- d$px - mean(range(f$px)) + Sx / 2
  d$py <- d$py - mean(range(f$py)) + Sy / 2
  dats[[p]]$df <- d
}
MAX_IN  <- 5.6
upi     <- max(Sx, Sy) / MAX_IN
panel_w <- Sx / upi                       # same for every panel
panel_h <- Sy / upi

spatial_plot <- function(df, lab, colvar, cols, order, legend_title, ncol_leg,
                         drop_lv = TRUE) {
  d <- df[is.finite(df$px) & is.finite(df$py), ]
  lv <- if (drop_lv) intersect(order, unique(d[[colvar]])) else order
  d[[colvar]] <- factor(d[[colvar]], levels = lv)
  d <- d[order(d[[colvar]] != "Low quality"), ]     # draw Low quality underneath
  ncell <- nrow(d)
  ptsize <- if (ncell > 1e5) 0.4 else 0.6
  ggplot(d, aes(x = px, y = py, color = .data[[colvar]])) +
    rasterise(geom_point(size = ptsize, stroke = 0), dpi = 400, dev = "ragg") +
    scale_color_manual(values = cols, drop = drop_lv, limits = lv, name = legend_title) +
    coord_fixed(ratio = 1, xlim = c(0, Sx), ylim = c(0, Sy), expand = FALSE) +
    guides(color = guide_legend(ncol = ncol_leg, override.aes = list(size = 2.2, alpha = 1))) +
    labs(title = lab, x = "x", y = "y") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(color = plat_cols[[lab]], face = "bold", hjust = 0.5),
          axis.text = element_text(color = "black"),
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 8))
}

ct_plots <- list()
for (p in names(dats)) {
  o <- dats[[p]]; df <- o$df
  nstate <- length(unique(df$Cell_states[!is.na(df$Cell_states)]))
  leg_w  <- if (nstate > 16) 2.6 else 1.6         # legend width allowance (inches)
  pS <- spatial_plot(df, o$lab, "Cell_states", state_cols, state_order, "Cell state",
                     if (nstate > 16) 2 else 1)
  save_pdf(pS, paste0("Fig_spatial_", p, "_Cell_states"), panel_w + leg_w, panel_h + 1.0)

  # cell_types: harmonize to canonical (palette) names and drop low-quality/Other
  df$cell_types <- harmonize_ct(df$cell_types)
  dfT <- df[df$cell_types %in% ct_order, ]
  pT <- spatial_plot(dfT, o$lab, "cell_types", ct_cols, ct_order, "Cell type", 1)
  save_pdf(pT, paste0("Fig_spatial_", p, "_cell_types"), panel_w + 1.6, panel_h + 1.0)
  # shared-legend version (all types listed) for the combined figure
  ct_plots[[p]] <- spatial_plot(dfT, o$lab, "cell_types", ct_cols, ct_order, "Cell type",
                                1, drop_lv = FALSE)
  cat("wrote", p, "cell-state + cell-type maps\n")
}

# combined 3-panel cell-type figure, stacked VERTICALLY with ONE shared legend.
# Legend kept on the atera panel: it contains all types (incl. Low quality), so
# every rasterized key glyph has data to draw a colored swatch.
comb <- (ct_plots[["Xenium280"]] + theme(legend.position = "none")) /
        (ct_plots[["Xenium5K"]]  + theme(legend.position = "none")) /
        (ct_plots[["atera"]]     + theme(legend.position = "right")) +
        plot_layout(ncol = 1)
TARGET_H <- 12                                   # shared height across all 3 composites
nat_w <- panel_w + 1.8; nat_h <- 3 * panel_h + 1.6
save_pdf(comb, "Fig_spatial_all3_celltype_vertical", nat_w * TARGET_H / nat_h, TARGET_H)

cat("\nAll-3-platform cell-state/type spatial maps under", OUT, "\n")
