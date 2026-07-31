#!/usr/bin/env Rscript
# Vertical 3-panel H&E figure (Xenium 280 / Xenium 5K / Atera).
suppressMessages({library(ggplot2); library(png); library(grid); library(patchwork)})

CACHE <- "/path/to/platform_comparison/results/manuscript_cache"
OUT   <- "/path/to/atera_demo/results/infercnv_run/manuscript_panels/figure1_panelG_HE"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

plat_labels <- c(atera = "Atera", xenium_280 = "Xenium 280", xenium_5k = "Xenium 5K")
plat_cols   <- c(Atera = "#4FA381", "Xenium 280" = "#D95F52", "Xenium 5K" = "#E5A552")
order_plat  <- c("xenium_280", "xenium_5k", "atera")

# ---- plotted cell coords per platform (SAME transforms as the spatial maps) ----
finite_xy <- function(px, py) { ok <- is.finite(px) & is.finite(py); list(px = px[ok], py = py[ok]) }
d280 <- readRDS(file.path(CACHE, "cellstate_Xenium280.rds"))
d5k  <- readRDS(file.path(CACHE, "cellstate_Xenium5K.rds"))
atr  <- read.csv("/path/to/atera_demo/cell_annotations.csv")
coords <- list(
  xenium_280 = finite_xy(-d280$x, -d280$y),   # both axes flipped
  xenium_5k  = finite_xy( d5k$y,   d5k$x),     # transposed to landscape
  atera      = finite_xy( atr$x,   atr$y)      # as-is
)
spanx <- sapply(coords, function(c) diff(range(c$px)))
spany <- sapply(coords, function(c) diff(range(c$py)))
Sx <- max(spanx); Sy <- max(spany)             # shared window (identical to spatial maps)

he_panel <- function(p) {
  img <- readPNG(file.path(CACHE, paste0("he_oriented_", p, ".png")))
  ih <- dim(img)[1]; iw <- dim(img)[2]; img_asp <- ih / iw
  lab <- plat_labels[[p]]
  # fit the image (preserving its aspect) inside this platform's cell footprint,
  # centered in the shared 0..Sx / 0..Sy window -> no distortion, true rel. size
  bw <- spanx[[p]]; bh <- spany[[p]]
  if (img_asp > bh / bw) { dh <- bh; dw <- bh / img_asp } else { dw <- bw; dh <- bw * img_asp }
  xr <- Sx / 2 + c(-1, 1) * dw / 2
  yr <- Sy / 2 + c(-1, 1) * dh / 2
  ggplot() +
    annotation_raster(as.raster(img), xmin = xr[1], xmax = xr[2],
                      ymin = yr[1], ymax = yr[2], interpolate = TRUE) +
    coord_fixed(ratio = 1, xlim = c(0, Sx), ylim = c(0, Sy), expand = FALSE) +
    labs(title = lab, x = "x", y = "y") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(color = plat_cols[[lab]], face = "bold", hjust = 0.5),
          axis.text = element_text(color = "black"))
}

plots <- lapply(order_plat, he_panel)
comb  <- plots[[1]] / plots[[2]] / plots[[3]] + plot_layout(ncol = 1)

# same page height as the spatial / UMAP composites; panels share scale so each
# occupies the full Sx x Sy window -> equal-size panels stacked vertically
TARGET_H <- 12
panel_asp <- Sy / Sx
nat_w <- Sx / Sx * 5.0 + 1.4                    # base width incl. axis/label margin
nat_h <- 3 * (5.0 * panel_asp) + 1.4
w <- nat_w * TARGET_H / nat_h
ggsave(file.path(OUT, "Fig_HE_all3_vertical.pdf"), comb,
       width = w, height = TARGET_H, device = "pdf", useDingbats = FALSE, dpi = 300)
ggsave(file.path(OUT, "Fig_HE_all3_vertical.png"), comb,
       width = w, height = TARGET_H, dpi = 200)
cat("H&E vertical figure written to", OUT, "\n")
