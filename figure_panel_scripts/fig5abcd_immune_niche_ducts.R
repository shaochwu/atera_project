#!/usr/bin/env Rscript
# Fig 5a-d -- immune-niche / clonal co-evolution across ducts (Atera). Output: fig5_abcd.pdf/.png

QUICK  <- FALSE
set.seed(42)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
  library(RANN)      # v2.6.2
  library(igraph)    # v2.3.2
})

## ---- parameters ------------------------------------------------------------
## The immune niche is niche 4 in the original numbering; in the published CSV
## it is relabelled "Immune niche". Both are recognised on load (see below).
IMMUNE_NICHE <- "Immune niche"
RADIUS       <- 200
KMAX         <- 2500
LINK_DIST    <- 40
KNN_LINK     <- 20
MIN_DUCT     <- 30
N_PERM       <- if (QUICK) 99 else 999

CLONES <- c("Normal Luminal", "Clone A", "Clone B", "Clone C", "Clone D")
XLAB   <- c("Normal", "Clone A", "Clone B", "Clone C", "Clone D")   # short axis labels
PAL    <- c("Normal Luminal" = "#8B5A2B",   # dark brown
            "Clone A"        = "#7B4FA3",   # dark purple
            "Clone B"        = "#D62728",   # red
            "Clone C"        = "#228B22",   # forest green
            "Clone D"        = "#FF8C00")   # orange

## ---- input -----------------------------------------------------------------
## Reads either the cached Seurat metadata (cell_metadata.rds) or the published
## all-cells annotation table (all_cells_annotations.csv). The CSV holds
## every cell's id, centroid, clonal/state label and niche -- the complete
## input needed to regenerate the figure, since the 200 um neighbourhoods span
## all cell types, not only the tumour compartment.
if (file.exists("cell_metadata.rds")) {
  md <- readRDS("cell_metadata.rds")
} else {
  md <- read.csv("all_cells_annotations.csv", check.names = FALSE)
  names(md)[names(md) == "cell_state_w_clone"] <- "cell_states_w_clones"
  names(md)[names(md) == "niche"]              <- "niches"
}
stopifnot(all(c("x_centroid","y_centroid","niches","cell_states_w_clones") %in% names(md)))
## normalise the immune-niche label: accept the original "4" or the CSV's "Immune niche"
md$niches <- ifelse(as.character(md$niches) == "4", "Immune niche", as.character(md$niches))
md$is_immune <- as.integer(md$niches == IMMUNE_NICHE)
stopifnot(sum(md$is_immune) > 0)
XY  <- as.matrix(md[, c("x_centroid","y_centroid")])
sel <- which(md$cell_states_w_clones %in% CLONES)

## ---- immune-niche fraction within RADIUS of each tumour cell ---------------
QXY <- XY[sel, , drop = FALSE]
nn  <- nn2(XY, QXY, k = KMAX, searchtype = "radius", radius = RADIUS)
tot <- rowSums(nn$nn.idx > 0) - 1L
stopifnot(max(tot) + 1L < KMAX)
hit <- matrix(0L, nrow(nn$nn.idx), ncol(nn$nn.idx))
present <- nn$nn.idx > 0
hit[present] <- md$is_immune[nn$nn.idx[present]]
immune_local <- (rowSums(hit) - md$is_immune[sel]) / pmax(tot, 1)

cells <- data.frame(
  x = md$x_centroid[sel], y = md$y_centroid[sel],
  clone = factor(as.character(md$cell_states_w_clones[sel]), levels = CLONES),
  immune = immune_local)
cells$stage <- as.integer(cells$clone) - 1L

## ---- duct segmentation -----------------------------------------------------
segment_ducts <- function(df, link, knn, min_n) {
  N <- nrow(df); q <- as.matrix(df[, c("x","y")])
  g <- nn2(q, k = knn)
  i <- rep(seq_len(N), knn - 1); j <- as.vector(g$nn.idx[, -1]); d <- as.vector(g$nn.dists[, -1])
  e <- d <= link
  graph <- add_edges(make_empty_graph(n = N, directed = FALSE), as.vector(rbind(i[e], j[e])))
  df$duct <- components(graph)$membership
  ducts <- df %>% group_by(duct) %>%
    summarise(n_cells = n(),
              purity  = max(table(clone)) / length(clone),   # before `clone` shadows
              clone   = factor(names(which.max(table(clone))), levels = CLONES),
              immune  = mean(immune), x = mean(x), y = mean(y), .groups = "drop") %>%
    filter(n_cells >= min_n) %>% mutate(stage = as.integer(clone) - 1L)
  list(cells = df, ducts = ducts)
}
seg   <- segment_ducts(cells, LINK_DIST, KNN_LINK, MIN_DUCT)
cells <- seg$cells; ducts <- seg$ducts
cat(sprintf("ducts %d | purity %.2f | cells %d\n",
            nrow(ducts), median(ducts$purity), sum(ducts$n_cells)))

## ---- association + torus-shift null (stage and x) --------------------------
rho_stage <- suppressWarnings(cor(ducts$stage, ducts$immune, method = "spearman"))
rho_x     <- suppressWarnings(cor(ducts$x,     ducts$immune, method = "spearman"))
xr <- range(XY[,1]); yr <- range(XY[,2]); W <- diff(xr); H <- diff(yr)
immune_xy <- XY[md$is_immune == 1, , drop = FALSE]
by_duct   <- split(seq_len(nrow(cells)), cells$duct)[as.character(ducts$duct)]
dmean <- function(v, idx) vapply(idx, function(k) mean(v[k]), numeric(1))
null <- vapply(seq_len(N_PERM), function(b) {
  s <- cbind(((immune_xy[,1] - xr[1] + runif(1,0,W)) %% W) + xr[1],
             ((immune_xy[,2] - yr[1] + runif(1,0,H)) %% H) + yr[1])
  f <- rowSums(nn2(s, QXY, k = KMAX, searchtype = "radius", radius = RADIUS)$nn.idx > 0) / pmax(tot,1)
  dm <- dmean(f, by_duct)
  suppressWarnings(c(stage = cor(ducts$stage, dm, method="spearman"),
                     x     = cor(ducts$x,     dm, method="spearman")))
}, numeric(2))
p_stage <- (1 + sum(null["stage",] >= rho_stage, na.rm=TRUE)) / (1 + N_PERM)
p_x     <- (1 + sum(null["x",]     >= rho_x,     na.rm=TRUE)) / (1 + N_PERM)
cat(sprintf("stage rho=%.3f P=%.4f | x rho=%.3f P=%.4f\n", rho_stage, p_stage, rho_x, p_x))

## ============================================================================
## FIGURE
## ============================================================================
keep <- cells %>% filter(duct %in% ducts$duct)
set.seed(7); keep$duct_f <- factor(keep$duct, levels = sample(unique(keep$duct)))
duct_cols <- rep(c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6",
                   "#bfef45","#fabed4","#469990","#dcbeff","#9A6324","#800000","#aaffc3",
                   "#808000","#000075","#a9a9a9","#ffd8b1","#00a5a5","#b03060"),
                 length.out = nrow(ducts))

## a: duct identity
pa <- ggplot(keep, aes(x, -y, colour = duct_f)) +
  geom_point(size = .3) + scale_colour_manual(values = duct_cols, guide = "none") +
  coord_fixed() + theme_void(base_size = 11)

## b: immune-niche fraction map, horizontal colourbar at bottom
pb <- ggplot(keep, aes(x, -y)) +
  geom_point(colour = "grey90", size = .25) +
  geom_point(data = ducts, aes(x, -y, size = n_cells, fill = immune),
             shape = 21, colour = "grey20") +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Immune niche fraction",
                       guide = guide_colourbar(title.position = "top", title.hjust = .5,
                                               barwidth = 9, barheight = .6)) +
  scale_size_area(max_size = 9, guide = "none") +
  coord_fixed() + theme_void(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_text(size = 9))

## c: immune-niche fraction vs clonal stage, linear fit, torus-shift P
pc <- ggplot(ducts, aes(stage, immune)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey20", fill = "grey80", linewidth = .6) +
  geom_point(aes(size = n_cells, fill = clone), shape = 21, colour = "grey25", alpha = .85) +
  scale_fill_manual(values = PAL, guide = "none", drop = FALSE) +
  scale_size_area(max_size = 6, name = "Cells per duct", breaks = c(1000, 2000)) +
  scale_x_continuous(breaks = 0:4, labels = XLAB, expand = expansion(add = .4)) +
  annotate("text", x = 4, y = 0, label = sprintf("p = %.3f", p_stage),
           hjust = 1, vjust = 0, size = 3.4) +
  labs(x = NULL, y = "Immune niche fraction within 200 um") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = c(.28, .9), legend.background = element_blank(),
        legend.title = element_text(size = 9), legend.key.size = unit(.5, "lines"))

## d: immune-niche fraction vs x position, loess fit
pd <- ggplot(ducts, aes(x, immune)) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.9, se = TRUE,
              colour = "grey20", fill = "grey80", linewidth = .6) +
  geom_point(aes(size = n_cells, fill = clone), shape = 21, colour = "grey25", alpha = .85) +
  scale_fill_manual(values = PAL, guide = "none", drop = FALSE) +
  scale_size_area(max_size = 6, guide = "none") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("  Spearman rho = %.2f, torus-shift P = %.3f", rho_x, p_x),
           hjust = 0, vjust = 1.6, size = 3.4) +
  labs(x = "x centroid (um)", y = "Immune niche fraction (200 um)") +
  theme_classic(base_size = 11)

fig <- pa + pb + pc + pd +
  plot_layout(design = "ACD\nBCD") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave("fig5_abcd.pdf", fig, width = 15, height = 8)
ggsave("fig5_abcd.png", fig, width = 15, height = 8, dpi = 200)
cat("wrote fig5_abcd.pdf / .png\n")
