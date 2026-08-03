#!/usr/bin/env Rscript
# Combined myeloid marker DotPlot across Atera / Xenium 5K / Xenium 280.

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

## ---- Inputs ------------------------------------------------------------------
data_dir <- "/path/to/data"

files <- c(
  atera = file.path(data_dir, "atera_myeloid_noC3TRMs.rds"),
  X5K   = file.path(data_dir, "X5K_myeloid_cell_state_annotated.rds"),
  X280  = file.path(data_dir, "DCIS_s1_bot_Myeloid_annots.rds")
)

## Marker genes (x-axis order); grouped roughly by cell state.
genes <- c(
  ## monocytes
  "VCAN","FCN1","CCR2","S100A8","S100A9","C3","CX3CR1","FCGR3A",
  ## macrophages / tissue-resident macrophages
  "FOLR2","MRC1","PLTP","C1QC","CXCR4","CD163","COLEC12","SLC40A1",
  "LYVE1","SELENOP",
  ## SPP1 / TAM
  "SPP1","CHI3L1","APOC1","MMP9","MMP12","CHIT1","CTSK",
  ## CXCL9 TAM
  "CXCL9","CXCL11",
  ## cDC1
  "XCR1","CLEC9A","IDO1","CADM1",
  ## cDC2
  "CD1C","FCER1A","CLEC10A","CD1E",
  ## activation / mregDC
  "CD83","IL2RG","BIRC3","CCR7","LAMP3",
  ## pDC
  "CLEC4C","LILRA4","IL3RA","GZMB"
)

## Cell-state display order per dataset (top -> bottom within each block).
state_order <- list(
  atera = c("Classical MONOs","Non-classical MONOs","FOLR2 TRMs",
            "SPP1 TAMs","CXCL9 TAMs","cDC1s","cDC2s","mregDC","pDCs"),
  X5K   = c("Classical MONOs","FOLR2 TRMs","Unannot Macrophages",
            "SPP1 TAMs","CXCL9 TAMs","cDC1","cDC2"),
  X280  = c("Classical MONOs","MONOs","LYVE1 TRMs","TAMs","DCs")
)

## ---- 1. Extract per-cell-state expression statistics from each object --------
## Seurat::DotPlot computes, per gene per identity, using the "data" layer of the
## default assay: avg.exp = mean(expm1(x)) and pct.exp = 100 * mean(x > 0).
## Genes not present in an object's panel are queried out first (blank in plot).
extract_stats <- function(nm, path) {
  obj     <- readRDS(path)
  present <- intersect(genes, rownames(obj))
  df      <- DotPlot(obj, features = present)$data
  data.frame(
    dataset       = nm,
    id            = as.character(df$id),
    features.plot = as.character(df$features.plot),
    avg.exp       = df$avg.exp,
    pct.exp       = df$pct.exp,
    stringsAsFactors = FALSE
  )
}

dp <- do.call(rbind, Map(extract_stats, names(files), files))

## ---- 2. Scale mean expression per gene, WITHIN each dataset -----------------
## Matches Seurat's scaling: z-score across a dataset's cell states, then clamp
## to [-2.5, 2.5]. Genes with a single state / zero variance -> 0.
col.min <- -2.5; col.max <- 2.5
dp <- dp %>%
  group_by(dataset, features.plot) %>%
  mutate(avg.exp.scaled = {
    z <- as.numeric(scale(avg.exp)); z[is.na(z)] <- 0
    pmax(pmin(z, col.max), col.min)
  }) %>%
  ungroup()

## ---- 3. Axis ordering -------------------------------------------------------
## y-axis label: "<dataset>: <cell state>"
dp$group <- paste0(dp$dataset, ": ", dp$id)

row_order <- unlist(lapply(names(state_order),
                           function(d) paste0(d, ": ", state_order[[d]])))
stopifnot(setequal(row_order, unique(dp$group)))   # guard against label typos

dp$features.plot <- factor(dp$features.plot, levels = genes[genes %in% dp$features.plot])
dp$group         <- factor(dp$group, levels = rev(row_order))  # ggplot y is bottom-up

## faint separators between the three dataset blocks
blk    <- sub(":.*", "", rev(row_order))
sep_at <- which(head(blk, -1) != tail(blk, -1)) + 0.5

## ---- 4. Plot ----------------------------------------------------------------
p <- ggplot(dp, aes(features.plot, group)) +
  geom_hline(yintercept = sep_at, colour = "grey85", linewidth = 0.4) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
  scale_size(range = c(0, 6), name = "% expressed") +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                         name = "Scaled avg exp") +
  theme_cowplot() +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y  = element_text(size = 8),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA)) +
  xlab("") + ylab("") +
  ggtitle("Myeloid marker expression across datasets")

## ---- 5. Save ----------------------------------------------------------------
W <- max(8, 0.28 * nlevels(dp$features.plot))   # inches; wide enough that dots never overlap
H <- max(4, 0.28 * nlevels(dp$group))
ggsave(file.path(data_dir, "myeloid_combined_DotPlot.pdf"), p,
       width = W, height = H, limitsize = FALSE, bg = "white")
ggsave(file.path(data_dir, "myeloid_combined_DotPlot.png"), p, dpi = 300,
       width = W, height = H, limitsize = FALSE, bg = "white")
