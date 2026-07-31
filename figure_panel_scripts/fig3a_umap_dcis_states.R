#!/usr/bin/env Rscript
# UMAP of the reannotated epithelial object coloured by DCIS1-4 + Merged Normal.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table)
  library(ggrastr)
})

EP_RDS  <- "/path/to/data/atera_w_reannotated_Ep_states.rds"
MAP_CSV <- "/path/to/atera_demo/cell_types_mapping.csv"
OUTDIR  <- "results/infercnv_run/biological_insight"
OUT     <- file.path(OUTDIR, "umap_ep_dcis_states")

# phylo-tree palette (plot_phylo_tree_ggtree.R fixed_pal), R colour names -> used
# directly. NB: the tree maps Basal to skyblue1 (== Luminal DCIS1), so to keep the
# two distinguishable here Basal gets a distinct neutral grey instead.
STATE_PAL <- c(
  "Luminal DCIS1" = "skyblue1",
  "Luminal DCIS2" = "#CD2626",
  "Luminal DCIS3" = "#DAA520",
  "Luminal DCIS4" = "darkslategray3",
  "Merged Normal" = "#9ACD32",
  "Basal"         = "grey55"
)
order  <- names(STATE_PAL)
NORMAL <- c("Normal Luminal" = "Merged Normal", "Ductal Luminal" = "Merged Normal")

cat("reading", EP_RDS, "...\n"); t0 <- Sys.time()
obj <- readRDS(EP_RDS)
cat("loaded in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# locate a UMAP reduction
reds <- names(obj@reductions)
cat("reductions:", paste(reds, collapse = ", "), "\n")
umap_name <- reds[grepl("umap", reds, ignore.case = TRUE)][1]
if (is.na(umap_name)) stop("no UMAP reduction found in the Ep object")
emb <- as.data.frame(Embeddings(obj, umap_name))
colnames(emb)[1:2] <- c("UMAP1", "UMAP2")

# cell ids: prefer the cell_id meta column (matches the mapping), else rownames
md <- obj@meta.data
cid <- if ("cell_id" %in% colnames(md)) as.character(md$cell_id) else rownames(md)
emb$cell_id <- cid

# attach DCIS labels from the broad-object mapping, collapse normals -> Merged
# Normal and all basal states (Basal DCIS1/2/3, Normal Basal) -> Basal
mp <- fread(MAP_CSV)[, .(cell_id = as.character(cell_id), Cell_states = as.character(Cell_states))]
lab <- mp$Cell_states
lab[lab %in% names(NORMAL)] <- NORMAL[lab[lab %in% names(NORMAL)]]
lab[grepl("Basal", lab)] <- "Basal"
mp$state <- lab
df <- merge(emb, mp[, .(cell_id, state)], by = "cell_id", all.x = TRUE)

n_tot <- nrow(df)
df <- df[df$state %in% order, ]
df$state <- factor(df$state, levels = order)
cat(sprintf("%d of %d epithelial cells carry a DCIS1-4/Merged-Normal label\n",
            nrow(df), n_tot))
print(table(df$state))

# cache coords+labels so the plot can be re-rendered without reloading the RDS
CACHE <- file.path(OUTDIR, "umap_ep_dcis_states_coords.csv")
fwrite(df[, c("cell_id", "UMAP1", "UMAP2", "state")], CACHE)
cat("cached coords ->", CACHE, "\n")

# plot order: draw Basal/reference first (underneath), rare DCIS states last
draw_order <- c("Basal", "Luminal DCIS2", "Luminal DCIS3", "Luminal DCIS1",
                "Merged Normal", "Luminal DCIS4")
df <- df[order(match(df$state, draw_order)), ]

p <- ggplot(df, aes(UMAP1, UMAP2, colour = state)) +
  rasterise(geom_point(size = 0.35, stroke = 0), dpi = 400, dev = "ragg") +
  scale_colour_manual(values = STATE_PAL, name = "Cell state", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  coord_equal() +
  labs(title = "Epithelial UMAP - Luminal DCIS states + Merged Normal + Basal",
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(size = 12),
        axis.text = element_blank(), axis.ticks = element_blank())

ggsave(paste0(OUT, ".png"), p, width = 8, height = 6, dpi = 400)
# editable-text PDF: base pdf() device keeps text as real (editable) text in
# Illustrator; rasterise() already flattened the points to a raster layer
ggsave(paste0(OUT, ".pdf"), p, width = 8, height = 6, device = "pdf")
cat("saved", paste0(OUT, ".png"), "and .pdf\n")
