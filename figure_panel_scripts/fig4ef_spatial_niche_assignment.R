#!/usr/bin/env Rscript
# Atera spatial niche assignment (BuildNicheAssay): niche map + cell-state x niche enrichment heatmap.

library(Seurat)
library(ggplot2)
library(pheatmap)

atera2 = readRDS('/path/to/data/atera_w_niche_assignment.rds')

cents <- CreateCentroids(atera2@meta.data[, c("x_centroid", "y_centroid")])
coords <- CreateFOV(cents, type = "centroids", assay = "RNA")
atera2[["fov"]] <- coords

table(atera2$cell_state)

# 1. compute niches using cell state assignments and DCIS clone assignment
atera2 <- BuildNicheAssay(object = atera2, fov = 'fov', group.by = "cell_states_w_clones",
                          niches.k = 5, neighbors.k = 30)

table(atera2$niches)
table(atera2$clone)
table(atera2$cell_states_w_clones)

# 2. plot sptail niche distribution
ggplot(atera2@meta.data, aes(x_centroid,y = 1- y_centroid, colour = as.character(niches))) + 
  geom_point(size = 0.01) + 
  scale_colour_manual(values = c('deeppink3','dodgerblue3','goldenrod1','seagreen','gray79')) +
  theme_void(base_size = 11) + coord_fixed() #+ NoLegend()


# 3. Cross-tabulate cell states vs niches
tab <- table(atera2$cell_states_w_clones, atera2$niches)

# 4. Normalize to column proportions (composition within each niche)
comp <- sweep(tab, 2, colSums(tab), "/")

# 5. Z-score across niches (per cell type) so "enrichment" pops out
mat <- t(scale(t(comp)))   # rows = cell types, scaled across niches

# 6. Heatmap
pheatmap(mat,
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("dodgerblue", "white", "firebrick1"))(100),
         main = "Cell type enrichment across niches",
         angle_col = 0,
         fontsize_row = 10)
