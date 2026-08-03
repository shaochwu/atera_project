#!/usr/bin/env Rscript
# Myeloid cell-state UMAPs (DimPlot) for Atera / Xenium 280 / Xenium 5K.
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})

#Atera
Myeloid_atera = readRDS('/path/to/data/atera_myeloid_for_publication.rds')

colors = c('dodgerblue1', 'goldenrod2','red2',
           'deeppink3','seagreen2','darkorange', 
           'green4','orchid4','turquoise3')

DimPlot(Myeloid_atera, cols = colors) + coord_fixed()


#X280

X280_myeloid_clean2 = readRDS('/path/to/data/DCIS_s1_bot_Myeloid_annots.rds')

DimPlot(X280_myeloid_clean2,
        cols = c('green4', 'red2','orchid4', 'turquoise3', 'darkorange')) + coord_fixed()

#X5K

Myeloid_5K_clean = readRDS('/path/to/data/X5K_myeloid_cell_state_annotated.rds')

DimPlot(Myeloid_5K_clean,
        cols = c('turquoise3','gray', 'red2', 'green4','seagreen2','deeppink3', 'darkorange')) + coord_fixed()
