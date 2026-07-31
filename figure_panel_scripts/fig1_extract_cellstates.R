#!/usr/bin/env Rscript
# Cache per-cell coords + Cell_states/cell_types from the Xenium 280 & 5K objects.
suppressMessages({library(Seurat); library(Matrix)})

CACHE <- "/path/to/platform_comparison/results/manuscript_cache"
files <- c(
  Xenium280 = "/path/to/DCIS_12_xenium/results/only_s1_bot/label_transfer/DCIS_s1_bot_final_annotated.rds",
  Xenium5K  = "/path/to/platform_comparison/data/Xenium5K_Cell_type_and_state.rds"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) files <- files[intersect(names(files), args)]

for (p in names(files)) {
  cat("\n====", p, "====\n"); flush.console()
  obj <- readRDS(files[[p]])
  md  <- obj@meta.data

  # cell-state column (name varies): prefer exact "Cell_states"
  state_col <- if ("Cell_states" %in% colnames(md)) "Cell_states" else
    grep("state", colnames(md), ignore.case = TRUE, value = TRUE)[1]
  type_col  <- if ("cell_types" %in% colnames(md)) "cell_types" else
    grep("cell_type", colnames(md), ignore.case = TRUE, value = TRUE)[1]
  cat("state col:", state_col, " | type col:", type_col, "\n")

  # spatial coords
  if ("x_centroid" %in% colnames(md)) {
    gx <- md$x_centroid; gy <- md$y_centroid
  } else {
    tc <- tryCatch(GetTissueCoordinates(obj), error = function(e) NULL)
    gx <- rep(NA_real_, nrow(md)); gy <- gx
    if (!is.null(tc)) {
      ids <- if ("cell" %in% colnames(tc)) as.character(tc$cell) else rownames(tc)
      m <- match(rownames(md), ids); gx <- tc$x[m]; gy <- tc$y[m]
    }
  }

  df <- data.frame(
    cell_id     = rownames(md),
    x           = as.numeric(gx),
    y           = as.numeric(gy),
    Cell_states = if (!is.na(state_col)) as.character(md[[state_col]]) else NA_character_,
    cell_types  = if (!is.na(type_col))  as.character(md[[type_col]])  else NA_character_,
    stringsAsFactors = FALSE
  )
  saveRDS(df, file.path(CACHE, paste0("cellstate_", p, ".rds")))
  cat("cached", nrow(df), "cells\n")
  cat("Cell_states:\n"); print(sort(table(df$Cell_states), decreasing = TRUE))
  cat("cell_types:\n");  print(sort(table(df$cell_types),  decreasing = TRUE))
  rm(obj, md); gc()
}
cat("\nDONE\n")
