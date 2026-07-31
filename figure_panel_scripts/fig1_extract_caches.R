#!/usr/bin/env Rscript
# Cache small extracts from each Seurat object for the platform-comparison panels.
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)

CACHE <- "/path/to/platform_comparison/results/manuscript_cache"
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

files <- c(
  atera   = "/path/to/platform_comparison/data/atera_Cell_type_and_state.rds",
  Xenium280 = "/path/to/DCIS_12_xenium/results/only_s1_bot/label_transfer/DCIS_s1_bot_final_annotated.rds",
  Xenium5K  = "/path/to/platform_comparison/data/Xenium5K_Cell_type_and_state.rds"
)
assays <- c(atera = "RNA", Xenium280 = "RNA", Xenium5K = "Xenium")

# optional: Rscript manuscript_extract_caches.R Xenium280  -> only redo that platform
only <- commandArgs(trailingOnly = TRUE)
if (length(only) > 0) files <- files[intersect(names(files), only)]

# Canonical cell-type harmonisation across platforms.
harmonize <- function(ct) {
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
    "Endothelial"="Endothelial","Enothelial"="Endothelial","Endothelium"="Endothelial"
  )
  unname(map[ct])
}

# 3 canonical markers per cell type (from celltype_marker_ridges.R).
markers <- list(
  "Basal"          = c("KRT5","KRT14","TP63"),
  "Normal Luminal" = c("KRT8","KRT18","FOXA1"),
  "Cancer Cells"   = c("ERBB2","EPCAM","MKI67"),
  "Myeloid"        = c("CD68","CD163","LYZ"),
  "Tcells"         = c("CD3D","CD3E","CD8A"),
  "Bcells"         = c("MS4A1","CD19","CD79A"),
  "PlasmaCells"    = c("MZB1","JCHAIN","IGHG1"),
  "FBs"            = c("COL1A1","PDGFRA","DCN"),
  "Pericytes"      = c("RGS5","ACTA2","MCAM"),
  "Endothelial"    = c("PECAM1","VWF","CDH5")
)
marker_genes <- unique(unlist(markers))
N_PER_CT <- 3000   # subsample cap per cell type for marker-expression cache

for (p in names(files)) {
  cat("\n================", p, "================\n"); flush.console()
  t0 <- Sys.time()
  obj <- readRDS(files[[p]])
  cat("loaded in", round(difftime(Sys.time(), t0, units="mins"),1), "min; cells:", ncol(obj), "\n"); flush.console()
  a <- assays[[p]]
  DefaultAssay(obj) <- a

  counts <- GetAssayData(obj, assay = a, layer = "counts")   # genes x cells
  md <- obj@meta.data

  tot <- Matrix::colSums(counts)
  ng  <- Matrix::colSums(counts > 0)
  ct_raw <- if ("cell_types" %in% colnames(md)) as.character(md$cell_types) else NA_character_
  cat("raw cell_types:\n"); print(sort(table(ct_raw), decreasing=TRUE))
  ct <- harmonize(ct_raw)

  # Spatial coords: x_centroid meta cols (atera/Xen280) or FOV image (Xen5K).
  if ("x_centroid" %in% colnames(md)) {
    gx <- md$x_centroid; gy <- md$y_centroid
  } else {
    tc <- tryCatch(GetTissueCoordinates(obj), error = function(e) NULL)
    gx <- rep(NA_real_, ncol(counts)); gy <- gx
    if (!is.null(tc)) {
      ids <- if ("cell" %in% colnames(tc)) as.character(tc$cell) else rownames(tc)
      mtc <- match(colnames(counts), ids)
      gx <- tc$x[mtc]; gy <- tc$y[mtc]
    }
  }

  meta <- data.frame(
    platform   = p,
    cell_id    = colnames(counts),
    x          = as.numeric(gx),
    y          = as.numeric(gy),
    n_count    = as.numeric(tot),
    n_gene     = as.numeric(ng),
    cell_type_raw = ct_raw,
    cell_type  = ct,
    stringsAsFactors = FALSE
  )
  saveRDS(meta, file.path(CACHE, paste0("meta_", p, ".rds")))

  # per-gene mean transcript per cell (all genes)
  gmean <- Matrix::rowMeans(counts)
  saveRDS(data.frame(gene = names(gmean), mean_count = as.numeric(gmean),
                     stringsAsFactors = FALSE),
          file.path(CACHE, paste0("genemean_", p, ".rds")))

  # marker expression (log-normalised "data" layer), subsampled per cell type
  dat <- GetAssayData(obj, assay = a, layer = "data")
  mg  <- intersect(marker_genes, rownames(dat))
  keep <- !is.na(ct)
  idx  <- which(keep)
  # subsample per harmonized cell type
  sub <- unlist(lapply(split(idx, ct[idx]),
                       function(i) if (length(i) > N_PER_CT) sample(i, N_PER_CT) else i),
                use.names = FALSE)
  em <- as.matrix(t(dat[mg, sub, drop = FALSE]))
  mexpr <- data.frame(cell_type = ct[sub], em, check.names = FALSE)
  saveRDS(mexpr, file.path(CACHE, paste0("markerexpr_", p, ".rds")))
  cat("cached: meta(", nrow(meta), ") genemean(", length(gmean),
      ") markerexpr(", nrow(mexpr), "x", length(mg), ")\n")

  rm(obj, counts, dat, meta, gmean, em, mexpr); gc()
}
cat("\nALL CACHES WRITTEN to", CACHE, "\n")
