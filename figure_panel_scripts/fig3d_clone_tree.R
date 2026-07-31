#!/usr/bin/env Rscript

# ggtree rendering of the inferCNV clone phylogeny (tip cell-state pies + branch CN-change labels).

suppressPackageStartupMessages({
  library(ggtree); library(treeio); library(ape)
  library(ggplot2); library(data.table); library(ggtext)
})

args <- commandArgs(trailingOnly = TRUE)
ga <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[[i + 1]]
}
newick_path <- ga("--newick")
events_path <- ga("--events")
comp_path   <- ga("--composition")
labels_path <- ga("--labels")
out_prefix  <- ga("--out")
root_tip    <- ga("--root-tip", "Merged_Normal")
pie_size    <- as.numeric(ga("--pie-size", "0.16"))   # inset width/height fraction
bar_w       <- as.numeric(ga("--bar-width", "0.9"))   # tip-barplot width  (x-units)
bar_h       <- as.numeric(ga("--bar-height", "0.8"))  # tip-barplot height (y-units)
bar_dx      <- as.numeric(ga("--bar-dx", "0"))        # extra rightward shift (x-units)
fig_w       <- as.numeric(ga("--fig-width", "12"))    # output figure width (inches)
fig_h       <- as.numeric(ga("--fig-height", "9"))    # output figure height (inches)
tip_glyph   <- ga("--tip-glyph", "bar")               # "bar" = composition bars; "none" = text only
lab_off     <- if (tip_glyph == "none") 0.34 else 0   # lift branch labels above the branch
# coord_fixed ratio (physical y-unit / x-unit). 1 = coord_equal (default). Values
# < 1 stretch the tree horizontally (wider); the tip barplot height is rescaled by
# 1/aspect so the bars keep the same physical shape.
aspect      <- as.numeric(ga("--aspect", "1"))
if (aspect != 1) bar_h <- bar_h / aspect

tr <- read.tree(newick_path)
ntip <- length(tr$tip.label)

# optional per-clone palette (group,color) shared with the heatmap so the clone
# colours match across figures; applied as direct tip-marker + tip-label colours
clone_pal_path <- ga("--clone-palette")
clone_pal <- NULL
if (!is.null(clone_pal_path) && file.exists(clone_pal_path)) {
  cp <- as.data.frame(fread(clone_pal_path)); colnames(cp) <- tolower(colnames(cp))
  clone_pal <- setNames(as.character(cp$color), as.character(cp$group))
}

lab <- as.data.frame(fread(labels_path)); rownames(lab) <- as.character(lab[[1]])

# composition only needed for the bar glyphs
comp <- NULL; states <- character(0)
if (tip_glyph != "none" && !is.null(comp_path) && file.exists(comp_path)) {
  comp <- as.data.frame(fread(comp_path), check.names = FALSE)
  rownames(comp) <- as.character(comp[[1]]); comp[[1]] <- NULL
  states <- colnames(comp)
}

# ---- state palette ---------------------------------------------------------
# fixed colours for the known cell states (skyblue1 = basal/non-DCIS, reused for
# the merged normal); anything else falls back to the generic pool
fixed_pal <- c("Luminal DCIS1" = "#87CEFF", "Luminal DCIS2" = "#DAA520",
               "Luminal DCIS3" = "#CD2626", "Luminal DCIS4" = "#7D26CD",
               "Merged Normal" = "#9ACD32", "Normal Luminal" = "#9ACD32",
               "Ductal Luminal" = "#9ACD32", "Basal" = "#87CEFF")
pal_pool <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
              "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD")
state_pal <- setNames(
  vapply(seq_along(states), function(i)
    if (states[i] %in% names(fixed_pal)) unname(fixed_pal[states[i]])
    else pal_pool[((i - 1) %% length(pal_pool)) + 1], character(1)),
  states)

# ---- base tree -------------------------------------------------------------
p <- ggtree(tr, size = 0.7, color = "grey30")
# optional clade rotations: --rotate "tipA,tipB;tipC,tipD" rotates the MRCA of
# each comma-separated tip group (used to reorder tips vertically)
rotate_arg <- ga("--rotate", "")
if (nzchar(rotate_arg)) {
  for (grp in strsplit(rotate_arg, ";")[[1]]) {
    tips <- trimws(strsplit(grp, ",")[[1]])
    nd <- if (length(tips) == 1) which(tr$tip.label == tips) else getMRCA(tr, tips)
    p <- ggtree::rotate(p, nd)
  }
}
# stretch tip spacing vertically (coord_equal locks the branch scale, so a taller
# canvas alone only adds whitespace; scaling y spreads the tips apart instead)
yscale <- as.numeric(ga("--y-scale", "1"))
if (yscale != 1) p$data$y <- p$data$y * yscale
dat <- p$data                      # node coordinates (x, y, branch, label, ...)

# tip labels: "clone X (n=...)" for clones, or the Cell_state name (underscores -> spaces)
tip_text <- sapply(tr$tip.label, function(t) {
  n <- if (t %in% rownames(lab)) as.integer(lab[t, "n_cells"])
       else if (!is.null(comp) && t %in% rownames(comp)) as.integer(sum(comp[t, ])) else NA
  disp <- gsub("_", " ", t)
  pre <- if (grepl("^[0-9_]+$", t)) "clone " else ""
  tag <- if (identical(t, root_tip)) "  [root]" else ""
  sprintf("%s%s%s\n(n=%s)", pre, disp, tag, format(n, big.mark = ","))
})
tipdat <- data.frame(node = seq_len(ntip),
                     x = dat$x[match(seq_len(ntip), dat$node)],
                     y = dat$y[match(seq_len(ntip), dat$node)],
                     label = tip_text[tr$tip.label])
# per-tip clone colour (from the shared palette); default grey when unmapped
tipdat$clonecol <- if (!is.null(clone_pal)) {
  ifelse(tr$tip.label %in% names(clone_pal),
         clone_pal[tr$tip.label], "grey30")
} else "black"

# ---- branch CN-change labels (map clade -> node) ---------------------------
ev <- as.data.frame(fread(events_path))
# --arm-labels: remap region tokens that are arm-restricted to carry their arm
# (e.g. 11 -> 11q, 12 -> 12q). Arm calls were derived from the genomic span of
# each region's inferCNV bins vs the hg38 centromere. Whole-chromosome events
# (span both p and q) keep their bare number. Off by default so other figures
# built from whole-chromosome tokens are unaffected.
arm_labels <- nzchar(ga("--arm-labels", ""))
# arm-restricted regions -> single arm; whole-chromosome regions (event spans
# both arms) -> "pq". Calls derived from each region's inferCNV-bin gene span
# vs the hg38 centromere.
ARM <- c("11" = "11q", "12" = "12q", "13" = "13q", "14" = "14q", "15" = "15q",
         "4" = "4pq", "5" = "5pq", "17" = "17pq", "18" = "18pq", "20" = "20pq")
armlab <- function(tok) {
  tok <- sub("^chr", "", tok)
  if (arm_labels && tok %in% names(ARM)) ARM[[tok]] else tok
}
chrnum <- function(tok) as.integer(sub("[pq].*$", "", sub("^chr", "", tok)))
armrank <- function(tok) if (grepl("p$", tok)) 1L else if (grepl("q$", tok)) 2L else 3L
# split a gains/losses field -> arm-labelled tokens sorted by chromosome number
# (then p before q, whole-chromosome last)
fmt_tokens <- function(s) {
  if (is.na(s) || s == "") return(character(0))
  toks <- trimws(strsplit(s, ",")[[1]]); toks <- toks[nzchar(toks)]
  labs <- vapply(toks, armlab, character(1))
  labs[order(vapply(labs, chrnum, integer(1)), vapply(labs, armrank, integer(1)))]
}
node_of <- function(tips_str) {
  tips <- trimws(strsplit(tips_str, ",")[[1]])
  if (length(tips) == 1) which(tr$tip.label == tips) else getMRCA(tr, tips)
}
# one label per branch (as before), but colour the chromosomes inline:
# gains red, losses blue (no purple) via ggtext richtext markup.
GAIN_COL <- "#b2182b"; LOSS_COL <- "#2166ac"
lab_rows <- list()
for (i in seq_len(nrow(ev))) {
  nd <- node_of(ev$descendant_tips[i])
  g <- fmt_tokens(ev$gains[i]); l <- fmt_tokens(ev$losses[i])
  gtxt <- if (length(g)) sprintf("<span style='color:%s'>%s</span>",
                                 GAIN_COL, paste0("+", g, collapse = " ")) else ""
  ltxt <- if (length(l)) sprintf("<span style='color:%s'>%s</span>",
                                 LOSS_COL, paste0("-", l, collapse = " ")) else ""
  txt <- paste(c(gtxt, ltxt)[nzchar(c(gtxt, ltxt))], collapse = "  ")
  # lift only TERMINAL-branch labels (to clear the tip text); keep internal-node
  # labels centred on their branch so clade events sit on the correct stem
  off <- if (nd <= ntip) lab_off else 0
  lab_rows[[i]] <- data.frame(node = nd, x = dat$branch[dat$node == nd],
                              y = dat$y[dat$node == nd] + off, label = txt)
}
labdf <- do.call(rbind, lab_rows)

# ---- assemble plot ---------------------------------------------------------
xmax <- max(dat$x)
draw_bars <- tip_glyph != "none" && !is.null(comp)
tip_nudge <- if (draw_bars) bar_w + xmax * 0.07 + bar_dx else xmax * 0.07
subtitle  <- if (draw_bars) {
  "tip barplots = cell-state composition; branch labels = CN changes"
} else {
  "tips = Cell_states; branch labels = CN changes"
}

p2 <- p +
  # solid clone-coloured marker at each tip (matches the heatmap clone palette)
  geom_point(data = tipdat, aes(x = x, y = y), colour = tipdat$clonecol,
             size = 3, inherit.aes = FALSE) +
  geom_text(data = tipdat, aes(x = x, y = y, label = label),
            colour = tipdat$clonecol, fontface = "bold",
            hjust = 0, nudge_x = tip_nudge, size = 3.1) +
  geom_richtext(data = labdf, aes(x = x, y = y, label = label),
                fill = "white", label.colour = "grey60", label.size = 0.15,
                size = 2.5, label.padding = unit(0.1, "lines")) +
  ggtitle(sprintf("inferCNV phylogeny (ggtree)\n%s", subtitle)) +
  xlim(0, xmax * 1.7 + bar_dx) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# dummy off-canvas layer to build the branch CN-change colour legend
cn_levels <- c("gain (+chr)", "loss (-chr)")
cn_cols   <- c("gain (+chr)" = "#b2182b", "loss (-chr)" = "#2166ac")
cn_df <- data.frame(x = rep(-1, 2), y = rep(-1, 2),
                    cn = factor(cn_levels, levels = cn_levels))

p3 <- p2
if (draw_bars) {
  # ---- tip barplots (stacked-by-state proportions) at the EXACT tip x,y -----
  keep <- tr$tip.label %in% rownames(comp)
  tipnodes <- which(keep)
  tipx <- dat$x[match(tipnodes, dat$node)]
  tipy <- dat$y[match(tipnodes, dat$node)]
  props <- as.matrix(comp[tr$tip.label[keep], states, drop = FALSE])
  props <- props / rowSums(props)
  G <- length(states); slot <- bar_w / G
  shift <- bar_w / 2 + xmax * 0.03 + bar_dx          # push bars right of branch labels
  x0 <- tipx + shift - bar_w / 2
  bars <- do.call(rbind, lapply(seq_along(tipnodes), function(i) {
    data.frame(state = factor(states, levels = states),
               xmin = x0[i] + (0:(G - 1)) * slot,
               xmax = x0[i] + (0:(G - 1)) * slot + slot * 0.85,
               ymin = tipy[i] - bar_h / 2,
               ymax = tipy[i] - bar_h / 2 + props[i, ] * bar_h)
  }))
  base_df <- data.frame(x = x0, xend = x0 + bar_w, y = tipy - bar_h / 2)
  p3 <- p3 +
    geom_segment(data = base_df, aes(x = x, xend = xend, y = y, yend = y),
                 colour = "grey70", linewidth = 0.25, inherit.aes = FALSE) +
    geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                               fill = state), inherit.aes = FALSE) +
    scale_fill_manual(values = state_pal, name = "cell state (tip barplots)",
                      guide = guide_legend(override.aes = list(size = 3.5)))
}

p3 <- p3 +
  geom_point(data = cn_df, aes(x = x, y = y, colour = cn), size = 0,
             inherit.aes = FALSE) +
  scale_colour_manual(values = cn_cols, name = "branch CN change",
                      guide = guide_legend(override.aes = list(size = 3.5))) +
  coord_fixed(ratio = aspect, clip = "off") +
  theme_tree2() +                    # adds the x-axis (distance from root)
  xlab(if (nzchar(ga("--xlab", ""))) ga("--xlab", "") else
       sprintf("CNV distance from %s root (Euclidean on inferCNV profile)",
               gsub("_", " ", root_tip))) +
  theme(legend.position = "left",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 11, margin = margin(t = 6)),
        axis.text.x = element_text(size = 9),
        axis.ticks.x = element_line(),
        plot.margin = margin(6, 6, 30, 6))

ggsave(paste0(out_prefix, ".png"), p3, width = fig_w, height = fig_h, dpi = 300)
ggsave(paste0(out_prefix, ".pdf"), p3, width = fig_w, height = fig_h)
cat("saved", paste0(out_prefix, ".png"), "\n")

# export the plotted tip order (top -> bottom) so heatmaps can match the tree
tipd <- dat[dat$isTip, ]
writeLines(tipd$label[order(-tipd$y)], paste0(out_prefix, "_tiporder.txt"))
