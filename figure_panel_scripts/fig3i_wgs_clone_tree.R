#!/usr/bin/env Rscript

# ggtree rendering of the WGS CNV phylogeny (16 samples; text-only branch CN-change labels).

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
tips_path   <- ga("--tips")
out_prefix  <- ga("--out")
outgroup    <- trimws(strsplit(ga("--outgroup", "cap_1_DL,cap_2_DL"), ",")[[1]])
fig_w       <- as.numeric(ga("--fig-width", "12"))
fig_h       <- as.numeric(ga("--fig-height", "8"))
aspect      <- as.numeric(ga("--aspect", "1.0"))   # coord_fixed ratio, as in example

tr <- read.tree(newick_path)
ntip <- length(tr$tip.label)

# ---- tip CNV-burden classes + palette --------------------------------------
tips <- as.data.frame(fread(tips_path)); rownames(tips) <- as.character(tips[[1]])
class_pal <- c("Diploid-like" = "#9ACD32", "Intermediate" = "#DAA520",
               "Tumour (high CNV)" = "#CD2626")
tip_cls <- function(t) {
  if (t %in% rownames(tips)) as.character(tips[t, "cnv_class"]) else "Intermediate"
}
tip_n <- function(t) if (t %in% rownames(tips)) as.integer(tips[t, "n_altered_arms"]) else NA

p <- ggtree(tr, size = 0.7, color = "grey30")

# reorder tips vertically so `group_top` sits at the TOP and the outgroup at the
# BOTTOM, keeping the default (ladderized) order for everything in between. This
# reassigns y directly (rotation only; topology, branch lengths, x unchanged).
group_top <- intersect(trimws(strsplit(
  ga("--group-top", "cap_16_DL,cap_11_DL,cap_12_DL,cap_14_DL"), ",")[[1]]),
  tr$tip.label)
if (length(group_top) > 0) {
  og  <- intersect(outgroup, tr$tip.label)
  d0  <- p$data[p$data$isTip, ]
  cur <- d0$label[order(d0$y)]                        # current bottom -> top
  mid <- setdiff(cur, c(group_top, og))
  bot2top <- c(og, mid, group_top)                    # outgroup bottom, group top
  rank <- setNames(seq_along(bot2top), bot2top)
  ch <- split(tr$edge[, 2], tr$edge[, 1])             # parent -> children
  desc_tips <- function(nd) if (nd <= ntip) nd else
    unlist(lapply(ch[[as.character(nd)]], desc_tips))
  newy <- p$data$y
  for (i in seq_len(nrow(p$data))) {
    nd <- p$data$node[i]
    newy[i] <- if (nd <= ntip) rank[[tr$tip.label[nd]]]
               else mean(rank[tr$tip.label[desc_tips(nd)]])
  }
  p$data$y <- newy
}
dat <- p$data

# tip labels: bold sample name (root gets a [root/normal] tag)
tip_text <- sapply(tr$tip.label, function(t) {
  disp <- gsub("_", " ", t)
  tag  <- if (t %in% outgroup) "  [root/normal]" else ""
  sprintf("%s%s", disp, tag)
})
tipdat <- data.frame(node = seq_len(ntip),
                     x = dat$x[match(seq_len(ntip), dat$node)],
                     y = dat$y[match(seq_len(ntip), dat$node)],
                     label = tip_text[tr$tip.label],
                     cls = vapply(tr$tip.label, tip_cls, character(1)))
tipdat$col <- class_pal[tipdat$cls]

# ---- branch CN-change labels (map clade -> node), arm-labelled --------------
ev <- as.data.frame(fread(events_path))
# arm remap, matching the example's --arm-labels: whole chromosomes -> "pq",
# arm-restricted regions -> single arm (calls derived from each region's span).
ARM <- c("11" = "11q", "12" = "12q", "13" = "13q", "14" = "14q", "15" = "15q",
         "4" = "4pq", "5" = "5pq", "17" = "17pq", "18" = "18pq", "20" = "20pq")
armlab  <- function(tok) { tok <- sub("^chr", "", tok)
                           if (tok %in% names(ARM)) ARM[[tok]] else tok }
chrnum  <- function(tok) as.integer(sub("[pq].*$", "", sub("^chr", "", tok)))
armrank <- function(tok) if (grepl("p$", tok)) 1L else if (grepl("q$", tok)) 2L else 3L
# collapse a same-direction {Np, Nq} pair into "Npq" (as the example shows whole-
# chromosome-concordant changes); a lone arm stays "Np"/"Nq".
collapse_pq <- function(labs) {
  if (!length(labs)) return(labs)
  base <- sub("[pq]$", "", labs); arm <- sub("^[0-9]+", "", labs)
  unlist(lapply(unique(base), function(b) {
    a <- arm[base == b]
    if (all(c("p", "q") %in% a)) paste0(b, "pq") else paste0(b, a)
  }))
}
fmt_tokens <- function(s) {
  if (is.na(s) || s == "") return(character(0))
  toks <- trimws(strsplit(s, ",")[[1]]); toks <- toks[nzchar(toks)]
  labs <- vapply(toks, armlab, character(1))
  labs <- labs[order(vapply(labs, chrnum, integer(1)), vapply(labs, armrank, integer(1)))]
  collapse_pq(labs)
}
node_of <- function(tips_str) {
  tt <- trimws(strsplit(tips_str, ",")[[1]])
  if (length(tt) == 1) which(tr$tip.label == tt) else getMRCA(tr, tt)
}
GAIN_COL <- "#b2182b"; LOSS_COL <- "#2166ac"
lab_off <- 0.34   # lift terminal-branch labels above the branch (as in example)
lab_rows <- list()
for (i in seq_len(nrow(ev))) {
  nd <- node_of(ev$descendant_tips[i])
  g <- fmt_tokens(ev$gains[i]); l <- fmt_tokens(ev$losses[i])
  gtxt <- if (length(g)) sprintf("<span style='color:%s'>%s</span>",
                                 GAIN_COL, paste0("+", g, collapse = " ")) else ""
  ltxt <- if (length(l)) sprintf("<span style='color:%s'>%s</span>",
                                 LOSS_COL, paste0("-", l, collapse = " ")) else ""
  txt <- paste(c(gtxt, ltxt)[nzchar(c(gtxt, ltxt))], collapse = "  ")
  off <- if (nd <= ntip) lab_off else 0
  lab_rows[[i]] <- data.frame(node = nd, x = dat$branch[dat$node == nd],
                              y = dat$y[dat$node == nd] + off, label = txt)
}
labdf <- do.call(rbind, lab_rows)

# ---- assemble plot (following the example's layer order & theme) ------------
xmax <- max(dat$x)
tip_nudge <- xmax * 0.07
cn_levels <- c("gain (+chr)", "loss (-chr)")
cn_cols   <- c("gain (+chr)" = GAIN_COL, "loss (-chr)" = LOSS_COL)
cn_df <- data.frame(x = rep(-1, 2), y = rep(-1, 2),
                    cn = factor(cn_levels, levels = cn_levels))
cls_levels <- names(class_pal)

p2 <- p +
  # solid tip marker (fill drives the burden legend), then bold tip-coloured label
  geom_point(data = tipdat, aes(x = x, y = y, fill = factor(cls, levels = cls_levels)),
             shape = 21, colour = "grey20", size = 3, inherit.aes = FALSE) +
  geom_text(data = tipdat, aes(x = x, y = y, label = label),
            colour = tipdat$col, fontface = "bold", hjust = 0,
            nudge_x = tip_nudge, size = 3.1, lineheight = 0.9, inherit.aes = FALSE) +
  geom_richtext(data = labdf, aes(x = x, y = y, label = label),
                fill = "white", label.colour = "grey60", label.size = 0.15,
                size = 2.5, label.padding = unit(0.1, "lines"), inherit.aes = FALSE) +
  scale_fill_manual(values = class_pal, name = "WGS CNV burden (tips)", drop = FALSE,
                    guide = guide_legend(override.aes = list(size = 3.5))) +
  # dummy off-canvas layer to build the branch CN-change colour legend
  geom_point(data = cn_df, aes(x = x, y = y, colour = cn), size = 0,
             inherit.aes = FALSE) +
  scale_colour_manual(values = cn_cols, name = "branch CN change",
                      guide = guide_legend(override.aes = list(size = 3.5))) +
  ggtitle(paste0("WGS CNV phylogeny (ggtree)\n",
                 "tips = WGS samples; branch labels = CN changes")) +
  xlim(0, xmax * 1.7) +
  coord_fixed(ratio = aspect, clip = "off") +
  theme_tree2() +
  xlab("NJ tree — WGS log2 on inferCNV gene-windows (distance from cap_1/cap_2 normal root)") +
  theme(legend.position = "left",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 11, margin = margin(t = 6)),
        axis.text.x = element_text(size = 9),
        axis.ticks.x = element_line(),
        plot.margin = margin(6, 6, 30, 6))

ggsave(paste0(out_prefix, ".png"), p2, width = fig_w, height = fig_h, dpi = 300)
ggsave(paste0(out_prefix, ".pdf"), p2, width = fig_w, height = fig_h)
cat("saved", paste0(out_prefix, ".png"), "\n")

# export the plotted tip order (top -> bottom) so heatmaps can match the tree
tipd <- dat[dat$isTip, ]
writeLines(tipd$label[order(-tipd$y)], paste0(out_prefix, "_tiporder.txt"))
