#!/usr/bin/env Rscript
# Fig 5a-d upstream: per-cell source table for fig5abcd_immune_niche_ducts.R. Output: fig5_cell_table.csv
suppressPackageStartupMessages({library(dplyr);library(RANN);library(igraph)})
IMMUNE_NICHE<-"Immune niche"; RADIUS<-200; KMAX<-2500; LINK<-40; KNN<-20; MIN<-30
CLONES<-c("Normal Luminal","Clone A","Clone B","Clone C","Clone D")
md<-readRDS("cell_metadata.rds")
md$niches<-ifelse(as.character(md$niches)=="4","Immune niche",as.character(md$niches)); md$is_immune<-as.integer(md$niches==IMMUNE_NICHE)
XY<-as.matrix(md[,c("x_centroid","y_centroid")])
sel<-which(md$cell_states_w_clones %in% CLONES)

## per-cell immune-niche fraction within 200um (index cell excluded both sides)
QXY<-XY[sel,,drop=FALSE]
nn<-nn2(XY,QXY,k=KMAX,searchtype="radius",radius=RADIUS)
tot<-rowSums(nn$nn.idx>0)-1L; stopifnot(max(tot)+1L<KMAX)
hit<-matrix(0L,nrow(nn$nn.idx),ncol(nn$nn.idx)); pr<-nn$nn.idx>0
hit[pr]<-md$is_immune[nn$nn.idx[pr]]
immune_local<-(rowSums(hit)-md$is_immune[sel])/pmax(tot,1)

## duct segmentation (single linkage)
q<-QXY; g<-nn2(q,k=KNN); N<-nrow(q)
i<-rep(seq_len(N),KNN-1); j<-as.vector(g$nn.idx[,-1]); d<-as.vector(g$nn.dists[,-1]); e<-d<=LINK
graph<-add_edges(make_empty_graph(n=N,directed=FALSE),as.vector(rbind(i[e],j[e])))
comp<-components(graph)$membership
tab<-data.frame(
  cell_id               = rownames(md)[sel],
  x_centroid            = md$x_centroid[sel],
  y_centroid            = md$y_centroid[sel],
  cell_state_w_clone    = as.character(md$cell_states_w_clones[sel]),
  niche                 = md$niches[sel],
  clonal_stage          = match(as.character(md$cell_states_w_clones[sel]),CLONES)-1L,
  n_neighbours_200um    = tot,
  immune_fraction_200um = round(immune_local,5),
  duct_id_raw           = comp
)
## keep only cells in retained ducts (>=30 cells); assign majority-clone label + duct means
sz<-table(comp); keepd<-as.integer(names(sz)[sz>=MIN])
tab$in_duct<-tab$duct_id_raw %in% keepd
dl<-tab %>% filter(in_duct) %>% group_by(duct_id_raw) %>%
  summarise(duct_n_cells=n(),
            duct_majority_clone=names(which.max(table(cell_state_w_clone))),
            duct_purity=round(max(table(cell_state_w_clone))/n(),3),
            duct_immune_fraction=round(mean(immune_fraction_200um),5),.groups="drop")
tab<-tab %>% left_join(dl,by="duct_id_raw")
## renumber ducts 1..128 for readability; NA for cells not in a retained duct
tab$duct_id<-match(tab$duct_id_raw,sort(keepd)); tab$duct_id_raw<-NULL
tab<-tab[order(tab$duct_id, tab$cell_id),
         c("cell_id","x_centroid","y_centroid","cell_state_w_clone","niche","clonal_stage",
           "n_neighbours_200um","immune_fraction_200um","in_duct","duct_id",
           "duct_majority_clone","duct_purity","duct_n_cells","duct_immune_fraction")]
write.csv(tab,"fig5_cell_table.csv",row.names=FALSE)
cat(sprintf("rows %d | in ducts %d | ducts %d\n",nrow(tab),sum(tab$in_duct),max(tab$duct_id,na.rm=TRUE)))
cat("cols:",paste(names(tab),collapse=", "),"\n")
