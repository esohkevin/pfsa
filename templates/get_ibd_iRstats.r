#!/usr/bin/env Rscript

require(SeqArray)
require(data.table)
require(isoRelate)
require(moimix)

pedfile <- "${ped}" #"ghana_pfsa4_Pf3D7_04_v3_1090896-1148289.ped" 
mapfile <- "${map}" #"ghana_pfsa4_Pf3D7_04_v3_1090896-1148289.map"
pedbase <- "${pedname}"
my_ibd <- read.table(
  "${ibdseg}",
  header = T
)

ped <- read.table(
  pedfile, 
  h = F
)

map <- read.table(
  mapfile, 
  h = F
)

ped.map <- list(
  ped, 
  map
)

my_genotypes <- getGenotypes(
  ped.map = ped.map,
  reference.ped.map = NULL,
  maf = 0.01,
  isolate.max.missing = 0.1,
  snp.max.missing = 0.1,
  chromosomes = NULL,
  input.map.distance = "M",
  reference.map.distance = "M"
)

# generate a binary IBD matrix
my_matrix <- getIBDmatrix(
  ped.genotypes = my_genotypes, 
  ibd.segments = my_ibd
)

my_iR <- getIBDiR(
  ped.genotypes = my_genotypes, 
  ibd.matrix = my_matrix, 
  groups = NULL
)

write.table(
  my_iR,
  paste0(pedbase, "_ibdiR.txt"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\\t"
)

png(
  paste0(pedbase, "_ibd_iR.png"),
  height = 12,
  width = 12,
  units = "cm",
  res = 300,
  pointsize = 7
)
plotIBDiR(
  ibd.iR = my_iR,
  interval = NULL,
  annotation.genes = NULL,
  annotation.genes.color = NULL,
  highlight.genes = NULL,
  highlight.genes.labels = FALSE,
  highlight.genes.color = NULL,
  highlight.genes.alpha = 0.1,
  point.size = 1,
  point.color = NULL,
  add.rug = FALSE,
  plot.title = "Significance of IBD sharing",
  add.legend = FALSE,
  facet.label = TRUE,
  facet.scales = "fixed"
)
dev.off()

# by group
my_groups <- ped[,c(1:2,1)]
colnames(my_groups) <- c("fid", "iid", "pid")

my_iR <- getIBDiR(
  ped.genotypes = my_genotypes,
  ibd.matrix = my_matrix,
  groups = my_groups
)

write.table(
  my_iR,
  paste0(pedbase, "_ibdiR_with-strat.txt"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\\t"
)

png(
  paste0(pedbase, "_ibd_iR_with-strat.png"),
  height = 28,
  width = 18,
  units = "cm",
  res = 300,
  pointsize = 7
)
plotIBDiR(
  ibd.iR = my_iR,
  interval = NULL,
  annotation.genes = NULL,
  annotation.genes.color = NULL,
  highlight.genes = NULL,
  highlight.genes.labels = FALSE,
  highlight.genes.color = NULL,
  highlight.genes.alpha = 0.1,
  point.size = 1,
  add.rug = FALSE,
  plot.title = "Significance of excess IBD sharing",
  add.legend = FALSE,
  facet.label = TRUE,
  facet.scales = "fixed"
)
dev.off()

