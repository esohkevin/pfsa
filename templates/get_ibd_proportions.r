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
  maf = ${params.maf},
  isolate.max.missing = ${params.mind},
  snp.max.missing = ${params.geno},
  chromosomes = NULL,
  input.map.distance = "M",
  reference.map.distance = "M"
)

# generate a binary IBD matrix
my_matrix <- getIBDmatrix(
  ped.genotypes = my_genotypes, 
  ibd.segments = my_ibd
)

# calculate the proportion of pairs IBD at each SNP
my_proportion <- getIBDproportion(
  ped.genotypes = my_genotypes, 
  ibd.matrix = my_matrix, 
  groups = NULL
)

write.table(
  my_proportion, 
  paste0(pedbase, "_ibd_prop.txt"),
  col.names = T, 
  row.names = F, 
  quote = F, 
  sep = "\\t"
)

# # plot the proportion of pairs IBD
# png(
#   paste0(pedbase, "_ibd_prop.png"),
#   height = 12,
#   width = 12,
#   units = "cm",
#   res = 300,
#   pointsize = 7
# )
#   plotIBDproportions(
#     ibd.proportions = my_proportion,
#     interval = NULL,
#     annotation.genes = NULL,
#     annotation.genes.color = NULL,
#     highlight.genes = NULL,
#     highlight.genes.labels = TRUE,
#     highlight.genes.color = NULL,
#     highlight.genes.alpha = 0.1,
#     add.rug = FALSE,
#     plot.title = "Proportion of pairs IBD",
#     add.legend = FALSE,
#     line.color = NULL,
#     facet.label = TRUE,
#     facet.scales = "fixed",
#     subpop.facet = FALSE
#   )
# dev.off()

my_groups <- ped[,c(1:2,1)]
colnames(my_groups) <- c("fid", "iid", "pid")

my_proportion <- getIBDproportion(
  ped.genotypes = my_genotypes,
  ibd.matrix = my_matrix,
  groups = my_groups
)

write.table(
  my_proportion,
  paste0(pedbase, "_ibd_prop_with-strat.txt"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\\t"
)

# png(
#   paste0(pedbase, "ibd_prop_with-strat.png"), 
#   height = 28,
#   width = 16,
#   units = "cm",
#   res = 600,
#   pointsize = 3
# )
#   plotIBDproportions(
#     ibd.proportions = my_proportion,
#     interval = NULL,
#     annotation.genes = NULL,
#     annotation.genes.color = NULL,
#     highlight.genes = NULL,
#     highlight.genes.labels = FALSE,
#     highlight.genes.color = NULL,
#     highlight.genes.alpha = 0.1,
#     line.color = NULL,
#     add.rug = FALSE,
#     plot.title = "Proportion of pairs IBD - with stratification",
#     add.legend = FALSE,
#     facet.label = TRUE,
#     facet.scales = "fixed",
#     subpop.facet = TRUE
#   )
# dev.off()

