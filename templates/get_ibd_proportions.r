#!/usr/bin/env Rscript

require(SeqArray)
require(data.table)
require(isoRelate)
require(moimix)

pedfile <- "${ped}"
mapfile <- "${map}"
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

