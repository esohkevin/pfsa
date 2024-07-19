#!/usr/bin/env Rscript

require(SeqArray)
require(data.table)
require(isoRelate)
require(moimix)

pedfile <- "${ped}"
mapfile <- "${map}"
pedbase <- "${pedname}"  

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

my_parameters <- getIBDparameters(
  ped.genotypes = my_genotypes, 
  number.cores = ${task.cpus}
)

write.table(
  my_parameters,
  paste0(pedbase, "_ibd_parameters.txt"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\\t"
)

my_ibd <- getIBDsegments(
  ped.genotypes = my_genotypes,
  parameters = my_parameters, 
  number.cores = ${task.cpus}, 
  minimum.snps = ${params.min_snps_persegment}, 
  minimum.length.bp = ${params.min_ibdsegment_length},
  error = 0.001
)

write.table(
  my_ibd, 
  paste0(pedbase, "_ibd_segments.txt"),
  col.names = T,
  row.names = F, 
  quote = F, 
  sep = "\\t"
)

