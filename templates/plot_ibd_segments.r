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

# plot IBD segments
png(
  paste0(pedbase, "_ibdseg.png"),
  height = 20,
  width = 16,
  units = "cm",
  res = 300,
  pointsize = 4
)
plotIBDsegments(
  ped.genotypes = my_genotypes,
  ibd.segments = my_ibd,
  interval = NULL,
  annotation.genes = NULL,
  annotation.genes.color = NULL,
  highlight.genes = NULL,
  highlight.genes.labels = FALSE,
  highlight.genes.color = NULL,
  highlight.genes.alpha = 0.15,
  segment.height = 0.5,
  number.per.page = NULL,
  fid.label = FALSE,
  iid.label = FALSE,
  ylabel.size = 9,
  add.rug = FALSE,
  plot.title = "Distribution of IBD segments",
  add.legend = TRUE,
  segment.color = NULL
)
dev.off()

plt_ibdseg_over_interval <- function(highlight.reg=NULL) {
  png(
    paste0(pedbase, "_ibdseg_over_interval.png"),
    height = 20,
    width = 16,
    units = "cm",
    res = 300,
    pointsize = 4
  )
  plotIBDsegments(
    ped.genotypes = my_genotypes,
    ibd.segments = my_ibd,
    interval = highlight.reg\$interval,
    annotation.genes = NULL,
    annotation.genes.color = NULL,
    highlight.genes = highlight.reg\$highlight.reg,
    highlight.genes.labels = TRUE,
    highlight.genes.color = "coral",
    highlight.genes.alpha = 0.2,
    segment.height = 0.5,
    number.per.page = NULL,
    fid.label = FALSE,
    iid.label = FALSE,
    ylabel.size = 9,
    add.rug = FALSE,
    plot.title = paste0(
      "Distribution of IBD segments over ", 
      highlight.reg\$interval[1],
      ":",
      highlight.reg\$interval[2],
      "-",
      highlight.reg\$interval[3]
    ),
    add.legend = TRUE,
    segment.color = NULL
  )
  dev.off()
}

chrom <- unique(my_ibd\$chr)

if(chrom == "Pf3D7_02_v3") {
  highlight.reg <- data.frame(
    name = c("ACS8", "PF3D7_0220300 (Plasmodium exported protein, unknown function)"),
    chr = c("Pf3D7_02_v3","Pf3D7_02_v3"),
    start = c("628091", "812892"),
    end = c("632681", "815853")
  )
  annotation <- list(
    highlight.reg = highlight.reg,
    interval = c("Pf3D7_02_v3", 628091-250000, 815853+250000)
  )
  plt_ibdseg_over_interval(
    highlight.reg = annotation
  )
} else if(chrom == "Pf3D7_04_v3") {
  highlight.reg <- data.frame(
    name = "FIKK4.2", 
    chr = "Pf3D7_04_v3", 
    start = "1115896", 
    end = "1123289"
  )
  annotation <- list(
    highlight.reg = highlight.reg,
    interval = c("Pf3D7_04_v3", 1121472-1000000, 1121472+100000)
  )
  plt_ibdseg_over_interval(
    highlight.reg = annotation
  )
} else if(chrom == "Pf3D7_11_v3") {
  highlight.reg <- data.frame(
    name = "PF3D7_1127000 (protein phosphatase, putative)",
    chr = "Pf3D7_11_v3",
    start = "1055701",
    end = "1058777"
  )  
  annotation <- list(
    highlight.reg = highlight.reg,
    interval = c("Pf3D7_11_v3", 1055701-250000, 1058777+250000)
  )
  plt_ibdseg_over_interval(
    highlight.reg = annotation
  )
}

