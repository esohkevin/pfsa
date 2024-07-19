#!/usr/bin/env Rscript

require(SeqArray)
require(data.table)
require(isoRelate)
require(moimix)

pedbase <- "${pedname}"

my_ibd <- read.table(
  "${ibdseg}",
  header = T
)

pops <- unique(my_ibd\$fid1)

if(length(pops) > 1) {
  # plot the proportion of pairs IBD: with stratification

  my_proportion <- read.table(
    "${ibdpropwithstrat}",
    header = T
  )

  png(
    paste0(pedbase, "_ibd_prop_with-strat.png"),
    height = 28,
    width = 16,
    units = "cm",
    res = 600,
    pointsize = 3
  )
    plotIBDproportions(
      ibd.proportions = my_proportion,
      interval = NULL,
      annotation.genes = NULL,
      annotation.genes.color = NULL,
      highlight.genes = NULL,
      highlight.genes.labels = FALSE,
      highlight.genes.color = NULL,
      highlight.genes.alpha = 0.1,
      line.color = NULL,
      add.rug = FALSE,
      plot.title = "Proportion of pairs IBD - with stratification",
      add.legend = FALSE,
      facet.label = TRUE,
      facet.scales = "fixed",
      subpop.facet = TRUE
    )
  dev.off()
} else {
  # plot the proportion of pairs IBD

  my_proportion <- read.table(
    "${ibdprop}",
    header = T
  )

  png(
    paste0(pedbase, "_ibd_prop.png"),
    height = 12,
    width = 12,
    units = "cm",
    res = 300,
    pointsize = 5
  )
    plotIBDproportions(
      ibd.proportions = my_proportion,
      interval = NULL,
      annotation.genes = NULL,
      annotation.genes.color = NULL,
      highlight.genes = NULL,
      highlight.genes.labels = TRUE,
      highlight.genes.color = NULL,
      highlight.genes.alpha = 0.1,
      add.rug = FALSE,
      plot.title = "Proportion of pairs IBD",
      add.legend = FALSE,
      line.color = NULL,
      facet.label = TRUE,
      facet.scales = "fixed",
      subpop.facet = FALSE
    )
  dev.off()
}

plt_ibdprop_over_interval <- function(highlight.reg=NULL) {
  png(
    paste0(pedbase, "_ibdprop_over_interval.png"),
    height = 16,
    width = 16,
    units = "cm",
    res = 300,
    pointsize = 1
  )
  plotIBDproportions(
    ibd.proportions = my_proportion,
    interval = highlight.reg\$interval,
    annotation.genes = NULL,
    annotation.genes.color = NULL,
    highlight.genes = highlight.reg\$highlight.reg,
    highlight.genes.labels = FALSE,
    highlight.genes.color = "coral",
    highlight.genes.alpha = 0.4,
    add.rug = FALSE,
    add.legend = FALSE,
    line.color = NULL,
    facet.label = TRUE,
    facet.scales = "fixed",
    subpop.facet = TRUE,
    plot.title = paste0(
      "Proportion of pairs IBD over ", 
      highlight.reg\$interval[1],
      ":",
      highlight.reg\$interval[2],
      "-",
      highlight.reg\$interval[3]
    )
  )
  dev.off()
}

chrom <- unique(my_ibd\$chr)

if(chrom == "Pf3D7_02_v3") { #pfsa1 and pfsa2
  highlight.reg <- data.frame(
    name = c("ACS8", "PF3D7_0220300 (Plasmodium exported protein, unknown function)"),
    chr = c("Pf3D7_02_v3","Pf3D7_02_v3"),
    start = c("628091", "812892"),
    end = c("632681", "815853")
  )
  annotation <- list(
    highlight.reg = highlight.reg,
    interval = c("Pf3D7_02_v3", 628091-200000, 815853+200000)
  )
  plt_ibdprop_over_interval(
    highlight.reg = annotation
  )
} else if(chrom == "Pf3D7_04_v3") { #pfsa4
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
  plt_ibdprop_over_interval(
    highlight.reg = annotation
  )
} else if(chrom == "Pf3D7_11_v3") { # pfsa3
  highlight.reg <- data.frame(
    name = "PF3D7_1127000 (protein phosphatase, putative)",
    chr = "Pf3D7_11_v3",
    start = "1055701",
    end = "1058777"
  )  
  annotation <- list(
    highlight.reg = highlight.reg,
    interval = c("Pf3D7_11_v3", 1055701-200000, 1058777+200000)
  )
  plt_ibdprop_over_interval(
    highlight.reg = annotation
  )
}

