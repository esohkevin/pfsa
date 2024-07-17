#!/usr/bin/env Rscript


#args <- commandArgs(TRUE)
#
#if(length(args) < 2) {
#  message("Usage: run_isorelate.r [pedfile] [mapfile]")
#  quit(save="no")
#} else {

  require(SeqArray)
  require(data.table)
  require(isoRelate)
  require(moimix)

  pedfile <- "${ped}" #"ghana_pfsa4_Pf3D7_04_v3_1090896-1148289.ped" 
  mapfile <- "${map}" #"ghana_pfsa4_Pf3D7_04_v3_1090896-1148289.map"

  pedbase <- sub(".ped", "", pedfile)  
  
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
    input.map.distance = "cM",
    reference.map.distance = "cM"
  )
  
  my_parameters <- getIBDparameters(
    ped.genotypes = my_genotypes, 
    number.cores = ${task.cpus}
  )
  
  my_ibd <- getIBDsegments(
    ped.genotypes = my_genotypes,
    parameters = my_parameters, 
    number.cores = ${task.cpus}, 
    minimum.snps = 20, 
    minimum.length.bp = 50000,
    error = 0.001
  )
  
  write.table(
    my_ibd, 
    paste0(pedbase, "_ibd.txt"),
    col.names = T,
    row.names = F, 
    quote = F, 
    sep = "\\t"
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

  my_iR <- getIBDiR(
    ped.genotypes = my_genotypes, 
    ibd.matrix = my_matrix, 
    groups = NULL
  )

  # plot IBD segments
  png(paste0(pedbase, "_ibdseg.png"))
  plotIBDsegments(
    ped.genotypes = my_genotypes,
    ibd.segments = my_ibd,
    interval = NULL,
    annotation.genes = NULL,
    annotation.genes.color = NULL,
    highlight.genes = NULL,
    highlight.genes.labels = FALSE,
    highlight.genes.color = NULL,
    highlight.genes.alpha = 0.1,
    segment.height = 0.6,
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
 
  # plot the proportion of pairs IBD
  png(paste0(pedbase, "ibd_prop.png"))
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

  png(paste0(pedbase, "ibd_prop_with-strat.png"))
  plotIBDproportions(ibd.proportions = my_proportion,
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
    subpop.facet = TRUE)
  dev.off()
#}
