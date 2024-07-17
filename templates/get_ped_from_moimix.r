#!/usr/bin/env Rscript

require(SeqArray)
require(data.table)
require(isoRelate)
require(moimix)

vcf <- "${vcf}"
base <- "${vcf.simpleName}"
gds <- paste0(base, ".gds")
fws <- paste0(base, ".txt")

# convert VCF to GDS
seqVCF2GDS(
  vcf,
  gds
)

# load GDS
isolates <- seqOpen(
  gds
)

# load metadata and get MOI - NB: sample order is the same as in VCF/GDS file
gh.pheno <- read.table(
  fws, 
  h=T
)

gh.fws <- gh.pheno$Fws

# export plink PED/MAP file with heterozygote calls
extractPED(
  gdsfile = isolates, 
  moi.estimates = gh.fws, 
  use.hets = T, 
  outfile = base
)
