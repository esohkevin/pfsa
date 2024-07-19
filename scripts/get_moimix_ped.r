#!/usr/bin/env Rscript


args <- commandArgs(TRUE)

if(length(args) < 2) {
  message("Usage: afr.r [fws-file] [vcf-file]")
  quit(save="n")
} else {
  require(SeqArray)
  require(data.table)
  require(isoRelate)
  require(moimix)

  f <- args[1]
  vcf <- args[2]
  outbase <- sub(".vcf.gz", "", vcf)
  
  pheno <- read.table(f, h=T)
  
  seqVCF2GDS(
    vcf,
    paste0(outbase, ".gds")
  )

  isolates <- seqOpen(
    paste0(outbase, ".gds")
  )

  extractPED(
    isolates, 
    moi.estimates = pheno$Fws, 
    use.hets = F, 
    outfile = outbase
  )
}

# seqVCF2GDS(
#   "pf7_Pf3D7_01_v3_afr_vcf.vcf.gz",
#   "pf7_Pf3D7_01_v3_afr_vcf.gds"
# )
# 
# seqVCF2GDS(
#   "pf7_Pf3D7_02_v3_afr_vcf.vcf.gz",
#   "pf7_Pf3D7_02_v3_afr_vcf.gds"
# )
# 
# seqVCF2GDS(
#   "pf7_Pf3D7_04_v3_afr_vcf.vcf.gz",
#   "pf7_Pf3D7_04_v3_afr_vcf.gds"
# )
# 
# seqVCF2GDS(
#   "pf7_Pf3D7_06_v3_afr_vcf.vcf.gz",
#   "pf7_Pf3D7_06_v3_afr_vcf.gds"
# )
# 
# seqVCF2GDS(
#   "pf7_Pf3D7_11_v3_afr_vcf.vcf.gz",
#   "pf7_Pf3D7_11_v3_afr_vcf.gds"
# )
# 
# 
# pf7_Pf3D7_01_v3_afr_vcf_isolates <- seqOpen("pf7_Pf3D7_01_v3_afr_vcf.gds")
# pf7_Pf3D7_02_v3_afr_vcf_isolates <- seqOpen("pf7_Pf3D7_02_v3_afr_vcf.gds")
# pf7_Pf3D7_04_v3_afr_vcf_isolates <- seqOpen("pf7_Pf3D7_04_v3_afr_vcf.gds")
# pf7_Pf3D7_06_v3_afr_vcf_isolates <- seqOpen("pf7_Pf3D7_06_v3_afr_vcf.gds")
# pf7_Pf3D7_11_v3_afr_vcf_isolates <- seqOpen("pf7_Pf3D7_11_v3_afr_vcf.gds")
# 
# extractPED(pf7_Pf3D7_01_v3_afr_vcf_isolates, moi.estimates=pheno$Fws, use.hets=T, outfile="pf7_Pf3D7_01_v3_afr")
# extractPED(pf7_Pf3D7_02_v3_afr_vcf_isolates, moi.estimates=pheno$Fws, use.hets=T, outfile="pf7_Pf3D7_02_v3_afr")
# extractPED(pf7_Pf3D7_04_v3_afr_vcf_isolates, moi.estimates=pheno$Fws, use.hets=T, outfile="pf7_Pf3D7_04_v3_afr")
# extractPED(pf7_Pf3D7_06_v3_afr_vcf_isolates, moi.estimates=pheno$Fws, use.hets=T, outfile="pf7_Pf3D7_06_v3_afr")
# extractPED(pf7_Pf3D7_11_v3_afr_vcf_isolates, moi.estimates=pheno$Fws, use.hets=T, outfile="pf7_Pf3D7_11_v3_afr")
