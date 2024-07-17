#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {

  println "\nPf3D7 Popgen Analysis\n"

  if(params.pedmap_dir == 'NULL') {
    ped = getPed()
    map = getMap()
    map.combine(ped)
       .map { map, ped -> tuple(map.simpleName, map, ped) }
       .set { pedmap }
  }
  else {
    getPedMap()
      .map { name, file -> tuple(name, file.first(), file.last()) }
      .set { pedmap }
     
  }
  
  get_ibd_segments(pedmap)
    .map { pedname, ibdparam, seg -> tuple(pedname, seg) }
    .set { ibdseg }

  pedmap
    .join(ibdseg)
    .set { ibdprop_input }

  plot_ibd_segments(ibdprop_input)

  get_ibd_proportions(ibdprop_input)
    .set { ibdprop }

  ibdprop
    .join(ibdseg)
    .set { plotprops_input }

  plot_ibd_proportions(plotprops_input)

  ibdiR = get_ibd_iRstats(ibdprop_input)
}


def getPedMap() {
  return channel.fromFilePairs( params.pedmap_dir + "/*.{ped,map}", size: 2 )
}

def getPed() {
  return channel.fromPath( params.ped )
}

def getMap() {
  return channel.fromPath( params.map )
}

process get_ibd_segments() {
  tag "calculating stats ..."
  label 'isorelate'
  label 'rbase'
  publishDir \
    path: "${params.output_dir}/results/", \
    mode: 'copy'
  input:
    tuple \
      val(pedname), \
      path(map), \
      path(ped)
  output:
    tuple \
      val(pedname), \
      path("${pedname}_ibd_parameters.txt"), \
      path("${pedname}_ibd_segments.txt")
  script:
    template 'get_ibd_segments.r'
}

process get_ibd_proportions() {
  tag "calculating stats ..."
  label 'isorelate'
  label 'rbase'
  publishDir \
    path: "${params.output_dir}/results/", \
    mode: 'copy'
  input:
    tuple \
      val(pedname), \
      path(map), \
      path(ped), \
      path(ibdseg)
  output:
    tuple \
      val(pedname), \
      path("${pedname}_ibd_prop.txt"), \
      path("${pedname}_ibd_prop_with-strat.txt")
  script:
    template 'get_ibd_proportions.r'
}

process get_ibd_iRstats() {
  tag "calculating stats ..."
  label 'isorelate'
  label 'rbase'
  publishDir \
    pattern: "*.txt", \
    path: "${params.output_dir}/results/", \
    mode: 'copy'
  publishDir \
    pattern: "*.png", \
    path: "${params.output_dir}/results/figs/", \
    mode: 'copy'
  input:
    tuple \
      val(pedname), \
      path(map), \
      path(ped), \
      path(ibdseg)
  output:
    path("*.{png,txt}")
  script:
    template 'get_ibd_iRstats.r'
}

process plot_ibd_segments() {
  tag "calculating stats ..."
  label 'isorelate'
  label 'rbase'
  publishDir \
    path: "${params.output_dir}/results/figs/", \
    mode: 'copy'
  input:
    tuple \
      val(pedname), \
      path(map), \
      path(ped), \
      path(ibdseg)
  output:
    path("*.png")
  script:
    template 'plot_ibd_segments.r'
}

process plot_ibd_proportions() {
  tag "calculating stats ..."
  label 'isorelate'
  label 'rbase'
  publishDir \
    path: "${params.output_dir}/results/figs/", \
    mode: 'copy'
  input:
    tuple \
      val(pedname), \
      path(ibdprop), \
      path(ibdpropwithstrat), \
      path(ibdseg)
  output:
    path("*.png")
  script:
    template 'plot_ibd_proportions.r'
}

process get_moimix_ped() {
  tag "calculating stats ..."
  label 'isorelate'
  label 'rbase'
  publishDir \
    path: "${params.output_dir}/results/", \
    mode: 'copy'
  input:
    path(vcf)
  output:
    path("${vcf.simpleName}.bed")
  script:
    template 'get_ped_from_moimix.r'
}

