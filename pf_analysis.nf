#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
  getPed;
  getMap;
  getPedMap;
  get_ibd_segments;
  plot_ibd_segments;
  get_ibd_proportions;
  plot_ibd_proportions;
  get_ibd_iRstats;
  get_moimix_ped
} from "${projectDir}/modules/pf_analysis.mdl"

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

