#!/bin/bash

container=/scratch/eshkev001/containers/sickleinafrica-rbase-rehh-latest.img
scriptdir=/home/eshkev001/projects/mgen/scripts/
metadata=/home/eshkev001/projects/mgen/metadata/
pfile=/scratch/eshkev001/projects/mgen/plink/pf3D7_biallelic_snps
pf7meta=${metadata}Pf7_fws_samples_passqc_pfsa1_status_MOI_1_global_select.ped
pedtemp=${metadata}Pf7_fws_samples_passqc_pfsa1_status_MOI_1_global_select.ped
pedpop=${metadata}Pf7_fws_samples_passqc_MOI_select_global.pops

mkdir -p {ped,fws,bed,vcf}

plink2 \
  --allow-extra-chr \
  --geno 0.10 \
  --keep ${pedtemp} \
  --maf 0.01 \
  --make-bed \
  --mind 0.10 \
  --out bed/pf7_global \
  --pfile ${pfile}

join \
  -1 1 \
  -2 1 \
  <(sed 's/\t/ /g' bed/pf7_global.fam | cut -f1-2 -d' ' | sort -k1) \
  <(sort -k1 ${pedtemp} | cut -f2- -d' ') \
  > pf7_global.fam.country.txt

for pop in $(cat ${pedpop}) ; do 
  grep -w ${pop} pf7_global.fam.country.txt \
    > ped/Pf7_fws_samples_passqc_MOI_${pop}.ped

  popped="ped/Pf7_fws_samples_passqc_MOI_${pop}.ped"

  for sample in $(awk '{print $1}' ${popped}); do
    grep -w ${sample} ${pf7meta}
  done | \
  sed '1i FID IID PID MID MOI Fws' \
  > fws/Pf7_fws_samples_passqc_Fws_MOI_${pop}.txt

  popfws="fws/Pf7_fws_samples_passqc_Fws_MOI_${pop}.txt"

  for chr in {01..14}; do
    plink2 \
      --allow-extra-chr \
      --bfile pf7_global \
      --chr Pf3D7_${chr}_v3 \
      --geno 0.05 \
      --indiv-sort f ${popped} \
      --keep ${popped} \
      --maf 0.02 \
      --make-bed \
      --out bed/pf7_Pf3D7_${chr}_v3_${pop}
  
    plink2 \
      --allow-extra-chr \
      --bfile bed/pf7_Pf3D7_${chr}_v3_${pop} \
      --export vcf-4.2 bgz id-paste=iid \
      --out vcf/pf7_Pf3D7_${chr}_v3_${pop}

    popvcf="vcf/pf7_Pf3D7_${chr}_v3_${pop}.vcf.gz"

    singularity \
      exec ${container} \
      Rscript \
        ${scriptdir}get_moimix_ped.r \
          ${popfws} \
          ${popvcf}

    join -1 1 -2 1 \
      <(awk '{print $1,$4}' ${popped}) \
      <(cut -f2- -d' ' vcf/pf7_Pf3D7_${chr}_v3_${pop}.ped) | \
      awk '{iid=$1; $1=$2; $2=iid; print}' \
      > vcf/pf7_Pf3D7_${chr}_v3_${pop}.ped.temp;

    mv \
      vcf/pf7_Pf3D7_${chr}_v3_${pop}.ped.temp \
      vcf/pf7_Pf3D7_${chr}_v3_${pop}.ped;

    mv \
      vcf/pf7_Pf3D7_${chr}_v3_${pop}.{ped,map} \
      ../pedmap/
  done
done

