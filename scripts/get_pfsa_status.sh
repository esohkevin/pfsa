#!/bin/bash

vcf=/scratch/GeneMAP/genemap/shared/cameroon/malariagen/pf3D7_biallelic_snps.vcf.gz

bcftools \
  view \
  -r Pf3D7_02_v3:631190 \
  --threads 40 ${vcf} | \
bcftools \
  query \
  -f '[%SAMPLE %GT\n]' | \
awk '{if($2=="0/0" || $2 == "0|0") {$3="pfsa1-"} else if($2=="1/1" || $2 == "1|1") {$3="pfsa1+hom"} else if($2=="0/1" || $2 == "0|1" || $2=="1/0" || $2 == "1|0") {$3="pfsa1+het"} else {$3="NA"} {print $1,$2,$3}}' \
> pf7_pfsa1_Pf3D7_02_v3_631190_T_A_status.txt

grep -v -e 'NA' -e '+het' \
  pf7_pfsa1_Pf3D7_02_v3_631190_T_A_status.txt | \
  sed 's/hom//g' \
  > pf7_pfsa1_Pf3D7_02_v3_631190_T_A_status_clean.txt

bcftools \
  view \
  -r Pf3D7_02_v3:814288 \
  --threads 40 ${vcf} | \
bcftools \
  query \
  -f '[%SAMPLE %GT\n]' | \
awk '{if($2=="0/0" || $2 == "0|0") {$3="pfsa2-"} else if($2=="1/1" || $2 == "1|1") {$3="pfsa2+hom"} else if($2=="0/1" || $2 == "0|1" || $2=="1/0" || $2 == "1|0") {$3="pfsa2+het"} else {$3="NA"} {print $1,$2,$3}}' \
> pf7_pfsa2_Pf3D7_02_v3_814288_C_T_status.txt

grep -v -e 'NA' -e '+het' \
  pf7_pfsa2_Pf3D7_02_v3_814288_C_T_status.txt | \
  sed 's/hom//g' \
  > pf7_pfsa2_Pf3D7_02_v3_814288_C_T_status_clean.txt

bcftools \
  view \
  -r Pf3D7_11_v3:1058035 \
  --threads 40 ${vcf} | \
bcftools \
  query \
  -f '[%SAMPLE %GT\n]' | \
awk '{if($2=="0/0" || $2 == "0|0") {$3="pfsa3-"} else if($2=="1/1" || $2 == "1|1") {$3="pfsa3+hom"} else if($2=="0/1" || $2 == "0|1" || $2=="1/0" || $2 == "1|0") {$3="pfsa3+het"} else {$3="NA"} {print $1,$2,$3}}' \
> pf7_pfsa3_Pf3D7_11_v3_1058035_T_A_status.txt

grep -v -e 'NA' -e '+het' \
  pf7_pfsa3_Pf3D7_11_v3_1058035_T_A_status.txt | \
  sed 's/hom//g' \
  > pf7_pfsa3_Pf3D7_11_v3_1058035_T_A_status_clean.txt

bcftools \
  view \
  -r Pf3D7_04_v3:1121472 \
  --threads 40 ${vcf} | \
bcftools \
  query \
  -f '[%SAMPLE %GT\n]' | \
awk '{if($2=="0/0" || $2 == "0|0") {$3="pfsa4-"} else if($2=="1/1" || $2 == "1|1") {$3="pfsa4+hom"} else if($2=="0/1" || $2 == "0|1" || $2=="1/0" || $2 == "1|0") {$3="pfsa4+het"} else {$3="NA"} {print $1,$2,$3}}' \
> pf7_pfsa4_Pf3D7_04_v3_1121472_T_A_status.txt

grep -v -e 'NA' -e '+het' \
  pf7_pfsa4_Pf3D7_04_v3_1121472_T_A_status.txt | \
  sed 's/hom//g' \
  > pf7_pfsa4_Pf3D7_04_v3_1121472_T_A_status_clean.txt

