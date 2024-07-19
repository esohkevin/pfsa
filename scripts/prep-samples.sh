#!/bin/bash

countries="Ghana Gambia Mali Cameroon DRCongo Kenya Tanzania Malawi Colombia Cambodia Vietnam PNG"

join -1 1 -2 1 <(sed '1d' Pf7_fws.txt | sort -k1) <(awk '$14 == "True"' Pf7_samples.txt | sort -k1) | sed '1i Sample Fws Study Country Admin_level_1 Country_latitude Country_longitude Admin_level_1_latitude Admin_level_1_longitude Year ENA All_samples_same_case Population %_callable QC_pass Exclusion_reason Sample_type Sample_was_in_Pf6' > Pf7_fws_samples_passqc.txt

join -1 1 -2 1 <(sed '1d' Pf7_fws_samples_passqc.txt | sort -k1) <(cut -f1,3 -d' ' pf7_pfsa1_Pf3D7_02_v3_631190_T_A_status_clean.txt | sort -k1) > temp1.txt

join -1 1 -2 1 <(sort -k1 temp1.txt) <(cut -f1,3 -d' ' pf7_pfsa2_Pf3D7_02_v3_814288_C_T_status_clean.txt | sort -k1) > temp2.txt 

join -1 1 -2 1 <(sort -k1 temp2.txt) <(cut -f1,3 -d' ' pf7_pfsa3_Pf3D7_11_v3_1058035_T_A_status_clean.txt | sort -k1) > temp3.txt 

join -1 1 -2 1 \
  <(sort -k1 temp3.txt) \
  <(cut -f1,3 -d' ' pf7_pfsa4_Pf3D7_04_v3_1121472_T_A_status_clean.txt | sort -k1) | \
sed '1i Sample Fws Study Country Admin_level_1 Country_latitude Country_longitude Admin_level_1_latitude Admin_level_1_longitude Year ENA All_samples_same_case Population %_callable QC_pass Exclusion_reason Sample_type Sample_was_in_Pf6 pfsa1_status pfsa2_status pfsa3_status pfsa4_status' \
> Pf7_fws_samples_passqc_pfsa_status_clean.txt

rm temp*

pfmeta=Pf7_fws_samples_passqc_pfsa_status_clean.txt

sed -i 's/Democratic_Republic_of_the_Congo/DRCongo/g' ${pfmeta}
sed -i "s|Côte_d'Ivoire|CotedIvoire|g" ${pfmeta}
sed -i 's/Papua_New_Guinea/PNG/g' ${pfmeta}
sed -i 's/Burkina_Faso/BurkinaFaso/g' ${pfmeta}

cut -f1-2,4,19- -d' ' ${pfmeta} | \
  awk '(NR==1){$8="MOI"}; (NR>1) { if($2<0.95) {$8=2} else{$8=1} } {print $0}' | \
  awk '{print $1,$2,$8,$3,$4,$5,$6,$7}' \
  > ${pfmeta/.txt/}_MOI.txt

for i in {5..8}; do
  pfstatus=$(cut -f1-4,${i} -d' ' ${pfmeta/.txt/}_MOI.txt | head -1 | cut -f5 -d' ')
  cut -f1-4,${i} -d' ' ${pfmeta/.txt/}_MOI.txt | \
  awk '{print $1,$1,$4,$5,$3,$2}' \
  > ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI.ped.txt
  sed '1d' ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI.ped.txt \
  > ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI.ped

  awk '!($5 == 2)' ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI.ped.txt \
  > ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI_1.ped.txt
  sed '1d' ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI_1.ped.txt \
  > ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI_1.ped

  for country in ${countries}; do
    grep \
      -w ${country} \
      ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI.ped;
  done | \
  awk '!($4 == "NA")' \
  > ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI_global_select.ped

  for country in ${countries}; do                                    
    grep \
      -w ${country} \
      ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI_1.ped; 
  done | \
  awk '!($4 == "NA")' \
  > ${pfmeta/_pfsa_status_clean.txt/}_${pfstatus}_MOI_1_global_select.ped
done

for country in ${countries}; do
  echo ${country}
done > Pf7_fws_samples_passqc_MOI_select_global.pops



# awk '{print $1,$1,$3,"0",$4,"0"}' Pf7_fws_samples_passqc_MOI.txt | sed '1d' | sed '1i FID IID PID MID MOI PHENO' > Pf7_fws_samples_passqc_MOI.ped.txt
# 
# awk '!($5 == 2) ' Pf7_fws_samples_passqc_MOI.ped.txt > Pf7_fws_samples_passqc_MOI_1.ped.txt
# 
# sed '1d' Pf7_fws_samples_passqc_MOI.ped.txt > Pf7_fws_samples_passqc_MOI.ped
# 
# sed '1d' Pf7_fws_samples_passqc_MOI_1.ped.txt > Pf7_fws_samples_passqc_MOI_1.ped
# 
# cut -f3 -d' ' Pf7_fws_samples_passqc_MOI.ped.txt | sed '1d' | sort | uniq > Pf7_fws_samples_passqc_MOI.ped.countries.txt
# 
# for country in ${countries}; do 
#   grep \
#     -w ${country} \
#     Pf7_fws_samples_passqc_MOI.ped; 
# done > Pf7_fws_samples_passqc_MOI_select_global.ped
# 
# cut \
#   -f3 \
#   -d' ' \
#   Pf7_fws_samples_passqc_MOI_select_global.ped | \
# sort | \
# uniq \
# > Pf7_fws_samples_passqc_MOI_select_global.pops
