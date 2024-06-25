#!/bin/bash



metadata=fire_metadata_ecodata_Akagi_continents_withAVB_v2.ini

OUTD=Firelist-FIRMS-MASfit; addopt=""
#OUTD=Firelist-FIRMS-RMSE; addopt="--rmse-fit"

mkdir -p $OUTD

for YY in `seq 2002 2023`; do
  preYY=`expr $YY - 1`
  nexYY=`expr $YY + 1`
  infiles="MODIS_ALL_COUNTRIES/modis/modis_${preYY}_all_countries.csv.gz MODIS_ALL_COUNTRIES/modis/modis_${YY}_all_countries.csv.gz MODIS_ALL_COUNTRIES/modis/modis_${nexYY}_all_countries.csv.gz"

  if [ $YY -ge 2023 ]; then
    mask=masks-5w/MOD61-${YY}-mask-5w.nc
  elif [ $YY -ge 2018 ]; then
    mask=masks-5w/Mask_20-oper-${YY}.nc
  else
    mask=masks-5w/Mask_21-${YY}.nc
  fi

  echo python3 FIRMS_cluster.py -v --metadata=$metadata --mask=$mask --datepref=$YY- $addopt  $infiles $OUTD/Firmslist-masked-$YY.nc
  
done #|xargs -P 24 -I{} sh -c '{}'

#usage: FIRMS_cluster.py [-h] [--metadata METADATA] [--datepref DATEPREF]
#                        [--mask MASK] [--extend EXTEND] [--rmse-fit]
#                        [--verbose]
#                        FIRMS-list.csv [FIRMS-list.csv ...] SILAM-list

