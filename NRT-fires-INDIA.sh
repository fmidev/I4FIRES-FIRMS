#!/bin/bash



metadata=fire_metadata_ecodata_Akagi_continents_withAVB_v2.ini
mask=masks-5w/MOD61-2023-mask-5w.nc



region=South_Asia



infiles="NRT-ASIA/MODIS_C6_1_South_Asia_7d.csv"

julday=date

for d in `seq -7 4`; do
  outd=`date -u  -d "$d days" +%Y-%m-%d`
  echo python3 FIRMS_cluster.py -v --metadata $metadata --mask=$mask --extend=7 $infiles --datepref=$outd  out.india-$outd.csv
done | xargs -P 20 -I{} sh -c '{}'

#usage: FIRMS_cluster.py [-h] [--metadata METADATA] [--datepref DATEPREF]
#                        [--mask MASK] [--extend EXTEND] [--rmse-fit]
#                        [--verbose]
#                        FIRMS-list.csv [FIRMS-list.csv ...] SILAM-list

