#!/bin/bash



metadata=fire_metadata_ecodata_Akagi_continents_withAVB_v2.ini
mask=masks-5w/MOD61-2023-mask-5w.nc


infiles="NRT-ASIA/MODIS_C6_1_SouthEast_Asia_MCD14DL_NRT_202417?.txt"


python3 FIRMS_cluster.py -v --metadata $metadata --mask=masks-5w/MOD61-2023-mask-5w.nc --extend=5  NRT-ASIA/MODIS_C6_1_South_Asia_MCD14DL_NRT_202417?.txt out.india-ext5.csv

#usage: FIRMS_cluster.py [-h] [--metadata METADATA] [--datepref DATEPREF]
#                        [--mask MASK] [--extend EXTEND] [--rmse-fit]
#                        [--verbose]
#                        FIRMS-list.csv [FIRMS-list.csv ...] SILAM-list

