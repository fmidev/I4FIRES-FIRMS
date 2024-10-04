#!/bin/bash

set -e 
set -u

metadata=fire_metadata_ecodata_Akagi_continents_withAVB_v2.ini
mask=masks-5w/MOD61-2023-mask-5w.nc

#regions="Europe South_Asia SouthEast_Asia"

fires_arch=${fires_arch:-FIRMS-SOURCE}
regions="${regions:-Europe}"
gzip="${gzip:-gzip}"

iprogress=0
nprogress=`echo $regions | wc -w`

date=`date -u +%Y%m%d`
tmpdir=$fires_arch/tmp


mkdir -p $tmpdir

for reg in $regions; do
  URL="https://firms.modaps.eosdis.nasa.gov/data/active_fire/modis-c6.1/csv/MODIS_C6_1_${reg}_7d.csv"
  filebase=`basename $URL .csv`
  outd=$fires_arch/$filebase
  outf=$outd/${filebase}_${date}.csv.gz
  tmpf=$tmpdir/$filebase.csv
 # [ -f ${outf} ] && continue
  mkdir -p $outd
  [  -f ${outf} ] || ( wget $URL -O $tmpf && $gzip -c $tmpf > $outf.tmp && mv $outf.tmp $outf && rm $tmpf )


  firelistd=$fires_arch/$filebase-firelist
  mkdir -p $firelistd
  
  firelistbase=`basename $filebase _7d`
  ##Generate firelists
  catfiles=""
  firelistcmds=""
  for d in `seq -4 4`; do ## Should start at least 2 days after file start
    outd=`date -u  -d "$d days" +%Y-%m-%d`
    firelistnc=$firelistd/${firelistbase}_$outd.nc
    catfiles="$catfiles $firelistnc"
    firelistcmds="$firelistcmds $python3 FIRMS_cluster.py -v --metadata $fires_arch/$metadata --mask=$mask --extend=7 ${outf} --datepref=$outd  $firelistnc \n"
  done 
  
  ## generate firelists
  echo -e $firelistcmds | xargs -P 1 -I{} -t sh -c '{}'

  latestbase=${firelistbase}_latest

  echo Concatenating $reg to latest

  ncrcat -O $catfiles $fires_arch/$latestbase.nc

  cat > $fires_arch/$latestbase.src <<EOF
FIRE_SOURCE_V2
source_name = Firelist-${reg}_latest_$date
source_sector_name = fire_PM

fire_metadata_file = ^$metadata
mode_distribution_type = FIXED_DIAMETER   ! later also: GAMMA_FUNCTION

aerosol_mode = 1  0.01 0.7  0.17 mkm ! PM 1      mode_number Dmin, Dmax, Daver D_unit

fire_list_file = ^$latestbase.nc

END_FIRE_SOURCE_V2
EOF

 echo $reg done!
 iprogress=`expr $iprogress + 1`; $progresscmd `expr $iprogress \* 100 / $nprogress`
done

$progresscmd 100

#usage: FIRMS_cluster.py [-h] [--metadata METADATA] [--datepref DATEPREF]
#                        [--mask MASK] [--extend EXTEND] [--rmse-fit]
#                        [--verbose]
#                        FIRMS-list.csv [FIRMS-list.csv ...] SILAM-list

