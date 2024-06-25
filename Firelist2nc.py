#!/usr/bin/env python3
# encoding: utf-8
# Copyleft Rostislav Kouznetsov 07.2022- 

descr="""Converts .csv firelists to .nc format
    assigning durations
"""

import pandas as pd
import netCDF4 as nc4
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import argparse


def fu_daylength(lat, utctimestamp): #daylength in seconds

    julday =  utctimestamp /(24*3600) % 365.2425   ## Approximate julian day
    # http://mathforum.org/library/drmath/view/56478.html
    #   D = daylength
    #   L = latitude
    #   J = day of the year
    #
    #   P = asin[.39795*cos(.2163108 + 2*atan{.9671396*tan[.00860(J-186)]})]
    #
    #                          _                                         _
    #                         / sin(0.8333*pi/180) + sin(L*pi/180)*sin(P) \
    #   D = 24 - (24/pi)*acos{  -----------------------------------------  }
    #                         \_          cos(L*pi/180)*cos(P)           _/
    #
    #Roux:  added check for saturation, asin removed from P: very small difference
    #
    ftmp = .39795*np.cos(.2163108 + 2.*(.00860*(julday-186.)))
    ftmp = (np.sin(0.8333/57.6) + np.sin(lat/57.6)*np.sin(ftmp))/(np.cos(lat/57.6)*np.cos(ftmp))
    ftmp = np.clip(ftmp, -1., 1.) #min(max(ftmp, -1.), 1.)

    return  24 * 3600 * ( 1. -  np.arccos(ftmp) / 3.1415927)


def checkdaylength():
    lat = 20.
    jd = np.arange(366) 
    dl = fu_daylength(lat, jd * (24*3600)) / 3600
    print(np.amin(dl), np.amax(dl))
    plt.plot(jd, dl)
    plt.grid()
    plt.show()
    

def trimPrecision(a, max_abs_err, max_rel_err):
   assert (a.dtype == np.float32)
  # return None
   if max_rel_err > 0:
       keepbits = int(np.ceil(np.log2(1./min(1,max_rel_err)))) ## -1 for rounding
       maskbits = max(20 -keepbits, 0)
       mask=(0xFFFFFFFF >> maskbits)<<maskbits
       b = a.view(dtype=np.int32)
       b &= mask
   
   if max_abs_err > 0:
       log2maxerr = np.log(max_abs_err)/np.log(2.)
       quantum = 2.**(np.floor(log2maxerr)+1.)

       a_maxtrimmed = quantum * 2**24 ## Avoid overflow
       idx = (np.abs(a) < a_maxtrimmed)
       a[idx] = quantum * np.around(a[idx] / quantum)
        



 


def FireList2NC(FireDF, outname, trimDayOnly=False, description=None):
        ##head = "FireStartZ,lat,lon,FRPeffMW,FRPmeanOBS,size,luidx,lutype,mask,obssats,daynight\n" 

    #
    # All fires are whole-day except for "day-only" fires
    # that start 2 hours after sunrize and end at sunset

    # Drop masked values, sort by starttime

    df = FireDF[FireDF['mask']==0].sort_values('FireStartZ')
    df = df.reset_index(drop=True)

    timestamps = np.array([int(dt.datetime.strptime( s, "%Y-%m-%dT%H:%M").replace(tzinfo=dt.timezone.utc).timestamp()) for s in df.FireStartZ]) ## Utc start of local day
    durations = np.zeros((len(timestamps),), dtype = np.int32 ) + 24*3600 ## 24 hours
    #

    ### LUTs
    luts = np.sort(df.lutype.unique())
    nluts_in = np.amax(df.luidx)
    nluts = len(luts)
    lu2idx = {}
    lutlen = 0
    for il, l in enumerate(luts):
        lu2idx[l] = il
        lutlen = max(lutlen, len(l))




    dflons = df.lon.to_numpy()
    dflats = df.lat.to_numpy()

    dffrps = df.FRPmeanOBS.to_numpy()
    
    # Sort starting times
    idxorder = np.argsort(timestamps)
    if not np.all(idxorder >= 0 ):
        import pdb; pdb.set_trace()

    with  nc4.Dataset(outname, "w", format="NETCDF4_CLASSIC") as dst:

      dst.note="All fires last for 24 hours from local midnight"
      if not description is None:
          dst.description=description

      dst.createDimension("time", None)

#      dst.createDimension("lut", nluts)
#      dst.createDimension("lut_str", lutlen)

#      ## Lut table
#      lutname = dst.createVariable('lut_name', 'c', ('lut', 'lut_str'))
#      for il, l in enumerate(luts):
#          lutname[il,:] =  nc4.stringtoarr(l, lutlen, dtype='S') 



      t = dst.createVariable('time', np.int32, ("time",), complevel = 5, zlib=True)
      t.standard_name = "time" 
      t.long_name = "Fire start time"
      t.units = "seconds since 1970-01-01 00:00:00 UTC" ;
      t.calendar = "standard" ;
      t.axis = "T" ;
      t[:] = timestamps[idxorder]

      dur = dst.createVariable('duration', np.int32, ("time",), complevel = 5, zlib=True)
      dur.long_name = "Fire duration"
      dur.units = 's'
      dur[:] = durations[idxorder]


      lo =  dst.createVariable('lon', np.float32, ("time",), complevel = 5, zlib=True)
      lo.standard_name = "longitude" ;
      lo.long_name = "longitude" ;
      lo.units = "degrees_east" ;
      lo[:] = df.lon.to_numpy()[idxorder]

      la =  dst.createVariable('lat', np.float32, ("time",), complevel = 5, zlib=True)
      la.standard_name = "latitude" ;
      la.long_name = "latitude" ;
      la.units = "degrees_north" ;
      la[:] = df.lat.to_numpy()[idxorder]

#      luidx = dst.createVariable('lutidx', np.int8, ("time",), complevel = 5, zlib=True)
#      luidx[:] = [ lu2idx[t] for t in df.lutype[idxorder] ]

      vals = np.float32(df.FRPeffPerFireMW.to_numpy()[idxorder])
      trimPrecision(vals, 0.001, 0.01) ## 1kW, 1%
      frp = dst.createVariable('FRPeffPerFire', np.float32, ("time",), complevel = 5, zlib=True)
      frp.long_name = 'Fire radiative power weighted with "_per_fire" diurnal cycle,  used for plume height'
      frp.units = 'MW'
      frp[:] = vals

      frp = dst.createVariable('FRPeffTotal', np.float32, ("time",), complevel = 5, zlib=True)
      vals = np.float32(df.FRPeffTotMW.to_numpy()[idxorder])
      trimPrecision(vals, 0.001, 0.01) ## 1kW, 1%
      frp.long_name = 'Fire radiative power weighted with "_total" diurnal cycle,  used for amount'
      frp.units = 'MW'
      frp[:] = vals

      vals = np.float32(df.FRPmeanOBS.to_numpy()[idxorder])
      trimPrecision(vals, 0.001, 0.01) ## 1kW, 1%
      frp = dst.createVariable('FRPmeanOBS', np.float32, ("time",), complevel = 5, zlib=True)
      frp.long_name = 'Fire radiative power mean of observations'
      frp.units = 'MW'
      frp[:] = vals

      vals = np.float32(df["size"].to_numpy()[idxorder])
      trimPrecision(vals, 0.01, 0.01) ## 10m, 1%
      fsize = dst.createVariable('firesize', np.float32, ("time",), complevel = 5, zlib=True)
      fsize.long_name = 'Fire size'
      fsize.units = 'km'
      fsize[:] = vals


def main():
    

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('--output', '-o', action="store", default='out.nc',  help='Output file')
    parser.add_argument('--trim_day', '-d', action="store_true", help='Trinm day-only fires to daytime')

    parser.add_argument("filenames", metavar='FireList(s)',  nargs='*', help='Input file(s)')

    #args = parser.parse_args("FIRES/TAN1/FireList-TAN1-202001.csv".split())
    args = parser.parse_args()



    if len(args.filenames) == 0:
        FireDF = pd.read_csv(sys.stdin)
    if len(args.filenames) == 1:
        FireDF = pd.read_csv(args.filenames[0])
    else:
        dflist = []
        for f in args.filenames:
            print("Reading ", f)
            df = pd.read_csv(f)
            dflist.append(df)


        FireDF = pd.concat(dflist,  ignore_index=True)


    FireList2NC(FireDF, args.output, trimDayOnly = args.trim_day, message=" ".join(sys.argv))


if __name__ == '__main__':
  # checkdaylength()
  main()
