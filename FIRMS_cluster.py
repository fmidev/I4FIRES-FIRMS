#!/usr/bin/env python3
# encoding: utf-8


#
#   Aggregate pixels to clusters for each granule
#   aggregate clusters to "Fires", assigns landuses and masking status
#   Fire has location, landuse, starttime (UTC of solar midnight) and 
#   "Effective FRP" -- best-fit of prescribed diurnal variation
#   Each fire assumed to live 24hours from its starttime
#   "incomplete fire" -- one whose lifetime is not fully covered with observations



import numpy as np
import subprocess
import io
import datetime as dt
from io import StringIO
import sys
import os
import pandas as pd


import land_use
import ClusterFires
import FireMask
import Firelist2nc

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', -1)

sats="A T 1 N".split() ##N= Suomi National Polar-orbiting Partnership (Suomi NPP), 1=NOAA-20 (designated JPSS-1 prior to launch)
                       ## Aqua Terra

deg2rad = np.pi / 180.
Rearth =6.4e3 #Kilometers
deg2km = Rearth * deg2rad

def xyzFromLonLat(lons,lats):
    #geographic degrees to 3D cartesian in km, assuming sperical earth
    x = Rearth * np.cos(lats * deg2rad) * np.cos(lons * deg2rad)
    y = Rearth * np.cos(lats * deg2rad) * np.sin(lons * deg2rad)
    z = Rearth * np.sin(lats * deg2rad)
    return x,y,z


def FireClustersFromFIRMS(df):

    ## Returns array of strings describing clusters

    clusters = []
    Verbose = False
#    Verbose = True



    dflons = df.longitude.to_numpy()
    dflats = df.latitude.to_numpy()
    dffrps = df.frp.to_numpy()
    dfsats = df.satellite.astype(str).to_numpy()
    dfdaynight = df.daynight.astype(str).to_numpy()
#    print(dfsats)
#    import pdb; pdb.set_trace()
   
    dfsats[dfsats == "Terra"] = 'T'
    dfsats[dfsats == "Aqua"] = 'A'
    dfsats[dfsats == "VJ114IMG_NRT"] = '1'
    dfsats[dfsats == "VNP14IMG_NRT"] = 'N'

#    if dfsats[0] in 'Terra Aqua'.split(): ## Can be Terra or Aqua instead of T and A
#        dfsats = np.array([ l[0] for l in dfsats ]) ## First letter
#    elif dfsats[0][0:3] in 'VJ1 VNP'.split():
#        d = {"VJ1":"1", "VNP":"N"}
#        dfsats = np.array([ d[l[0:3]] for l in dfsats ]) ##translate first three letters
    dfsize = np.sqrt((df.scan*df.track).to_numpy()) ## in km


    addGlint = ("glintAz" in df.keys()) 

    try: ## Integer  Comes from FIRMS file
        timestamps = np.array([int(dt.datetime.strptime( "%s %04d"%(d,t), "%Y-%m-%d %H%M").replace(tzinfo=dt.timezone.utc).timestamp()) for d,t in zip(df.acq_date,df.acq_time)])
    except TypeError:
        if df["acq_time"][1][2] == ':': ## HH:MM Comes from NRT laads
          timestamps = np.array([int(dt.datetime.strptime( "%s %s"%(d,t), "%Y-%m-%d %H:%M").replace(tzinfo=dt.timezone.utc).timestamp()) for d,t in zip(df.acq_date,df.acq_time)])
        else:
            raise



    ## Clasrterize pixels to "fires" for each satellite:
    ## Pixels closer than 0.5 size -- duplicate and average
    ## Pixels closer than 1.5 size -- aggregate
    ## Method:
    ## Loop over pixels of given satellite 
    ## and aggregate within overlapping time windows for each satellite
    ## to ensure that no fire is cut between timetamp
    ## Window -- one or two timestamps

    sortedTindex = np.argsort(timestamps, kind='stable') ## Sort stuff over time

    sortedsteps = timestamps[sortedTindex]
    utimestamps = np.sort(np.unique(timestamps))
    nutimestamps = len(utimestamps)

    WindowClusters = {}

    barktimestep =  utimestamps[0] + 24*3600

    for iut in range(nutimestamps): ## Loop over sorted unique timestamps

        if (utimestamps[iut] > barktimestep):
                print("Processing pixels %s, end %s"%(
                    dt.datetime.utcfromtimestamp( utimestamps[iut] ).strftime("%Y-%m-%dT%H:%M"),
                    dt.datetime.utcfromtimestamp( utimestamps[-1] ).strftime("%Y-%m-%dT%H:%M"),
                    ))
                barktimestep += 24*3600
       
        if iut == 0 or utimestamps[iut] - utimestamps[iut-1] > 120: 
              ## No tail from previous timestamp expected. Aggregate only this timstamp 
 #           idx, = np.where(sortedsteps == utimestamps[iut])
            itstartsorted = np.searchsorted(sortedsteps, utimestamps[iut], side='left')
            itendsorted = np.searchsorted(sortedsteps, utimestamps[iut], side='right')
        else:
            # Aggregate with previous
#            idx, = np.where((sortedsteps == utimestamps[iut-1]) + (sortedsteps == utimestamps[iut]) )
            itstartsorted = np.searchsorted(sortedsteps, utimestamps[iut-1], side='left')
            itendsorted = np.searchsorted(sortedsteps, utimestamps[iut], side='right')

        if iut == nutimestamps-1 or utimestamps[iut+1] - utimestamps[iut] > 120:
            dumpAll = True  ## No need to store: no clusters will extend to the next timestamp
        else:
            dumpAll = False 

        ## Regardless a satellite
        prevWindowClusters = WindowClusters.copy()  ##
        WindowClusters = {} ## Remember clusters found in this window
       
        sliceidx =  sortedTindex[itstartsorted:itendsorted]
        for sat in sats:
            myidx = sliceidx[dfsats[sliceidx] == sat ] ##Subset frame for this satellite
            if len(myidx) < 1:
                continue


            x,y,z = xyzFromLonLat(dflons[myidx], dflats[myidx])
            size = dfsize[myidx] ## in km
            frp = dffrps[myidx]

        #    icluster = clusterize(x,y,z,size)
            npixels = np.int32(len(x)) # or whatever 
                
            icluster = np.zeros((npixels,),dtype=np.int32) - 1  ## Cluster no
            idup = icluster.copy()                           ## Duplicate pixel
            dupdist = np.zeros((npixels,),dtype=np.float32) + np.nan ## Distance from duplicate
            if Verbose:
                print("# Aggregating %d pixels( %d --%d)"%(npixels, myidx[0], myidx[-1]))

            ClusterFires.cluster(x,y,z, size, frp, myidx, icluster, idup, dupdist,   npixels)
            nclusters = np.amax(icluster)+1  ## I could nnot figure out how to output nclusters frpom Fortran

            ## Aggregate found clusters and dump original info
            for icl in range(nclusters):
                clusteridx = myidx[ icluster == icl ]  ## Select this cluster
                n = len(clusteridx)
                if n < 1:
                    continue ## No empty clusters (could appear from merging)
                
                if not dumpAll: ## Need to store clusters we have seen
                    k = " ".join([ "%d"%i for i in sorted(clusteridx) ]) ## Key for a cluster
                    WindowClusters[k] = 1
                if dumpAll or (k in prevWindowClusters):  ##  This window added nothing
                    clusteridxnodup = myidx[ (icluster == icl) * (idup < 0)] 
                    seldup =  (icluster == icl) * (idup >= 0)  ## Boolean array
                    Verbose = (np.sum(seldup) > 100)
                    if Verbose:
                        # Report luster content
                        print("------------------------")
                        print(df.loc[clusteridxnodup])
                        for i, ref, dist in zip(myidx[seldup], idup[seldup], dupdist[seldup]):
                            print("#" + df.loc[i:i].to_string(header=False) + "#Duplicate of %d, dist=%7.3fkm! "%(ref,dist))
                    lons = dflons[clusteridxnodup] 
                    lats = dflats[clusteridxnodup] 
                    frps = dffrps[clusteridxnodup] 

                    # These two should sum up to the observed FRP
                    frptot = np.sum(frps)                 ## Cluster FRP
                    dupfrp = np.sum(dffrps[idup[seldup]]) ## FRP in dups  

                    if (frptot > 1e-3):
                        lomean  = np.sum(lons*frps)/frptot
                        lamean  = np.sum(lats*frps)/frptot

                        pixsize = np.mean(dfsize[clusteridxnodup])  ## ~"typical" pix size
                        daynight = dfdaynight[clusteridxnodup[0]] # Daynight flag

                        clustersize =  np.sum((((lons - lomean)*np.cos(lamean*deg2rad))**2 + ((lats - lamean))**2)*frps) / frptot  ## radius deg2
                        clustersize = 2*np.sqrt(clustersize) * deg2km + pixsize ## diameter to km , Add pixel size

                        clustertstamp =  timestamps[clusteridxnodup[0]] 
                        acq_time =  dt.datetime.utcfromtimestamp(clustertstamp).strftime("%Y-%m-%dT%H:%M:%SZ")
                        LTime =  dt.datetime.utcfromtimestamp(clustertstamp + int(lomean * 240) ).strftime("%Y-%m-%dT%H:%MLT")

                        if addGlint:
                            glintAz=df['glintAz'][clusteridxnodup[0]] ## Take one
                            glintZe=df['glintZe'][clusteridxnodup[0]]
                        
                            # Cluster summary
                            clusterstr =  "%s,%10.5f,%10.5f,%10.3f,%7.3f,%7.3f,%7.3f,%7.3f,%1s,%1s,%s,%d,%d,%10.3f"%(
                                          LTime, lamean, lomean, frptot, clustersize, pixsize,glintZe,glintAz,sat,daynight,acq_time, len(lons), np.sum(seldup), dupfrp
                                          )
                        else:
                            # Cluster summary
                            clusterstr =  "%s,%10.5f,%10.5f,%10.3f,%7.3f,%7.3f,%1s,%1s,%s,%d,%d,%10.3f"%(
                                          LTime, lamean, lomean, frptot, clustersize, pixsize, sat,daynight,acq_time, len(lons), np.sum(seldup), dupfrp
                                          )
                        clusters.append(clusterstr)

                        if Verbose:
                            clusterstr =  "ClusterFRP %10.5f %10.5f %10.1f %7.3f %1s %s # %d fires + %d dups"%(
                                      lamean, lomean, frptot, clustersize, sat, acq_time, len(lons), np.sum(seldup))
                            print(clusterstr)
                    elif Verbose: 
                         print("Zero FRP, skipping")
    if addGlint:
      head = "LTime,lat,lon,frptot,clsize,pixsize,glintZe,glintAz,sat,daynight,acq_timestr,nPix,nDup,dupFRP\n"
    else:
      head = "LTime,lat,lon,frptot,clsize,pixsize,sat,daynight,acq_timestr,nPix,nDup,dupFRP\n"

    # Return a dataframe of values
    return pd.read_csv(StringIO(head + '\n'.join(clusters)))


def ClustersToFires(df, firemetadata, firemask, ifFitRMSE=False):

    # Converts cluster observations to "Fire descriptions"
    # Here fire is something that lives from midnight local to midnight local
    # and has prescribed diurnal cycle (landuse specifric, from metadata)
    # Aggregates in space and time, assigns landuse
    #
    # firemetadata --  land_use object
    # firemask -- FireMask object
    # 
    #


    Verbose = False


    dflats = df.lat.to_numpy()
    dflons = df.lon.to_numpy()
    dffrps = df.frptot.to_numpy()
    dfsats = df.sat.astype(str).to_numpy()
    dfdaynight = df.daynight.astype(str).to_numpy()
    dfsize = df.clsize.to_numpy() ## in km
    timestamps = np.array([int(dt.datetime.strptime( s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc).timestamp()) for s in df.acq_timestr])
    localDays =  np.int32((timestamps + dflons * 240)//86400) # 240 sec/deg

    ldstart =  np.amin(localDays)
    ldend   =  np.amax(localDays)

    x,y,z = xyzFromLonLat(dflons, dflats)

    dumpedclusters =np.zeros((len(x),), dtype = np.int8)

    
    fires = []
    for localDay in range(ldstart,ldend+1):

        idxday, = np.where(localDays == localDay)
        sortedidx = np.argsort(z[idxday], kind='stable')  ### Sorting z
        unsortedidx = np.argsort(sortedidx, kind='stable') ## Undoing sorting

        myx = x[idxday][sortedidx]
        myy = y[idxday][sortedidx]
        myz = z[idxday][sortedidx]
        mys = dfsize[idxday][sortedidx]


        nclusters = np.zeros((1,),dtype=np.int32)

        iclusterzorder = np.zeros((len(myx),),dtype=np.int32)

        ClusterFires.clustersoffire(myx,myy,myz,mys, iclusterzorder, nclusters, len(myx)) 

        nclusters = nclusters[0]  #
        iclusterdforder = iclusterzorder[unsortedidx]

        ## Aggregate found clusters and dump original info
        for icl in range(nclusters):
            clusteridx = idxday[ iclusterdforder == icl ]  ## Select this cluster
            n = len(clusteridx)
            if n < 1:
                   continue ## No empty clusters (could appear from merging)
                
            if Verbose:
                # Report cluster content
                print("------------------------")
                print(df.loc[clusteridx])
            if np.any(dumpedclusters[clusteridx] > 0):
                print("Already dumped clusters")
                import pdb; pdb.set_trace()
            dumpedclusters[clusteridx] = 1

            
            lons = dflons[clusteridx] 
            lats = dflats[clusteridx] 
            frps = dffrps[clusteridx] 
            tstamps = timestamps[clusteridx]
            frptot = np.sum(frps)
            if not frptot > 0:
                print("-------FRP <= 0-----------------")
                print(df.reset_index(drop=True).loc[clusteridx])
                raise ValueError("nonpositive FRP")
                
            lomean  = np.sum(lons*frps)/frptot
            lamean  = np.sum(lats*frps)/frptot
            clustersize =  np.sum((((lons - lomean)*np.cos(lamean*deg2rad))**2 + ((lats - lamean))**2)*frps) / frptot  ## radius deg2
            clustersize = 2*np.sqrt(clustersize) * deg2km ## diameter to km


                
            [luidx] = firemetadata.get_LU_4_fires(np.array([lomean]), np.array([lamean]))
            lutype = firemetadata.LUtypes[luidx]


            hrinday = np.array(firemetadata.diurnal[:,luidx], dtype=np.float32)  ##total  (not per_fire)
            hrinday_per_fire = np.array(firemetadata.diurnal_per_fire[:,luidx], dtype=np.float32)  ##(per_fire)
            mask = firemask.get_mask([lomean],[lamean])

            localoff = int(lomean * 240)
            firestart = np.floor((tstamps[0] + localoff) / (24*3600)) * 3600 * 24  - localoff ## Utc time of local midnight, before first obs

            ## Satellites observed the fire
            obssats = "".join(dfsats[clusteridx])
            daynight = "".join(dfdaynight[clusteridx])

            frpmean = np.array([0,0,0], dtype = np.float32) ##
            if n==1:
              frpmean[0] = min(frptot*2, frptot /          hrinday[np.int32((tstamps[0]-firestart)/3600)%24] )
              frpmean[1] = min(frptot*2, frptot / hrinday_per_fire[np.int32((tstamps[0]-firestart)/3600)%24] )
              frpmean[2] = frptot
            else:
              torder = np.argsort(tstamps)
              localt = np.array(tstamps[torder] - firestart, dtype=np.int32) % (24*3600) #obs in seconds since local daystart
              frptorder = np.array(frps[torder], dtype=np.float32)
              ## cluster2fire handles only one timevar
              if ifFitRMSE: ## Minimize RMSE
                  ClusterFires.cluster2firermse(localt, frptorder, hrinday, frpmean[0:-1], n)
                  ClusterFires.cluster2firermse(localt, frptorder, hrinday_per_fire, frpmean[1:], n)
              else: ## MAS fit, as in SILAM
                  ClusterFires.cluster2firesilam(localt, frptorder, hrinday, frpmean[0:-1], n)
                  ClusterFires.cluster2firesilam(localt, frptorder, hrinday_per_fire, frpmean[1:], n)

            # Cluster summary
            clusterstr =  "%s,%10.5f,%10.5f,%10.3f,%10.3f,%10.3f,%7.3f,%2d,%s,%d,%s,%s"%(
                              dt.datetime.utcfromtimestamp(firestart).strftime("%Y-%m-%dT%H:%M"),
                              lamean, lomean, frpmean[0], frpmean[1], frpmean[2], clustersize, luidx, lutype, mask, obssats, daynight)
            fires.append(clusterstr)
#
    if np.any(dumpedclusters == 0):
        print("Lost clusters")
        raise ValueError()
    head = "FireStartZ,lat,lon,FRPeffTotMW,FRPeffPerFireMW,FRPmeanOBS,size,luidx,lutype,mask,obssats,daynight\n" 

    return pd.read_csv(StringIO(head + '\n'.join(fires)))


def expandClusterDF(clusterdf, expand_days):

    isoformatZ  = '%Y-%m-%dT%H:%M:%SZ'
    isoformatLT = '%Y-%m-%dT%H:%MLT'
    one_day = dt.timedelta(days=1)

    acqtimeZ = dt.datetime.strptime(clusterdf.iloc[-1]["acq_timestr"], isoformatZ)
    prevdayacktimeZ = (acqtimeZ - one_day).strftime(isoformatZ)
    repdf = clusterdf[clusterdf["acq_timestr"] > prevdayacktimeZ ] 
   
    ltimes = [ dt.datetime.strptime(t, isoformatLT)  for t in  repdf['LTime'] ]
    acqtimes = [ dt.datetime.strptime(t, isoformatZ)  for t in  repdf['acq_timestr'] ]

    dflist = [ clusterdf ]
    for d in range(expand_days):
        fakedf = repdf.copy()
        fakedf.set_index(repdf.index + (d+1)*len(ltimes), inplace=True)
        fakedf['LTime'] = [ (t + (d+1)*one_day).strftime(isoformatLT) for t in  ltimes ]
        fakedf['acq_timestr'] = [ (t + (d+1)*one_day).strftime(isoformatZ) for t in  acqtimes ]
        dflist.append(fakedf)

    return pd.concat(dflist)


def getFullDaysFires(dfFires, dfFiesExt ):
    # dfFires "observed fires"
    # dfFiesExt -- fires made from extended obs

    # Return a list of daily "complete" dataframes
    # First: UTC day strting at least 24 hours after the first fire start
    # Last: Last day unaffcted by extension
    # 

    assert (len(dfFires) < len(dfFiresExt)) 
    
    # Iterate through the rows and compare
    for i in len(dfFires):
        if not dfFires.iloc[i].equals(dfFiresExt.iloc[i]):
            break

    # i has  an index of the first tow that differ
    # Rows before it have not been affected by the replication.

    fireformat= '%Y-%m-%dT%H:%M'
    times = [ dt.datetime.strptime(t, isoformatLT)  for t in  repdf['FireStartZ'] ]
    


if __name__ == '__main__':

    import sys	
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", metavar="FIRMS-list.csv",  help="FIRMS hotspot lists (in)", nargs="+")
    parser.add_argument("outf", metavar="SILAM-list",  help="SILAM fire list (out) .nc or .csv")
    parser.add_argument("--metadata",  help="Fire metadata")
    parser.add_argument("--datepref",  help="Prefix to select for fire-start UTC datetime e.g. YYYY-MM-DDT", default=None)
    parser.add_argument("--mask",  help="Fire mask in ")
    parser.add_argument("--extend", type=int, default=0, help="extend dataframe with given number of days")
    parser.add_argument("--rmse-fit", "-r", action="store_true",  help="Use RMSE fit for temporal profile")
    parser.add_argument("--verbose", "-v", action="store_true",  help="Verbose")


    args = parser.parse_args(sys.argv[1:])

#    firemetadatafile = "/home/kouzne/00PASO/fires_v2_0/fire_metadata_ecodata_Akagi_PM_v5_7_cbm5.ini"
#    firemetadatafile = "fire_metadata_ecodata_Akagi_continents_withAVB_v2.ini"
    firemetadata = land_use.land_use(args.metadata)
    if args.verbose:
        print("Got metadata from "+args.metadata)


    firemask = FireMask.FireMask.from_nc(args.mask)
    print("got Firemask from "+ args.mask)

    dflist = []
    for f in args.infiles:
          dflist.append( FireClustersFromFIRMS(pd.read_csv(f)))
          if args.verbose:
              print(f, "Done")
    clusterdf = pd.concat(dflist,  ignore_index=True)
    if args.verbose:
        print ("Concat Done")

    if args.verbose: 
        print ("Clusters to fires")
    if args.extend == 0:
      Firelist = ClustersToFires(clusterdf, firemetadata, firemask)
    else:
      newdf = expandClusterDF(clusterdf, args.extend)  ## Extend for x days in the future
      Firelist =  ClustersToFires(newdf, firemetadata, firemask)

    if not args.datepref is None:
        if args.verbose: 
             print ("Selecting "+args.datepref)
        len_in = len(Firelist)
        Firelist = Firelist.loc[Firelist["FireStartZ"].str.startswith(args.datepref) ]
        len_out = len(Firelist)
        if args.verbose:
            print("Selected %d  out of %d fires "%(len_out, len_in))


    ### Dump output
    if args.outf.endswith(".csv"):
        if args.verbose:
            print("Dumping csv output to "+ args.outf)
        Firelist.sort_values(by=["FireStartZ"]).to_csv(args.outf, index=False)
    elif args.outf.endswith(".nc") or  args.outf.endswith(".nc4"):
        if args.verbose:
            print("Dumping nc output to "+ args.outf)
    
        Firelist2nc.FireList2NC(Firelist, args.outf, description=" ".join(sys.argv))


    

