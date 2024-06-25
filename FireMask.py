#!/usr/bin/env python3
# encoding: utf-8

import netCDF4 as nc4
import numpy as np
import png
import os

class FireMask:
         @classmethod
         def from_nc(cls,firemaskfile, varname=None):

             if varname is None:
                varname = 'mask'

             with nc4.Dataset(firemaskfile) as nc:
                 #get rid of a dummy dimension
                 fm =  cls(nc.variables['lat'][:].data,
                                nc.variables['lon'][:].data,
                                nc.variables[varname][0,:,:])
                 try: 
                     fm.long_name = nc.variables[varname].long_name 
                 except AttributeError:
                     pass

                 try: 
                     fm.comment = nc.comment
                 except AttributeError:
                     pass

                 return fm

         def __init__(self, lats, lons, mask):

             ny,nx = mask.shape
             self.shape = (ny,nx,)
             if nx == len(lons) and  ny == len(lats):
                  self.lons = lons
                  self.lats = lats 
                  dx = (lons[-1] - lons[0])/(nx-1)
                  dy = (lats[-1] - lats[0])/(ny-1)
                  self.lonbnds = np.concatenate( (lons[0:1]-0.5*dx,  0.5*(lons[1:] + lons[0:-1]), lons[nx-1:nx] + 0.5*dx,))
                  self.latbnds = np.concatenate( (lats[0:1]-0.5*dy,  0.5*(lats[1:] + lats[0:-1]), lats[ny-1:ny] + 0.5*dy,))
             elif  nx + 1 == len(lons) and  ny + 1  == len(lats):
                 self.lonbnds = lons
                 self.latbnds = lats
                 self.lons = 0.5*(self.lonbnds[1:] + self.lonbnds[0:-1])
                 self.lats = 0.5*(self.latbnds[1:] + self.latbnds[0:-1])
             else:
                 raise IndexError("Wrong shapes")
             try: 
                 self.mask = np.int8(mask.filled(fill_value=0))
             except AttributeError:
                 self.mask = np.int8(mask)


         def get_mask(self, lats, lons):
             # 
             # Get mask status for array of coordinates
             idxlon = np.searchsorted(self.lonbnds[1:-1], lons)
             idxlat = np.searchsorted(self.latbnds[1:-1], lats)
             return self.mask[idxlat,idxlon]


         def to_nc(self, filename, comment=None, long_name = None):

           with  nc4.Dataset(filename, "w", format="NETCDF4_CLASSIC") as dst:
             # copy attributes
             if not comment is None:
                  dst.setncattr("comment", comment)

             dst.createDimension("time", None)
             dst.createDimension("lon", len(self.lons))
             dst.createDimension("lat", len(self.lats))

             t = dst.createVariable('time', np.int32, ("time",))
             t.standard_name = "time" 
             t.long_name = "dummy time dimension"
             t.units = "hours since 1-1-1 00:00:00" ;
             t.calendar = "standard" ;
             t.axis = "T" ;
             t[:] = [0]

             lo =  dst.createVariable('lon', np.float32, ("lon",), complevel = 5, zlib=True)
             lo.standard_name = "longitude" ;
             lo.long_name = "longitude" ;
             lo.units = "degrees_east" ;
             lo.axis = "X" ;
             lo[:] = self.lons[:]

             la =  dst.createVariable('lat', np.float32, ("lat",), complevel = 5, zlib=True)
             la.standard_name = "latitude" ;
             la.long_name = "latitude" ;
             la.units = "degrees_north" ;
             la.axis = "Y" ;
             la[:] = self.lats[:]

             mask = dst.createVariable('mask', np.int8, ("time","lat","lon",), complevel = 5, zlib=True)
             if not long_name is None: 
                 mask.long_name = long_name
             mask[0,:,:] = self.mask[:]

         def to_KML(self,basename, pngrefname = None):

            kmlname = "%s.kml"%(basename)
            pngname = "%s.png"%(basename)
            if pngrefname is None:
                pngrefname = pngname
            
            try:
                title = self.long_name
            except AttributeError:
                 title = "Fire Mask"
            

            ##Debug
            if False:
                j = np.searchsorted(self.latbnds[1:],np.arange(-80,81,1))
                self.mask[j,:] = 1
                self.mask[j-1,:] = 1
                #print (np.int32(self.lats[j]*1000))

                i = np.searchsorted(self.lonbnds[1:], np.arange(-179,180,1))
                #print (np.int32(self.lons[i-1]*1000))
                self.mask[:,i-1] = 1
                self.mask[:,i] = 1
                print (self.lonbnds[0], self.lonbnds[-1])

            try:
                description = self.comment
            except AttributeError:
                description = "Fire Mask dummy description"

            kmltext = kml_template % dict(
                            name = title,
                            description = description, 
                            N = self.latbnds[-1],
                            S = self.latbnds[0],
                            E = self.lonbnds[-1],
                            W = self.lonbnds[0],
                            pngname = os.path.basename(pngrefname),
                        )

            with open(kmlname,'wt') as kml:
                kml.write(kmltext)



            ny, nx = self.mask.shape

            palette=[(0x00,0x00,0x00,0x00), (0xff,0x00,0x00,0x50)]
            nup = 1 ## Tile images
            w = png.Writer(nx*nup, ny*nup,  palette=palette, bitdepth=1)
            with open(pngname,'wb') as f:
                #flip upside-down (.png order) and replicate each pixel
                if np.amax(self.mask) == 1:
                    w.write(f, np.flip(self.mask, axis=0).repeat(nup*nup).reshape(ny,nx,nup,nup).transpose(0,2,1,3).reshape(ny*nup,nx*nup))
                else:
                    w.write(f, 1*(self.mask > 6))

         def to_KML2(self,basename):

            kmlname = "%s.kml"%(basename)

            try:
                title = self.long_name
            except AttributeError:
                 title = "Fire Mask"

            try:
                description = self.comment
            except AttributeError:
                description = "Fire Mask dummy description"

            with open(kmlname,'wt') as kml:
                head="""<?xml version="1.0" encoding="UTF-8"?>
                    <kml xmlns="http://earth.google.com/kml/2.0"> <Document>
                    <Placemark>
                        <name>%s</name>
                    <Style>
                      <LineStyle>
                        <color>ff000000</color>
                      </LineStyle>
                      <PolyStyle>
                        <color>7f9f7fe0</color>
                      </PolyStyle>
                    </Style>
                    <MultiGeometry>
                 """%(title)
                kml.write(head)
                z = np.where(self.mask > 0)
                for iy,ix in zip(z[0], z[1]):
                    kml.write(makePixel(ix,iy,self.lonbnds,self.latbnds))
                tail="\n</MultiGeometry></Placemark> </Document> </kml>\n"
                kml.write(tail)

def makePixel(ix,iy,xbnds,ybnds):
    return """
<Point> <coordinates>%(x0)f,%(y0)f</coordinates> </Point>
<Polygon> <outerBoundaryIs>  <LinearRing> <coordinates> %(x1)f,%(y1)f %(x2)f,%(y2)f %(x3)f,%(y3)f %(x4)f,%(y4)f %(x1)f,%(y1)f </coordinates> </LinearRing> </outerBoundaryIs> </Polygon>"""%dict(
         x0=0.5*(xbnds[ix]+xbnds[ix+1]), 
         y0=0.5*(ybnds[iy]+ybnds[iy+1]), 
         x1=xbnds[ix], y1=ybnds[iy],
         x2=xbnds[ix+1], y2=ybnds[iy],
         x3=xbnds[ix+1], y3=ybnds[iy+1],
         x4=xbnds[ix], y4=ybnds[iy+1],
          )





kml_template="""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Folder>
    <name>Fire Mask</name>
    <description>KML generated by SILAM team</description>
    <GroundOverlay>
      <name>%(name)s</name>
      <description>  %(description)s
          </description>
      <Icon>
        <href>%(pngname)s</href>
      </Icon>
      <LatLonBox>
        <north>%(N)f</north>
        <south>%(S)f</south>
        <east>%(E)f</east>
        <west>%(W)f</west>
      </LatLonBox>
    </GroundOverlay>
  </Folder>
</kml>
"""
             
if __name__ == '__main__':


    ncfile = "Mask20.nc"
    Newmask = FireMask.from_nc(ncfile)
    Newmask.to_KML2(ncfile)
    #Newmask.to_KML("/home/kouzne/mnt/silam.fmi.fi/roux/kml/mask", "https://silam.fmi.fi/roux/kml/mask.png" )
