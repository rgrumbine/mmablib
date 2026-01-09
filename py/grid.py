"""
module for grid classes. Includes only the mapping, not the data
#Robert Grumbine
#1 June 2018
"""

#from const import *
import const
from ijpt import *
from latpt import *

def ok(x, loc):
  """ return true if ijpt 'loc' is inside the range of the grid x """
  return (loc.i >= 0 and loc.i < x.shape[0] and
          loc.j >= 0 and loc.j < x.shape[1] )

#-----------------------------------------------------------
# RG: make ABC
#Includes only the mapping, not the data
#############################################################
class psgrid:
  """ class psgrid defines polar stereographic grids """

  def locate(self, i, j, z):
    ''' psgrid.locate(i,j,z) -- copy i,j in to z.lat, z.lon '''
    z.lat = i
    z.lon = j

  def inv_locate(self, z):
    ''' psgrid.inv_locate(i,j,z) -- copy z.lat, z.lon to i,j ''' 
    i = z.lon
    j = z.lat
    return (i,j)

class llgrid:
  """ class llgrid defines lat-lon grids """

  def __init__(self, dlat, dlon, firstlat, firstlon, nx, ny):
    #debug print("hello from llgrid.__init__", flush=True)
    self.dlat         = dlat
    self.dlon         = dlon
    self.firstlat     = firstlat
    self.firstlon     = firstlon
    self.nx           = nx
    self.ny           = ny
    #debug print("degree area = ",const.degree_area, flush=True)
    self.darea_base   = abs(self.dlat*self.dlon)*const.degree_area
    self.dlat_rad     = self.dlat * const.rpdg
    self.firstlat_rad = self.firstlat * const.rpdg

  def locate(self, i, j, z):
    ''' llgrid locate(i,j,z) -- return latpt z of point i,j '''
    z.lat = self.firstlat + j*self.dlat
    z.lon = self.firstlon + i*self.dlon
    while (z.lon > 360.):
      z.lon -= 360.

  def inv_locate(self, lat, lon):
    ''' llgrid inv_locate(lat, lon) -- return (i,j) of location lat,lon '''
    j = (lat - self.firstlat)/self.dlat
    if (lon < 0):
      i = self.nx - (self.firstlon - lon)/self.dlon
    else:
      i = (lon - self.firstlon)/self.dlon

    return(i,j)

  def cellarea(self, j, i):
    ''' llgrid cellarea(j,i) -- return cell area in km^2 of grid point j,i '''
    #Original:
    #tlat = self.firstlat + j*self.dlat
    #dy = (self.dlat) * 111.1
    #dx = (self.dlon) * 111.1 * cos(tlat*const.rpdg )
    #return abs(dx*dy)
    #Promote/pre-compute grid constants (darea_base) and
    #    pre-translate degrees to radians in the base class (*lat_rad)
    return self.darea_base * cos(self.firstlat_rad + j*self.dlat_rad)

#############################################################
class global_nthdeg(llgrid):
  """ 
  llgrid.global_nthdeg -- 1/n degree lat-lon grid
  x = global_nthdeg(8) for a 1/8th degree grid, forex.
  """
  def __init__(self, n = 1.0):
    #self.dlat = -1./float(n)
    #self.dlon =  1./float(n)
    #self.firstlat = 90 + self.dlat/2.
    #self.firstlon = self.dlon / 2.
    #self.nx = 360*n
    #self.ny = 180*n
    dlat = -1./float(n)
    dlon =  1./float(n)
    firstlat = 90 + dlat/2.
    firstlon = dlon / 2.
    nx = 360*n
    ny = 180*n
    super().__init__(dlat, dlon, firstlat, firstlon, nx, ny)



class global_5min(global_nthdeg):
  """ llgrid.global_5min -- 5 arcminute lat-lon grid """
  def __init__(self):
    super().__init__(12)


class global_halfdeg(global_nthdeg):
  """ llgrid.global_halfdeg -- half degree lat-lon grid """
  def __init__(self):
    super().__init__(2)

class global_qdeg(global_nthdeg):
    ''' global_qdeg -- global quarter degree grid '''
    def __init__(self):
        super().__init__(4)
