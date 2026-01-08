"""
Class latpt -- working with latitude longitude points
#Robert Grumbine
#1 June 2018
"""

from math import sin, cos, asin, acos, atan2

from const import *

class latpt:
  """
  Class latpt -- working with latitude longitude points
  #Robert Grumbine
  #1 June 2018
  """

  def __init__(self,lat = 0., lon = 0.):
    self.lat = float(lat)
    self.lon = float(lon)

  def show(self):
    """
    latpt.show prints out the longitude and latitude
    """
    print (str(self.lon),str(self.lat))

  def distance(self, x ):
    """ 
    latpt.distance computes the distance in km from self to argument x
    """
### These are in the const.py
#    earth_radius = 6371.2 #km
#    rpdg = 3.1416 / 180.
##
    ab = (90.-self.lat)*const.rpdg
    ac = (90.-x.lat)   *const.rpdg
    bc = abs(self.lon - x.lon)*const.rpdg

    arg = cos(ab)*cos(ac)+sin(ab)*sin(ac)*cos(bc)
    arg = min(1, arg)
    arg = max(-1, arg)

    return const.earth_radius * acos(arg)

#Direction in degrees, distance in km (const.earth_radius)
  def bearing_from(self, distance, direction) :
     """
     Compute the bearing given distance and direction in degrees and km, 
       respectively
     """
     final      = latpt()
     R          = const.earth_radius
     direction *= const.rpdg  #beware conventions on direction
     lat1 = self.lat*const.rpdg
     lon1 = self.lon*const.rpdg
     lat2 = asin( sin(lat1)*cos(distance/R) + cos(lat1)*sin(distance/R)*cos(direction))
     lon2 = lon1 + atan2(sin(direction)*sin(distance/R)*cos(lat1), \
                         cos(distance/R)-sin(lat1)*sin(lat2))
     final.lat = lat2 / const.rpdg
     final.lon = lon2 / const.rpdg
     return final
