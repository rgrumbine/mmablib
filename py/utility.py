"""
#Robert Grumbine
"""

import datetime
from math import *


from const import *

#--------------- Utility Functions --------------------------------
# Convert an 8 digit date to a python date object
def parse_8digits(tag):
  """ Convert an 8 digit int to a datetime.date object """
  tmp = int(tag)
  (yy,mm,dd) = (int(int(tmp)/10000),int((int(tmp)%10000)/100),int(tmp)%100)
  tag_out = datetime.date(int(yy), int(mm), int(dd))
  return tag_out
#------------------------------------------------

#----------------------- for mapping --------------------------
# approximating ellipsoidal flattening in WGS84
def rearth(lat):
  """ return the radius of the earth at the given latitude """
  return 6378.137 - 21.385*sin(lat*const.rpdg)

#haversine arcdis
#  http://www.movable-type.co.uk/scripts/gis-faq-5.1.html
#assumes lat lon in degrees, distance in km
def harcdis(pt1, pt2):
  """ haversine formula for distance between latpts """
  dlon = pt2.lon - pt1.lon
  dlat = pt2.lat - pt1.lat
  mlat = (pt1.lat + pt2.lat)/2.

  a = sin(dlat*const.rpdg/2)**2 + cos(pt1.lat*const.rpdg)*cos(pt2.lat*const.rpdg)*\
                                  sin(dlon*const.rpdg/2)**2
  c = 2.*asin(min(1.,sqrt(a)))

  return c*rearth(mlat)

#--------------------------------------------------------------

#line="46 45 16.4 - 48 46 32.1"
def dms_dpddd(line):
  """convert degree minute second lat-lon pair to decimal degrees"""
  words=line.split()
  lat = float(words[0]) + float(words[1])/60. + float(words[2])/3600.
  lon = float(words[4]) + float(words[5])/60. + float(words[6])/3600.
  lon = -lon
  return (lat, lon)
