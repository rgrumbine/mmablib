"""
Working with skiles points and drift forecasts
#Robert Grumbine
"""

#from math import *
import sys
from latpt import *
nmtokm = 1.852

#better to give skpt a latpt (has-a)
class skpt:
  ''' class skpt defines a class with the location and wmoid of a point '''

  def __init__(self, lat = 0., lon = 0., wmoid = 0):
    self.lat = lat
    self.lon = lon
    self.wmoid  = int(wmoid)

  def show(self, fout = sys.stdout):
    ''' skpt.show(fout = xxx) -- print out point information to (optional) file xxx '''
    print("skiles point ",self.wmoid, self.lon, self.lat, file = fout)


class skfcst:

   def __init__(self, wmoid=0, lon=0., lat=0., wdir=0., fdist=0.):
      self.wmoid   = wmoid
      self.lon  = float(lon)
      self.lat  = float(lat)
      self.wdir  = float(wdir)
      self.dist = float(fdist)

   def show(self):
      ''' skfcst.show -- print out an skfcst point '''
      print(self.wmoid, self.lon, self.lat, self.wdir, self.dist)

class skiles2_points :
  ''' class skiles2_points -- work on seaice_forecast.points file of points '''

  def __init__(self, fname='seaice_forecast.points'):
    self.model = []
    #
    fin = open(fname, encoding='utf-8')
    count = 0
    x     = skpt()
    self.model.append(x)
    self.model[count] = skpt()
    count = 1
    for line in fin:
      words = line.split()
      wmoid  = words[0]
      lat = float(words[1])
      lon = float(words[2])
      self.model.append(x)
      self.model[count] = skpt(lat, lon, wmoid)
      #self.model[count].show()
      count += 1
    fin.close()

class skiles2_forecast :
  days    =  int(16)
  header  =  int( 5)
  nskiles =  int(207)

  def __init__(self, fname):
    self.refpts = skiles2_points()
    self.fcst = []
    self.npoints = int(0)

    fin = open(fname, encoding='utf-8')
    j   = int(0)
    day = int(0)
    npts_fcst=int(0)
    blank_count = int(0)

    x = skfcst()

    for line in fin:

      if ( j < skiles2_forecast.header or j == (207 + 5) ):
        npts_fcst += 1
      else:
        if ( (j < 207+5)):
          words = line.split()
          wmoid = int(words[0])
          lon = self.refpts.model[wmoid].lon
          lat = self.refpts.model[wmoid].lat
          wdir = float(words[1])
          fdist = float(words[2])*nmtokm
        else:
          if ( ( j > 212 )):
            if (len(line) < 5):
              blank_count += 1
            else:
              words = line.split()
              wmoid = int(words[0])
              lon = float(words[1])
              lat = float(words[2])
              wdir = float(words[3])
              fdist = float(words[4])*nmtokm

        if (blank_count == 0):
          self.fcst.append(x)
          self.fcst[j+int(day*self.npoints)-npts_fcst] = skfcst(wmoid, lon, lat, wdir, fdist)
        elif (blank_count == 2):
          j = -1
          blank_count = 0
          npts_fcst = 0
          day += 1

      j += 1
      if (day != 0 and j == 0) :
        self.npoints = len(self.fcst) / day
# loop back to process next line.
#  n.b. if j = 0, we've reached the end of a day and now know how many
#    points are in the forecast
# all 16 days are in the one long vector

  def nearest(self, lat, lon, tolerance = 50) :
     ''' nearest(lat, lon, tolerance) -- return index to skiles point 
             closest to the given lat,lon, -1 if not within tolerance km ''' 
     fdist = 1.e6
     j = int(-1)
     ref = latpt(lat,lon)
     for i in range (int(0), int(self.npoints) ) :
       x = latpt(self.fcst[i].lat, self.fcst[i].lon)
       tmp = ref.distance(x)
       if (tmp < fdist):
         fdist = tmp
         j = i
     #print("dist ",j,dist,lat,lon)
     if (fdist < tolerance) :
       return j

     return -1
