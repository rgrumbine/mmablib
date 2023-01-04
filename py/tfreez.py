from math import *
import sys

salinity = float(sys.argv[1])

def tfreeze(salinity):
  a1 = -0.0575
  a2 =  1.710523E-3
  a3 = -2.154996E-4
  
  if (salinity >= 0):
    tfreez = salinity*(a1+a2*sqrt(salinity)+a3*salinity)
  
  return tfreez


print("salinity, freezing point ",salinity, tfreeze(salinity))
