"""
Compute freezing point as a function of salinity
"""

from math import sqrt
import sys

salinity = float(sys.argv[1])

def tfreeze(salt):
  """
  Compute freezing point as a function of salinity
  """
  a1 = -0.0575
  a2 =  1.710523E-3
  a3 = -2.154996E-4

  if (salt >= 0):
    tfreez = salt*(a1+a2*sqrt(salt)+a3*salt)
  else:
    print("tfreeze negative salinity ",salt)
    tfreez = 0.0

  return tfreez


print("salinity, freezing point ",salinity, tfreeze(salinity))
