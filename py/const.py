"""
class const to hold a number of constants useful to earth mapping et al.
#Robert Grumbine
#1 June 2018
"""
import math

class const:
    """
    class const to hold a number of constants useful to earth mapping et al.
    """

    rpdg           = math.pi / 180.
    dgpr           = 180. / math.pi

    earth_radius   = 6371.2   # km -- sphere
    earth_spheroid = 6378.137 # km -- oblate spheroid (WGS84)
    eccen2         = 0.00669438

    #Unit conversion:
    kmtonm = 1./1.852  #multiply km by this

    #degree_area    = earth_radius*earth_radius*4*math.pi*math.pi / 360./180.
    degree_area    = (2.*math.pi*earth_radius / 360.)**2
#Computed parameters
    eccen    = math.sqrt(eccen2)
    ps_chi1  =  eccen2/2. + 5.*eccen2*eccen2/24. + eccen2*eccen2*eccen2/12.
    ps_chi2  =  7.*eccen2*eccen2/48. + 29.*eccen2*eccen2*eccen2/240.
    ps_chi3  =  7.*eccen2*eccen2*eccen2/120.
    ps_ll    = math.pow( math.pow(1.+eccen,1.+eccen) * math.pow(1.-eccen,1.-eccen) , 1./2.)
