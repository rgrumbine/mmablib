REAL FUNCTION harcdis(lat1, lon1, lat2, lon2)
  USE constants
  IMPLICIT none
  REAL lat1, lon1, lat2, lon2
  REAL dlat, dlon, mlat
  REAL a, c

  dlon = lon2 - lon1
  dlat = lat2 - lat1
  mlat = (lat1 + lat2)/2.

  a = sin(dlat*rpd/2)**2 + cos(lat1*rpd)*cos(lat2*rpd)*sin(dlon*rpd/2)**2
  c = 2.*asin(min(1.,sqrt(a)))

! approximating ellipsoidal flattening RG WGS84
  harcdis = c * (6378.137 - 21.385*sin(mlat*rpd) )
  RETURN
END function harcdis

