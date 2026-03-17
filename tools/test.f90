PROGRAM test
  IMPLICIT none
  INTEGER i

  !functions:
  REAL tfreez, vapor, iagree, sumx, sumx2, sumxy
  REAL new_wdir, wdir, arcdis
  REAL(kind=selected_real_kind(14,300)) :: darcdis
  COMPLEX ciagree

  REAL u, v, d
  REAL lat1, lon1, lat2, lon2
  REAL tmp

  REAL(kind=selected_real_kind(14,300)) :: dlat1, dlon1, dlat2, dlon2
  REAL(kind=selected_real_kind(14,300)) :: dtmp

  INTEGER npts
  PARAMETER(npts = 63)
  REAL dir1(npts), dir2(npts), dist1(npts), dist2(npts)
  REAL x1(npts), y1(npts), x2(npts), y2(npts)
  REAL ia, r2, vcor
  REAL xbar, ybar, sig2x, sig2y
  COMPLEX ctmp, cdist1(npts), cdist2(npts)

  u = 1.
  v = 1.
  tmp = new_wdir(u,v,d)  
  PRINT *,'new_wdir',tmp,d
  tmp = wdir(u,v,d)
  PRINT *,'wdir',tmp,d

  lat1 = 45.0
  lon1 = 0.0
  lat2 = 60.0
  lon2 = 90.0
  tmp = arcdis(lon1,lat1, lon2, lat2)
  PRINT *,'arcdis ',tmp
  dlat1 = 45.0
  dlon1 = 0.0
  dlat2 = 60.0
  dlon2 = 90.0
  dtmp = darcdis(dlon1,dlat1, dlon2, dlat2)
  PRINT *,'darcdis ',tmp

  DO i = 1, npts
    dist1(i) = 23.+i
    dist2(i) = 32.+i
    cdist1(i) = CMPLX(dist1(i))
    cdist2(i) = CMPLX(dist2(i))
    dir1(i)  = 45.+i
    dir2(i)  = 90.+i
  ENDDO
  CALL ssanaly(dist1, dir1, dist2, dir2, npts, ia, r2, vcor)
  PRINT *,'ssanaly ia,r2,vcor ',ia, r2, vcor
  tmp = iagree(dist1, dist2, npts)
  PRINT *,'ia alt ',tmp
  ctmp =  ciagree(cdist1, cdist2, npts)
  PRINT *,'ciagree ',ctmp
  CALL correl(dist1, dist2, npts, r2, xbar, ybar, sig2x, sig2y)
  PRINT *,'correl ',r2,xbar,ybar,sig2x,sig2y
  CALL fit(dist1, dist2, npts, xbar, ybar, r2)
  PRINT *,'fit ',xbar, ybar, r2
  CALL vectorize(dist1, dir1, x1, y1, npts)
  CALL vectorize(dist2, dir2, x2, y2, npts)
  CALL VCC(x1, y1, x2, y2, npts, vcor)
  PRINT *,'vcc alt ',vcor

  tmp = tfreez(35.0)
  PRINT *,'tfreez of 35.0 ',tmp

  tmp = vapor(273.15,3)
  PRINT *,'saturation vapor pressure over water at 0 C ',tmp 

  CALL vectorize(dist1, dir1, dist2, dir2, npts)
  PRINT *,'vectorize x,y ',dist2, dir2

  PRINT *,'sumx of dist1 ',sumx(dist1,npts)
  PRINT *,'sumx2 of dist1 ',sumx2(dist1,npts)
  PRINT *,'sumxy of dist1 dist2',sumxy(dist1,dist2,npts)

  
END

! 4 -r--r--r-- 1 Robert.Grumbine climate 1139 Jun  9  2018 new_wdir.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate 1129 Jun  9  2018 wdir.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate 2704 Apr 25  2019 arcdis.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate 2022 Feb 23 15:03 darcdis.f90
! 4 -rw-r--r-- 1 Robert.Grumbine climate 1108 Sep 15  2022 ssanaly2.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate  605 Sep 15  2022 tfreez.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate 1265 Sep 15  2022 tvap.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate  510 Sep 15  2022 vectorize.f90
! 8 -r--r--r-- 1 Robert.Grumbine climate 6841 Feb 23 14:58 VCC.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate  809 Feb 23 14:58 ciagree.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate  941 Feb 23 14:58 correl.f90
! 4 -rw-r--r-- 1 Robert.Grumbine climate  803 Feb 23 14:58 fit.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate  788 Feb 23 14:58 iagree.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate  420 Feb 23 14:59 sumx.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate  453 Feb 23 14:59 sumx2.f90
! 4 -r--r--r-- 1 Robert.Grumbine climate  445 Feb 23 14:59 sumxy.f90


! No test
!! 4 -rw-r--r-- 1 Robert.Grumbine climate  796 Jun 28  2021 waldrop_mse.f90
!!12 -rw-r--r-- 1 Robert.Grumbine climate 9715 Sep 15  2022 GRIBIT.f90
!! 8 -r--r--r-- 1 Robert.Grumbine climate 4115 Sep 15  2022 mapxy.f90
!! 4 -r--r--r-- 1 Robert.Grumbine climate  820 Sep 15  2022 wmoout.f90
!! 4 -r--r--r-- 1 Robert.Grumbine climate 1521 Sep 15  2022 bsslz1.f90
!! 4 -r--r--r-- 1 Robert.Grumbine climate 1433 Feb 23 14:59 gaulat.f90
!! 8 -rw-r--r-- 1 Robert.Grumbine climate 5102 Feb 23 15:11 w3ft01.f90
