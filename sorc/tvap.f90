      REAL FUNCTION VAPOR(T,K1)
!=======================================================================
!  PROGRAMMED BY:
!     C.KOCH                  UNI, BONN                             1986
!  MODIFIED BY:
!     A.STOESSEL              MPI, HAMBURG                          1989
!     Robert Grumbine         NOAA/NWS                              1992
!  PURPOSE:
!     -CALCULATION OF SATURATION VAPOR PRESSURE FOR AIR TEMPERATURE
!       (K1=1), OVER ICE (K1=2) AND OVER WATER (K1=3)
!  INTERFACE:
!     -T:   TEMPERATURE OF ATMOSPHERE, ICE OR OCEAN
!     -K1:  INDEX FOR CHOICE OF QUANTITY TO BE CALCULATED
!=======================================================================
      IMPLICIT none
      REAL T
      INTEGER K1

!      GOTO(1,2,3),K1
      IF (K1 .EQ. 1) THEN
        VAPOR=611.21*EXP((18.729-(MIN(T,300.)-273.15)/227.3)*   &
              (MIN(T,300.)-273.15)/(MAX(T,200.)-273.15+257.87))
      ELSE IF (K1 .EQ. 2) THEN
       VAPOR=611.15*EXP((23.036-(MIN(T,273.15)-273.15)/333.7)*   &
              (MIN(T,273.15)-273.15)/(MAX(T,200.)-273.15+279.82))
      ELSE
       VAPOR=0.9815*611.21*EXP((18.729-(MIN(T,300.)-273.15)/227.3)*   &
              (MIN(T,300.)-273.15)/(MAX(T,260.)-273.15+257.87))
      ENDIF

      RETURN
      END
