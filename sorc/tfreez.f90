      REAL FUNCTION tfreez(salinity)
!     Constants taken from Gill, 1982.
!     Author: Robert Grumbine
!     LAST MODIFIED: 1 February 2002.
!          MODIFIED: 21 September 1994.

      IMPLICIT none

      REAL salinity
      REAL a1, a2, a3
      PARAMETER (a1 = -0.0575)
      PARAMETER (a2 =  1.710523E-3)
      PARAMETER (a3 = -2.154996E-4)

      IF (salinity .LT. 0.) THEN
        tfreez = 0.
!D        PRINT *,'tfreez was passed a negative salinity',salinity
        salinity = 1.e-5
      ELSE
        tfreez = salinity*(a1+a2*SQRT(salinity)+a3*salinity)
      ENDIF

      RETURN
      END
