      FUNCTION sumxy(x, y, n)
!     Function to compute the sum of the product of two vectors.
!       Use double precision summation, but return a single precision answer.  
!     Robert Grumbine 12/09/1993.

      IMPLICIT none

      INTEGER n
      REAL x(n), y(n), sumxy

      INTEGER i
      DOUBLE PRECISION sum

      sum = 0.0
      DO i = 1, n
        sum = sum + x(i)*y(i)
      ENDDO

      sumxy = REAL(sum)
      RETURN
      END
