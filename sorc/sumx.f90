      FUNCTION sumx(x,n)
!     Function to compute the sum of a vector.  Use double precision
!       summation, but return a single precision answer.  
!     Robert Grumbine 12/09/1993.

      IMPLICIT none

      INTEGER n
      REAL sumx, x(n)

      DOUBLE PRECISION sum
      INTEGER i
      
      sum = 0.
      DO i = 1, n
        sum = sum + DBLE(x(i))
      ENDDO
      sumx = SNGL(sum)
      RETURN
      END
