      SUBROUTINE ssanaly(dist1, dir1, dist2, dir2, npts, ia, r2, vcor)
!     Conduct simple analyses of displacement fields
!     1- correlation between 1 and 2 displacements
!     2  index of agreement
!     3  vector correlation.
!     Robert Grumbine 21 April 1994.
!     Revised for generality in working with buoys as well as forecast models.
!     Robert Grumbine 10 April 1995.
!     Variant ssanaly derived from sanaly for use with a single vector.
!     Robert Grumbine 10 April 1995.

      IMPLICIT none

      INTEGER, intent(in) :: npts
      REAL, intent(in) :: dir1(npts), dist1(npts)
      REAL, intent(in) :: dir2(npts), dist2(npts)
      REAL, intent(out) :: ia, r2, vcor

      REAL, allocatable :: x1(:), x2(:), y1(:), y2(:)
      REAL iagree
      REAL xbar, ybar, sig2x, sig2y
    
      ia = iagree(dist1, dist2, npts)
   
      CALL correl(dist1, dist2, npts, r2, & 
                            xbar, ybar, sig2x, sig2y)

      ALLOCATE(x1(npts), x2(npts), y1(npts), y2(npts))
      CALL vectorize(dist1, dir1, x1, y1, npts)
      CALL vectorize(dist2, dir2, x2, y2, npts)
      CALL vcc(x1, y1, x2, y2, npts, vcor)
      DEALLOCATE(x1, x2, y1, y2)
   
      RETURN
      END
