! Adapted from code written by John Burkardt after Fred Fritsch.

SUBROUTINE spline_pchip_set( n, x, f, d )

! Sets derivatives for a piecewise cubic Hermite interpolating polynomial (PCHIP).
! The resulting polynomial has desirable shape and monotonicity properties.

  IMPLICIT NONE

  INTEGER :: n    ! number of data points
  REAL    :: x(n) ! strictly increasing independent variable values
  REAL    :: f(n) ! dependent variable values to be interpolated
  REAL    :: d(n) ! derivative values at the data points

  INTEGER :: i
  INTEGER :: ierr
  INTEGER :: nless1
  REAL    :: del1
  REAL    :: del2
  REAL    :: dmax
  REAL    :: dmin
  REAL    :: drat1
  REAL    :: drat2
  REAL    :: dsave
  REAL    :: h1
  REAL    :: h2
  REAL    :: hsum
  REAL    :: hsumt3
  REAL    :: temp
  REAL    :: w1
  REAL    :: w2
  REAL    :: pchst

  ierr   = 0
  nless1 = n - 1
  h1     = x(2) - x(1)
  del1   = (f(2) - f(1)) / h1
  dsave  = del1

  h2 = x(3) - x(2)
  del2 = (f(3) - f(2)) / h2

! Set d(1) via non-centered three point formula, adjusted to be shape preserving.

  hsum = h1 + h2
  w1   = (h1 + hsum) / hsum
  w2   = -h1 / hsum
  d(1) = w1 * del1 + w2 * del2

  IF( pchst( d(1), del1 ) <= 0. )THEN

    d(1) = 0.

  ELSE IF( pchst( del1, del2 ) < 0. )THEN ! monotonicity switches direction

    dmax = 3.* del1

    IF( ABS( dmax ) < ABS( d(1) ) ) d(1) = dmax

  END IF
!
!  Loop through interior points.
!
  DO i = 2, nless1

    IF( i > 2 )THEN
      h1   = h2
      h2   = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = (f(i+1) - f(i)) / h2
    END IF

!  Set d(i) = 0 unless data are strictly monotonic

    d(i) = 0.
    temp = pchst( del1, del2 )

    IF( temp < 0. )THEN
      ierr  = ierr + 1
      dsave = del2
    ELSE IF( temp == 0. )THEN ! count number of changes in direction of monotonicity
      IF( del2 /= 0. )THEN
        IF( pchst( dsave, del2 ) < 0. )THEN
          ierr = ierr + 1
        END IF
        dsave = del2
      END IF
    ELSE ! use Brodlie modification of Butland formula
      hsumt3 = 3.0D+00 * hsum
      w1     = (hsum + h1) / hsumt3
      w2     = (hsum + h2) / hsumt3
      dmax   = MAX ( ABS( del1 ), ABS( del2 ) )
      dmin   = MIN ( ABS( del1 ), ABS( del2 ) )
      drat1  = del1 / dmax
      drat2  = del2 / dmax
      d(i)   = dmin / (w1 * drat1 + w2 * drat2)
    END IF

  END DO

!  Set d(n) via non-centered three point formula, adjusted to be shape preserving

  w1   = -h2 / hsum
  w2   = (h2 + hsum) / hsum
  d(n) = w1 * del1 + w2 * del2

  IF ( pchst( d(n), del2 ) <= 0. ) THEN
    d(n) = 0.
  ELSE IF ( pchst( del1, del2 ) < 0. ) THEN
    dmax = 3.* del2
    IF( ABS( dmax ) < ABS( d(n) ) )THEN ! monotonicity switches direction
      d(n) = dmax
    END IF
  END IF

END SUBROUTINE spline_pchip_set


SUBROUTINE spline_pchip_val( n, x, f, d, xe, fe, de )

! This routine evaluates the cubic Hermite polynomial and its derivative at xe
! If xe > x(n), fe is extrapolated

  IMPLICIT NONE

  INTEGER :: n    ! number of data points
  REAL    :: x(n) ! strictly increasing independent variable values
  REAL    :: f(n) ! function values
  REAL    :: d(n) ! derivative values
  REAL    :: xe   ! point at which the function is to be evaluated
  REAL    :: fe   ! function value at xe
  REAL    :: de   ! derivative at xe

  INTEGER :: klo
  INTEGER :: khi
  INTEGER :: k
  REAL    :: c2
  REAL    :: c3
  REAL    :: del1
  REAL    :: del2
  REAL    :: delta
  REAL    :: h
  REAL    :: xdiff

  IF( xe <= x(1) )THEN
    fe = f(1)
    de = d(1)
  ELSE IF( xe >= x(n) )THEN ! extrapolate
    fe = f(n) + (xe - x(n)) * d(n)
    de = d(n)
  ELSE

!   Find the interval in x that brackets xe by means of bisection

    klo = 1
    khi = n
    DO WHILE( khi - klo > 1 )
      k = (khi + klo) / 2
      IF( x(k) > xe )THEN
        khi = k
      ELSE
        klo = k
      ENDIF
    END DO

!   Evaluate cubic at xe

    h     = x(khi) - x(klo)
    delta = (f(khi) - f(klo)) / h
    del1  = (d(klo) - delta) / h
    del2  = (d(khi) - delta) / h
    c2    = -(del1 + del1 + del2)
    c3    = (del1 + del2) / h
    xdiff = xe - x(klo)

    fe = f(klo) + xdiff * (d(klo) + xdiff * (c2 + xdiff * c3))
    de = d(klo) + xdiff * (2.* c2 + xdiff * 3.* c3)

  END IF

END SUBROUTINE spline_pchip_val


REAL FUNCTION pchst( arg1, arg2 )

! Computes the sign of ARG1 * ARG2.
! Returns +/-1 (or 0 if either/both arguments is/are 0).

  REAL :: arg1
  REAL :: arg2

  IF( arg1 == 0. .OR. arg2 == 0. )THEN
    pchst = 0.
  ELSE
    pchst = SIGN( 1., arg1 ) * SIGN( 1., arg2 )
  END IF

END FUNCTION pchst
