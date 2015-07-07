MODULE coolmath

 implicit none

 CONTAINS

 ! =============================================================================
 REAL(8) FUNCTION inverf(x)

 ! Calculate the inverse error function using:
 ! Sergei Winitzki, 6 February 2008, "A handy approximation for the error
 ! function and its inverse"
 !
 ! a_winitzki = [8/(3pi)]*[(pi-3)/(4-pi)]
 ! b_winitzki = 2/(pi*a)
 !
 ! Accuracy is < 4e-3

  implicit none

  REAL(8), INTENT(IN) :: x
  REAL(8), PARAMETER :: a_winitzki = 0.1400122886866666D0
  REAL(8), PARAMETER :: b_winitzki = 4.546884979448284D0
  REAL(8) :: y, logxsq
  INTEGER :: factor

  IF( x .LT. 0.0D0 ) THEN
    factor = -1
  ELSE
    factor = 1
  END IF

  logxsq = DLOG(1.0D0 - x**2)
  y = b_winitzki + 0.5D0*logxsq
  y = DSQRT( y**2 - (logxsq/a_winitzki) ) - y
  inverf = factor*DSQRT(y)

 END FUNCTION
 ! =============================================================================

END MODULE coolmath
