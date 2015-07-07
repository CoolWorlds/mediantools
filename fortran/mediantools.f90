MODULE mediantools

use coolmath

 implicit none

 CONTAINS

 ! =============================================================================
 SUBROUTINE detrend(needs_sort,scale_errors,bin_in,n_in,x_in,n_out,x_out)

 implicit none

 INTEGER :: i, m
 INTEGER :: noutliers, nclean, nnorm
 INTEGER :: bin, midpoint
 LOGICAL, INTENT(IN) :: needs_sort, scale_errors
 INTEGER, INTENT(IN) :: bin_in, n_in
 REAL(8), DIMENSION(3,n_in), INTENT(IN) :: x_in
 REAL(8), DIMENSION(3,n_in) :: x
 REAL(8), DIMENSION(n_in) :: t
 LOGICAL, DIMENSION(n_in) :: x_outlier
 REAL(8), DIMENSION(:,:), ALLOCATABLE :: xclean, xnorm
 LOGICAL, DIMENSION(:), ALLOCATABLE :: xnorm_outlier
 INTEGER, INTENT(OUT) :: n_out
 REAL(8), DIMENSION(3,n_in), INTENT(OUT) :: x_out

 ! Do not allow even bin sizes
 IF( MOD(bin_in, 2) .EQ. 0 ) THEN
   bin = bin_in + 1
   write(*,*) 'Increasing bin size to ',bin
 ELSE
   bin = bin_in
 END IF

 ! Sort the array, if needed
 x(:,:) = x_in(:,:)
 IF( needs_sort ) THEN
   t(:) = x(1,:)
   call double_sort(t,x,n_in,3)
 END IF

 ! Find outliers
 call identify_outliers(scale_errors,.FALSE.,bin,n_in,x,x_outlier,noutliers)
 
 ! Define xclean
 nclean = n_in - noutliers
 ALLOCATE(xclean(3,nclean))
 m = 0
 DO i=1,n_in
   IF( x_outlier(i) ) THEN
     ! donothing
   ELSE
     m = m + 1
     xclean(:,m) = x(:,i)
   END IF
 END DO

 ! Execute median filter
 midpoint = mid_point(bin)
 nnorm = nclean+2-2*midpoint
 ALLOCATE(xnorm(3,nnorm))
 call median_filter(.FALSE.,bin,nclean,xclean,xnorm)

 ! Check for any remaining outliers
 ALLOCATE(xnorm_outlier(nnorm))
 call identify_outliers(scale_errors,.FALSE.,bin,nnorm,xnorm,xnorm_outlier,m)

 ! Define final xout array
 n_out = 0
 DO i=1,nnorm
   IF( xnorm_outlier(i) ) THEN
     ! donothing
   ELSE
     n_out = n_out + 1
     x_out(:,n_out) = xnorm(:,i)
   END IF
 END DO
 ! Fill the remaining entries, which are just dummy entries
 DO i=n_out+1,n_in
   x_out(:,i) = 0.0D0
 END DO

 END SUBROUTINE detrend
 ! =============================================================================

 ! =============================================================================
 SUBROUTINE identify_outliers(scale_errors,needs_sort,bin,n,x,outlier,noutliers)

 implicit none

 INTEGER :: i
 LOGICAL, INTENT(IN) :: scale_errors, needs_sort
 INTEGER, INTENT(IN) :: bin, n
 REAL(8), DIMENSION(3,n), INTENT(IN) :: x
 REAL(8), DIMENSION(2,n) :: x2, xmed
 REAL(8), DIMENSION(n) :: xres
 REAL(8) :: maderror, mederror, sigmafac, sigmas
 INTEGER, INTENT(OUT) :: noutliers
 LOGICAL, DIMENSION(n), INTENT(OUT) :: outlier
 REAL(8), PARAMETER :: roottwo = 1.414213562373095D0
 REAL(8), PARAMETER :: madfac = 1.482602218505602D0

 ! Define median model
 x2(1,:) = x(1,:)
 x2(2,:) = x(2,:)
 call median_model(needs_sort,bin,n,x2,xmed)
 DO i=1,n
   xres(i) = xmed(2,i) - x(2,i)
 END DO

 ! Scale errors
 IF( scale_errors ) THEN
   maderror = madfac*median_deviation(xres,n)
   mederror = median(x(3,:),n)
   sigmafac = maderror/mederror
   write(*,*) 'Scaling errors by ',sigmafac
 ELSE
   sigmafac = 1.0D0
 END IF

 ! Number of sigmas to clips
 sigmas = inverf(1.0D0-(1.0D0/n))

 ! Identify outliers
 noutliers = 0
 DO i=1,n
   IF( DABS(xres(i)) .GE. sigmas*sigmafac*x(3,i) ) THEN
     outlier(i) = .TRUE.
     noutliers = noutliers + 1
   ELSE
     outlier(i) = .FALSE.
   END IF
 END DO
 write(*,*) 'Number of outliers = ',noutliers,' [',&
            100.0D0*DBLE(noutliers)/DBLE(n),'%]'

 END SUBROUTINE identify_outliers
 ! =============================================================================

 ! =============================================================================
 REAL(8) FUNCTION median_deviation(x,n)

 implicit none

 INTEGER, INTENT(IN) :: n
 REAL(8), DIMENSION(n), INTENT(IN) :: x
 REAL(8), DIMENSION(n) :: res
 INTEGER :: i
 REAL(8) :: xmed

 xmed = median(x,n)
 DO i=1,n
   res(i) = DABS( x(i) - xmed )
 END DO
 median_deviation = median(res,n)

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 INTEGER FUNCTION mid_point(bin)

 implicit none

 INTEGER, INTENT(IN) :: bin

 mid_point = FLOOR( REAL(bin)/2.0D0 ) + 1

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 SUBROUTINE median_filter(needs_sort,bin,n,x,z)

 implicit none

 LOGICAL, INTENT(IN) :: needs_sort
 INTEGER, INTENT(IN) :: bin, n
 REAL(8), DIMENSION(3,n), INTENT(IN) :: x
 REAL(8), DIMENSION(3,n+1-bin), INTENT(OUT) :: z
 REAL(8), DIMENSION(2,n) :: x2
 REAL(8), DIMENSION(2,n+1-bin) :: y
 INTEGER :: i, midpoint

 x2(1,:) = x(1,:)
 x2(2,:) = x(2,:)
 call median_model(needs_sort,bin,n,x2,y)

 ! Mid-segment
 midpoint = mid_point(bin)
 DO i=1,n+2-2*midpoint
   z(1,i) = y(1,i+midpoint-1)
   z(2,i) = x(2,i+midpoint-1)/y(2,i+midpoint-1)
   z(3,i) = x(3,i+midpoint-1)/y(2,i+midpoint-1)
 END DO

 END SUBROUTINE median_filter
 ! =============================================================================

 ! =============================================================================
 SUBROUTINE median_model(needs_sort,bin,n,xin,y)

 implicit none

 LOGICAL, INTENT(IN) :: needs_sort
 INTEGER, INTENT(IN) :: bin, n
 REAL(8), DIMENSION(2,n), INTENT(IN) :: xin
 REAL(8), DIMENSION(2,n), INTENT(OUT) :: y
 REAL(8), DIMENSION(2,n) :: x
 REAL(8), DIMENSION(n) :: t
 INTEGER :: i, j, jstart, jend, jmax
 INTEGER :: midpoint
 REAL(8), DIMENSION(bin) :: t_temp, f_temp
 REAL(8), DIMENSION(n+1-bin) :: tmed, fmed

 x(:,:) = xin(:,:)

 ! Check if sorting needed
 IF( needs_sort ) THEN
   t(:) = x(1,:)
   call double_sort(t,x,n,2)
 END IF

 ! Compute median model
 jmax = n + 1 - bin
 midpoint = mid_point(bin)
 DO j=1,jmax
   jstart = j
   jend = j + bin - 1
   t_temp = x(1,jstart:jend)
   f_temp = x(2,jstart:jend)
   tmed(j) = t_temp(midpoint)
   fmed(j) = median(f_temp,bin)
 END DO

 ! First-segment
 DO i=1,midpoint-1
   y(:,i) = x(:,i)
 END DO
 ! Mid-segment
 DO i=1,n+2-2*midpoint
   y(1,i+midpoint-1) = tmed(i)
   y(2,i+midpoint-1) = fmed(i)
 END DO
 ! Last-segment
 DO i=n-midpoint+2,n
   y(:,i) = x(:,i)
 END DO

 END SUBROUTINE median_model
 ! =============================================================================

 ! =============================================================================
 REAL(8) FUNCTION median(xin,n)

 ! Find the median of X(1), ... , X(N), using as much of the quicksort
 ! algorithm as is needed to isolate it.
 ! N.B. On exit, the array X is partially ordered.
 ! Latest revision - 26 November 1996

 implicit none

 INTEGER, INTENT(IN) :: n
 REAL(8), INTENT(IN), DIMENSION(n) :: xin

 ! Local variables
 REAL(8), DIMENSION(n) :: x
 REAL(8) :: temp, xhi, xlo, xmax, xmin
 LOGICAL :: odd
 INTEGER :: hi, lo, nby2, nby2p1, mid, i, j, k

 x(:) = xin(:)

 nby2 = n / 2
 nby2p1 = nby2 + 1
 odd = .true.

 ! HI & LO are position limits encompassing the median.
 IF (n == 2 * nby2) odd = .false.
   lo = 1
   hi = n
 IF (n < 3) THEN
   IF (n < 1) THEN
      median = 0.0
      RETURN
   END IF
   median = x(1)
   IF (n == 1) RETURN
   median = 0.5*(median + x(2))
   RETURN
 END IF

 ! Find median of 1st, middle & last values.
 10 mid = (lo + hi)/2
 median = x(mid)
 xlo = x(lo)
 xhi = x(hi)
 IF (xhi < xlo) THEN          ! Swap xhi & xlo
   temp = xhi
   xhi = xlo
   xlo = temp
 END IF
 IF (median > xhi) THEN
   median = xhi
 ELSE IF (median < xlo) THEN
   median = xlo
 END IF

 ! The basic quicksort algorithm to move all values <= the sort key (median)
 ! to the left-hand end, and all higher values to the other end.
 i = lo
 j = hi
 50 DO
    IF (x(i) >= median) EXIT
      i = i + 1
    END DO
    DO
      IF (x(j) <= median) EXIT
      j = j - 1
    END DO
 IF (i < j) THEN
   temp = x(i)
   x(i) = x(j)
   x(j) = temp
   i = i + 1
   j = j - 1
  ! Decide which half the median is in.
  IF (i <= j) GO TO 50
 END IF

 IF (.NOT. odd) THEN
   IF (j == nby2 .AND. i == nby2p1) GO TO 130
   IF (j < nby2) lo = i
   IF (i > nby2p1) hi = j
   IF (i /= j) GO TO 100
   IF (i == nby2) lo = nby2
   IF (j == nby2p1) hi = nby2p1
 ELSE
   IF (j < nby2p1) lo = i
   IF (i > nby2p1) hi = j
   IF (i /= j) GO TO 100

   ! Test whether median has been isolated.
   IF (i == nby2p1) RETURN
 END IF
 100 IF (lo < hi - 1) GO TO 10

 IF (.NOT. odd) THEN
   median = 0.5*(x(nby2) + x(nby2p1))
   RETURN
 END IF
 temp = x(lo)
 IF (temp > x(hi)) THEN
   x(lo) = x(hi)
   x(hi) = temp
 END IF
 median = x(nby2p1)
 RETURN

 ! Special case, N even, J = N/2 & I = J + 1, so the median is
 ! between the two halves of the series.   Find max. of the first
 ! half & min. of the second half, then average.

 130 xmax = x(1)
 DO k = lo, j
   xmax = MAX(xmax, x(k))
 END DO
 xmin = x(n)
 DO k = i, hi
   xmin = MIN(xmin, x(k))
 END DO
 median = 0.5*(xmin + xmax)

 RETURN

 END FUNCTION
 ! =============================================================================

 ! =============================================================================
 SUBROUTINE double_sort(Xsort,Rsort,length,param)

 implicit none

 INTEGER, INTENT(IN) :: length, param
 INTEGER :: j
 REAL(8), DIMENSION(length), INTENT(INOUT) :: Xsort
 REAL(8), DIMENSION(param,length), INTENT(INOUT) :: Rsort
 INTEGER :: swaps_made, counts
 REAL(8) :: tempX
 REAL(8), DIMENSION(param) :: tempR

 !write(*,*) 'Xsort = ',Xsort
 DO ! Repeat this loop until we break out
  swaps_made=0  ! Initially, we've made no swaps
  ! Make one pass of the bubble sort algorithm
  DO counts=1,(length-1)
   ! If item is greater than the one after it, then we initiate a swap
   IF( Xsort(counts) .GT. Xsort(counts+1) ) THEN
    ! Define initial temp
    tempX = Xsort(counts)
    DO j=1,param
     tempR(j) = Rsort(j,counts)
    END DO
    ! Displace by one
    Xsort(counts) = Xsort(counts+1)
    DO j=1,param
     Rsort(j,counts) = Rsort(j,counts+1)
    END DO
    ! Finalize swap
    Xsort(counts+1) = tempX
    DO j=1,param
     Rsort(j,counts+1) = tempR(j)
    END DO
    ! Increase the swaps_made count
    swaps_made = swaps_made+1
   END IF
  END DO
  ! If no swaps, break loop
  IF( swaps_made == 0 ) exit
 END DO

 END SUBROUTINE double_sort
 ! =============================================================================

END MODULE mediantools
