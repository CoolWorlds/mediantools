PROGRAM example

use mediantools

 implicit none

 INTEGER :: i
 INTEGER, PARAMETER :: n = 60560
 INTEGER :: bin, nnice
 REAL(8), DIMENSION(3,n) :: x, xnice

 ! Read-in the raw.dat
 open(unit=10,file='raw.dat')
 DO i=1,n
   read(10,*) x(:,i)
 END DO
 close(10)

 bin = 51
 call detrend(.TRUE.,.TRUE.,bin,n,x,nnice,xnice)

 ! Write-out the raw.med.dat
 open(unit=11,file='raw.med.dat')
 DO i=1,nnice
   write(11,*) xnice(:,i)
 END DO
 close(11)

END PROGRAM example
