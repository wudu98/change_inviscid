subroutine fh_n_weno3_f(f0,f1,f2,fhn)	  ! negatively biased
  implicit none
  integer ::  i
  real*8  ::  f0,f1,f2,fhn 
  real*8  ::  w(0:1),beta(0:1),alfa(0:1),sumalfa
  real*8  ::  fh(0:1)
  real*8   ::  d(0:1)=(/1.0/3.0,2.0/3.0/)
  real,parameter :: epsilon=1e-6

  fh(0)=  1.5*f1  - 0.5*f2
  fh(1)=  0.5*f0  + 0.5*f1

  beta(0)=(f2-f1)**2
  beta(1)=(f1-f0)**2

  do i=0,1
    alfa(i)=d(i)/(epsilon+beta(i))**2
  end do
  sumalfa=alfa(0)+alfa(1)
  do i=0,1
  w(i)=alfa(i)/sumalfa
  end do

  fhn=0
  do i=0,1
   fhn=fhn+w(i)*fh(i)
  end do
 
return
end