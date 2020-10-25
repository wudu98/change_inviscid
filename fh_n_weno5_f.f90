subroutine fh_n_weno5_f(f_1,f0,f1,f2,f3,fhn)	  ! negatively biased
  implicit none
  integer ::  i
  real*8  ::  f_1,f0,f1,f2,f3,fhn 
  real*8  ::  w(0:2),beta(0:2),alfa(0:2),sumalfa
  real*8  ::  fh(0:2)
  real*8    ::  d(0:2)=(/0.1,0.6,0.3/)
  real*8,parameter :: epsilon=1e-6

  fh(0)=11.0/6.0*f1  - 7.0/6.0*f2 + 1.0/3.0*f3
  fh(1)= 1.0/3.0*f0  + 5.0/6.0*f1 - 1.0/6.0*f2
  fh(2)=-1.0/6.0*f_1 + 5.0/6.0*f0 + 1.0/3.0*f1

  beta(0)=13.0/12.0*(f1 - 2.0*f2 + f3)**2 + 0.25*(3.0*f1  - 4.0*f2 +    f3)**2
  beta(1)=13.0/12.0*(f0 - 2.0*f1 + f2)**2 + 0.25*(    f0  -     f2        )**2
  beta(2)=13.0/12.0*(f_1- 2.0*f0 + f1)**2 + 0.25*(  f_1   - 4.0*f0+ 3.0*f1)**2

  do i=0,2
	  alfa(i)=d(i)/(epsilon+beta(i))**2
  end do
  sumalfa=alfa(0)+alfa(1)+alfa(2)
  do i=0,2
	w(i)=alfa(i)/sumalfa
  end do

  fhn=0
  do i=0,2
    fhn=fhn+w(i)*fh(i)
  end do
 
 return
end