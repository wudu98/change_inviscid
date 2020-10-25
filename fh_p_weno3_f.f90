subroutine fh_p_weno3_f(f_1,f0,f1,fhp)   !positively biased
  implicit none
  integer ::  i
  real*8  ::  f_1,f0,f1,fhp  
  real*8  ::  w(0:1),beta(0:1),alfa(0:1),sumalfa
  real*8  ::  fh(0:1)
  real*8   ::  d(0:1)=(/2.0/3.0,1.0/3.0/)
  real,parameter :: epsilon=1e-6

  fh(0)= 0.5*f0  + 0.5*f1
  fh(1)=-0.5*f_1 + 1.5*f0

  beta(0)=(f1-f0)**2
  beta(1)=(f0-f_1)**2
  
  do i=0,1
    alfa(i)=d(i)/(epsilon+beta(i))**2
  end do
  sumalfa=alfa(0)+alfa(1)
  do i=0,1
  w(i)=alfa(i)/sumalfa
  end do

  fhp=0
  do i=0,1
    fhp=fhp+w(i)*fh(i)
  end do
 
 return
end