subroutine fh_p_weno5_f(f_2,f_1,f0,f1,f2,fhp)   !positively biased
  implicit none
  integer ::  i
  real*8  ::  f_2,f_1,f0,f1,f2,fhp  
  real*8  ::  w(0:2),beta(0:2),alfa(0:2),sumalfa
  real*8  ::  fh(0:2)
  real*8  ::  d(0:2)=(/0.3,0.6,0.1/)
  real,parameter :: epsilon=1e-6

  fh(0)= 1.0/3.0*f0  + 5.0/6.0*f1  -  1.0/6.0*f2
  fh(1)=-1.0/6.0*f_1 + 5.0/6.0*f0  +  1.0/3.0*f1
  fh(2)= 1.0/3.0*f_2 - 7.0/6.0*f_1 + 11.0/6.0*f0

  beta(0)=13.0/12.0*(f0 - 2.0*f1 + f2)**2  + 0.25*(3.0*f0- 4.0*f1 +   f2)**2
  beta(1)=13.0/12.0*(f_1- 2.0*f0 + f1)**2  + 0.25*(f_1     -   f1       )**2
  beta(2)=13.0/12.0*(f_2- 2.0*f_1+ f0)**2  + 0.25*(f_2 - 4.0*f_1+ 3.0*f0)**2

  do i=0,2
	  alfa(i)=d(i)/(epsilon+beta(i))**2
  end do
  sumalfa=alfa(0)+alfa(1)+alfa(2)
  do i=0,2
	w(i)=alfa(i)/sumalfa
  end do

  fhp=0
  do i=0,2
    fhp=fhp+w(i)*fh(i)
  end do
 
 return
end