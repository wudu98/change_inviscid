! subroutine inviscid(Lu1,Lu2,Lu3,Lu4,Lu5)
! !-----------------------------------------------------------------------------
! 	implicit none
! 	include "Declaration.f90"
! 	integer :: i,k,j,nn,n0
!     real*8  :: k1,k2,k3,dx
! 	real*8 :: dF1(na,nc,nr), dF2(na,nc,nr), dF3(na,nc,nr), dF4(na,nc,nr), dF5(na,nc,nr)
! 	real*8  :: Lu1(na,nc,nr),  Lu2(na,nc,nr),  Lu3(na,nc,nr),  Lu4(na,nc,nr),  Lu5(na,nc,nr)    !��ճ��ķ���ֵ
!     real*8,allocatable :: fp(:,:),fn(:,:),fpx(:,:),fnx(:,:)
! 	real*8 ::  ddda,   dpda,    duada,   ducda,   durda !������߽����
! 	real*8 ::  dddr,   dpdr,    duadr,   ducdr,   durdr !�������	
! 	real*8 :: phi1,phi2,phi3,phi4,phi5,d1,d2,d3,d4,d5
! 	integer :: case_LB,case_UB
! !***************************!***************************
! !������Ϊ��:��ȡ��һ��������ƽ�е�������,ȷ�������ұ߽��λ��n0��NN,
! !���жϳ���߽�λ���Ƿ�Ϊ�����߽�,����ǵĻ���ȡΪcase_LB����case_UB����1
! !Ȼ�����dfdx.f90�Ƚ���ʸͨ���ķ���,�ٵ���WENO��ʽ�����ɢ����
! dx=da
! k2=0
! k3=0
! nn=na
! allocate(fp(1-BE:na+BE,5),fn(1-BE:na+BE,5),Fpx(na,5),Fnx(na,5))
! do k=1,nr
!    do j=1,nc	 
! 	 if(nr0+k<nz2)then
! 			 if(nx1>na0.and.nx1<=na0+na)then
! 				do i=1,na
! 				  if(na0+i==nx1)n0=i
! 				end do
! 				case_lb=1
! 				if(na0+na==nx)then
! 				  case_ub=1
! 				else
! 				  case_ub=0			  
! 				endif
! 			 elseif(na0>nx1)then
! 			   n0=1
! 			   case_lb=0
! 			   if(na0+na==nx)then		  
! 				 case_ub=1
! 			   else
! 				 case_ub=0         !1�����������߽��ڴ�		  
! 			   endif
! 			 endif
! 	 else !(nr0+k>=nz2)then
! 			 n0=1		 
! 			 if(na0==0)then
! 			   case_lb=1
! 			   if(na0+na==nx)then 
! 				 case_ub=1
! 			   else
! 				 case_ub=0         !1�����������߽��ڴ�
! 			   endif
! 			 else
! 			   case_lb=0
! 			   if(na0+na==nx)then		  
! 				 case_ub=1
! 			   else
! 				 case_ub=0         !1�����������߽��ڴ�		  
! 			   endif		   
! 			 endif	   
! 	 endif
! 	 include "dfdx.f90"
!   end do
! enddo
! deallocate(fp,fn,Fpx,Fnx)
! !------------------------------------------------------------------------------------------------------------------------------------------

! if(na0==0)then	
! i=1
! do k=1,nr
!   if(nr0+k>=nz2)then
!   do j=1,nc
!      dpda=(-3.0d0* p(i,j,k)+4.0d0* p(i+1,j,k)- p(i+2,j,k))/(2.0d0*da)
! 	 ddda=(-3.0d0* d(i,j,k)+4.0d0* d(i+1,j,k)- d(i+2,j,k))/(2.0d0*da)
! 	duada=(-3.0d0*ua(i,j,k)+4.0d0*ua(i+1,j,k)-ua(i+2,j,k))/(2.0d0*da)
! 	ducda=(-3.0d0*uc(i,j,k)+4.0d0*uc(i+1,j,k)-uc(i+2,j,k))/(2.0d0*da)
! 	durda=(-3.0d0*ur(i,j,k)+4.0d0*ur(i+1,j,k)-ur(i+2,j,k))/(2.0d0*da)
! 	if((ua(i,j,k)-sv(i,j,k))<0.0)then
!         phi1=(dadx(i)*ua(i,j,k)-sv(i,j,k)*dadx(i))*(-d(i,j,k)*sv(i,j,k)*duada+dpda)
!     else
! 	    phi1=0.0d0
! 	endif
! 	if(ua(i,j,k)<0.0)then
! 	    phi2= dadx(i)*ua(i,j,k)*(sv(i,j,k)**2*ddda-dpda)
!         phi3= dadx(i)*ua(i,j,k)*(dadx(i)*ducda)
!         phi4= dadx(i)*ua(i,j,k)*(dadx(i)*durda)
! 	else
! 	    phi2=0.0d0
! 		phi3=0.0d0
! 		phi4=0.0d0
! 	endif
! 	if((ua(i,j,k)+sv(i,j,k))<0.0)then
!         phi5=(dadx(i)*ua(i,j,k)+sv(i,j,k)*dadx(i))*( d(i,j,k)*sv(i,j,k)*duada+dpda)
! 	else
! 	    phi5=0.0d0
! 	endif
! !*****************************************
! 	d1=(0.5d0*(phi1+phi5)+phi2)/sv(i,j,k)**2
! 	d2=(phi5-phi1)/(2.0d0*d(i,j,k)*sv(i,j,k))
! 	d3=phi3/dadx(i)
! 	d4=phi4/dadx(i)
! 	d5=0.5d0*(phi1+phi5)
! 	df1(i,j,k)=(          d1            )*jj(i,j,k)
! 	df2(i,j,k)=(ua(i,j,k)*d1+d(i,j,k)*d2)*jj(i,j,k)
! 	df3(i,j,k)=(uc(i,j,k)*d1+d(i,j,k)*d3)*jj(i,j,k)
! 	df4(i,j,k)=(ur(i,j,k)*d1+d(i,j,k)*d4)*jj(i,j,k)
! 	df5(i,j,k)=( 0.5d0*(ua(i,j,k)**2+uc(i,j,k)**2+ur(i,j,k)**2)*d1 &
! 	                                      +d(i,j,k)*ua(i,j,k)*d2 &
! 										  +d(i,j,k)*uc(i,j,k)*d3 &
! 										  +d(i,j,k)*ur(i,j,k)*d4 &
! 										  +               2.5d0*d5 )*jj(i,j,k)
!   end do
!   endif
! end do
! endif
! !**********************************
! !Զ���߽�  �����һ������(na0+na==nx)�����һ�������(i=na)
! if(na0+na==nx)then	
! i=na
! do k=1,nr
!   do j=1,nc
!      dpda=(3.0d0* p(i,j,k)-4.0d0* p(i-1,j,k)+ p(i-2,j,k))/(2.0d0*da)
! 	 ddda=(3.0d0* d(i,j,k)-4.0d0* d(i-1,j,k)+ d(i-2,j,k))/(2.0d0*da)
! 	duada=(3.0d0*ua(i,j,k)-4.0d0*ua(i-1,j,k)+ua(i-2,j,k))/(2.0d0*da)
! 	ducda=(3.0d0*uc(i,j,k)-4.0d0*uc(i-1,j,k)+uc(i-2,j,k))/(2.0d0*da)
! 	durda=(3.0d0*ur(i,j,k)-4.0d0*ur(i-1,j,k)+ur(i-2,j,k))/(2.0d0*da)
! 	if((ua(i,j,k)-sv(i,j,k))<=0.0)then
!                 phi1=0.0d0
!         else
! 	    phi1=(dadx(i)*ua(i,j,k)-sv(i,j,k)*dadx(i))*(-d(i,j,k)*sv(i,j,k)*duada+dpda)
! 	endif
! 	if(ua(i,j,k)<=0.0)then
!                 phi2=0.0d0
!         	phi3=0.0d0
! 		phi4=0.0d0
! 	else
! 	    phi2= dadx(i)*ua(i,j,k)*(sv(i,j,k)**2*ddda-dpda)
!             phi3= dadx(i)*ua(i,j,k)*(dadx(i)*ducda)
!             phi4= dadx(i)*ua(i,j,k)*(dadx(i)*durda)
! 	endif
! 	if((ua(i,j,k)+sv(i,j,k))<=0.0)then
!         phi5=0.0d0
! 	else
! 	    phi5=(dadx(i)*ua(i,j,k)+sv(i,j,k)*dadx(i))*( d(i,j,k)*sv(i,j,k)*duada+dpda)
! 	endif
! !*****************************************
! 	d1=(0.5d0*(phi1+phi5)+phi2)/sv(i,j,k)**2
! 	d2=(phi5-phi1)/(2.0d0*d(i,j,k)*sv(i,j,k))
! 	d3=phi3/dadx(i)
! 	d4=phi4/dadx(i)
! 	d5=0.5d0*(phi1+phi5)
! 	df1(i,j,k)=(          d1            )*jj(i,j,k)
! 	df2(i,j,k)=(ua(i,j,k)*d1+d(i,j,k)*d2)*jj(i,j,k)
! 	df3(i,j,k)=(uc(i,j,k)*d1+d(i,j,k)*d3)*jj(i,j,k)
! 	df4(i,j,k)=(ur(i,j,k)*d1+d(i,j,k)*d4)*jj(i,j,k)
! 	df5(i,j,k)=( 0.5d0*(ua(i,j,k)**2+uc(i,j,k)**2+ur(i,j,k)**2)*d1 &
! 	                                      +d(i,j,k)*ua(i,j,k)*d2 &
! 										  +d(i,j,k)*uc(i,j,k)*d3 &
! 										  +d(i,j,k)*ur(i,j,k)*d4 &
! 										  +               2.5d0*d5 )*jj(i,j,k)
!   end do
! end do
!  endif
! !**********************************

!     do k=1,nr
! 	  do j=1,nc
! 		do i=1,na
! 		    if(    (na0+i>=nx1.and.nr0+k>1).or.(nr0+k>=nz2)      )then
! 			    Lu1(i,j,k)=df1(i,j,k)
! 				Lu2(i,j,k)=df2(i,j,k)
! 				Lu3(i,j,k)=df3(i,j,k)
! 				Lu4(i,j,k)=df4(i,j,k)
! 				Lu5(i,j,k)=df5(i,j,k)
! 			endif
! 		end do
! 	  end do
! 	end do

! !**************************************************************************************************

! k1=0
! k3=0
! n0=1
! nn=nc
! case_lb=0
! case_ub=0
! allocate(fp(1-BE:nc+BE,5),fn(1-BE:nc+BE,5),Fpx(nc,5),Fnx(nc,5))
! do k=1,nr
! dx=z(k)*dc
! do i=1,na
!        if(    (na0+i>=nx1.and.nr0+k>1).or.(nr0+k>=nz2)      )then
! 		  do j=1-BE,nc+BE
! 		    k2=dcdy(j)*jj(i,j,k)
! 		    call steger_warming(k1,k2,k3,d(i,j,k),ua(i,j,k),uc(i,j,k),ur(i,j,k),sv(i,j,k),Fp(j,:),Fn(j,:))
! 		  end do
!   	      call dfdx_p(case_lb,case_ub,n0,nn,dx,fp,fpx)
!   		  call dfdx_n(case_lb,case_ub,n0,nn,dx,fn,fnx)
!           do j=n0,nn
! 		    df1(i,j,k)=fpx(j,1)+fnx(j,1)
! 			df2(i,j,k)=fpx(j,2)+fnx(j,2)
! 			df3(i,j,k)=fpx(j,3)+fnx(j,3)
! 			df4(i,j,k)=fpx(j,4)+fnx(j,4)
! 			df5(i,j,k)=fpx(j,5)+fnx(j,5)
!           end do
! 	    endif
! end do
! end do
! deallocate(fp,fn,Fpx,Fnx)

! !**************************************
!     do k=1,nr
! 	  do j=1,nc
! 		do i=1,na
! 		    if(    (na0+i>=nx1.and.nr0+k>1).or.(nr0+k>=nz2)      )then
! 			    Lu1(i,j,k)=Lu1(i,j,k)+df1(i,j,k)
! 				Lu2(i,j,k)=Lu2(i,j,k)+df2(i,j,k)
! 				Lu3(i,j,k)=Lu3(i,j,k)+df3(i,j,k)
! 				Lu4(i,j,k)=Lu4(i,j,k)+df4(i,j,k)
! 				Lu5(i,j,k)=Lu5(i,j,k)+df5(i,j,k)
! 			endif
! 		end do
! 	  end do
! 	end do
! !**************************************************************************************************
! !------------------------------------------------------------------------------------------------------------------------------------------
! dx=dr
! k1=0
! k2=0
! nn=nr
! allocate(fp(1-BE:nr+BE,5),fn(1-BE:nr+BE,5),Fpx(nr,5),Fnx(nr,5))
! do j=1,nc
!   do i=1,na
!     if(na0+i<nx1)then
! 		  if(nz2>nr0.and.nz2<=nr0+nr)then
! 			do k=1,nr
! 			  if(nr0+k==nz2)n0=k
! 			end do
! 			case_lb=1
! 			if(nr0+nr==nz)then		  
! 			  case_ub=1
! 			else
! 			  case_ub=0
! 			endif		
! 		  elseif(nr0>nz2)then
! 			n0=1
! 			case_lb=0
! 			if(nr0+nr==nz)then		  
! 			  case_ub=1
! 			else
! 			  case_ub=0         !1�����������߽��ڴ�
! 			endif		
! 		  endif
! 	else
! 		  n0=1
! 		  if(nr0==0)then
! 			case_lb=1
! 			if(nr0+nr==nz)then		  
! 			  case_ub=1
! 			else
! 			  case_ub=0         !1�����������߽��ڴ�
! 			endif		
! 		  else
! 			case_lb=0
! 			if(nr0+nr==nz)then		  
! 			  case_ub=1
! 			else
! 			  case_ub=0         !1�����������߽��ڴ�
! 			endif					
! 		  endif	  	  
! 	endif
! 	include "dfdz.f90"
!   end do
! enddo
! deallocate(fp,fn,Fpx,Fnx)
! !------------------------------------------------------------------------------------------------------------------------------------------
! !���ڱ߽�
! if(nr0+nr==nz)then
! k=nr
! do j=1,nc
!   do i=1,na
!          dddr=(3.0d0* d(i,j,k)-4.0d0* d(i,j,k-1)+ d(i,j,k-2))/(2.0d0*dr)
! 		 dpdr=(3.0d0* p(i,j,k)-4.0d0* p(i,j,k-1)+ p(i,j,k-2))/(2.0d0*dr)
! 		duadr=(3.0d0*ua(i,j,k)-4.0d0*ua(i,j,k-1)+ua(i,j,k-2))/(2.0d0*dr)
! 		ducdr=(3.0d0*uc(i,j,k)-4.0d0*uc(i,j,k-1)+uc(i,j,k-2))/(2.0d0*dr)
!         durdr=(3.0d0*ur(i,j,k)-4.0d0*ur(i,j,k-1)+ur(i,j,k-2))/(2.0d0*dr)
! 		if((ur(i,j,k)-sv(i,j,k))<=0.0)then
!             phi1=0.0d0
! 		else
! 		    phi1=(drdz(k)*ur(i,j,k)-sv(i,j,k)*drdz(k))*(-d(i,j,k)*sv(i,j,k)*durdr+dpdr)
! 		endif

!         if(ur(i,j,k)<=0.0)then
! 		    phi2=0.0d0
! 			phi3=0.0d0
! 			phi4=0.0d0
! 		else
! 			phi2=drdz(k)*ur(i,j,k)*(sv(i,j,k)**2*dddr-dpdr)
! 			phi3=drdz(k)*ur(i,j,k)*(drdz(k)*duadr)
! 			phi4=drdz(k)*ur(i,j,k)*(drdz(k)*ducdr)
! 		endif
! 		if((ur(i,j,k)+sv(i,j,k))<=0.0)then
! 		    phi5=0.0d0
! 		else
! 	        phi5=(drdz(k)*ur(i,j,k)+sv(i,j,k)*drdz(k))*(d(i,j,k)*sv(i,j,k)*durdr+dpdr)
! 		endif	
! !*****************************************
! 	d1=(0.5d0*(phi1+phi5)+phi2)/sv(i,j,k)**2
! 	d2=phi3/drdz(k)
! 	d3=phi4/drdz(k)
! 	d4=(phi5-phi1)/(2.0d0*d(i,j,k)*sv(i,j,k))
! 	d5=0.5d0*(phi1+phi5)
! 	df1(i,j,k)=(          d1            )*jj(i,j,k)
! 	df2(i,j,k)=(ua(i,j,k)*d1+d(i,j,k)*d2)*jj(i,j,k)
! 	df3(i,j,k)=(uc(i,j,k)*d1+d(i,j,k)*d3)*jj(i,j,k)
! 	df4(i,j,k)=(ur(i,j,k)*d1+d(i,j,k)*d4)*jj(i,j,k)
! 	df5(i,j,k)=( 0.5d0*(ua(i,j,k)**2+uc(i,j,k)**2+ur(i,j,k)**2)*d1 &
! 	                                      +d(i,j,k)*ua(i,j,k)*d2 &
! 										  +d(i,j,k)*uc(i,j,k)*d3 &
! 										  +d(i,j,k)*ur(i,j,k)*d4 &
! 										  +               2.5d0*d5 )*jj(i,j,k)	
! !*****************************************
!   end do
! end do
! endif
! !**************************************************************************************************

! !**************************************************************************************************
! !��ճ����ķ���ֵ
!     do k=1,nr
! 	  do j=1,nc
! 		do i=1,na
! 		    if(    (na0+i>=nx1.and.nr0+k>1).or.(nr0+k>=nz2)      )then
! 			    Lu1(i,j,k)=(Lu1(i,j,k)+df1(i,j,k))/jj(i,j,k)
! 				Lu2(i,j,k)=(Lu2(i,j,k)+df2(i,j,k))/jj(i,j,k)
! 				Lu3(i,j,k)=(Lu3(i,j,k)+df3(i,j,k))/jj(i,j,k)
! 				Lu4(i,j,k)=(Lu4(i,j,k)+df4(i,j,k))/jj(i,j,k)
! 				Lu5(i,j,k)=(Lu5(i,j,k)+df5(i,j,k))/jj(i,j,k)
! 			endif
! 		end do
! 	  end do
! 	end do

! !**************************************************************************************************
!   return
! end
! !------------------------------------------------------------------------------------------------------------------------------------------

!**********************************************************************************************************************************
subroutine steger_warming(k1,k2,k3,d,u,v,w,sv,Fp,Fn)
    implicit none
	real*8  :: k1,k2,k3
	real*8  :: k1b,k2b,k3b,kk
	real*8  :: d,u,v,w,sv                   !sv:sound velocity
	real*8  :: lambda(3),lamb(3)!lap(3),lan(3)   !��������ֵ
	real*8  :: Fp(5),Fn(5)
    integer :: i
	real*8,parameter :: gama   =1.4d0   !gama=1.4
	real*8,parameter :: gama2_2=0.8d0   !2(r-1)
	real*8,parameter :: gama2  =2.8d0   !2r
	real*8,parameter :: gama_1 =0.4d0   !r-1
    real*8,parameter :: gama_w =2.0d0   !(3-r)/2(r-1)
	real*8,parameter :: epsilon_2=1.0D-12    !epsilon**2

	kk=dsqrt(k1*k1+k2*k2+k3*k3)
	lambda(1)=k1*u+k2*v+k3*w
	lambda(2)=lambda(1)+sv*kk
	lambda(3)=lambda(1)-sv*kk              !����任�µ�����ֵ
    k1b=k1/kk
	k2b=k2/kk
	k3b=k3/kk
		
    do i=1,3
	    lamb(i)=0.5d0*(lambda(i)+(lambda(i)**2+epsilon_2)**0.5) !��lambplus
	end do

	Fp(1)=d/gama2*(gama2_2*lamb(1)  +lamb(2)           +lamb(3)           )
	Fp(2)=d/gama2*(gama2_2*lamb(1)*u+lamb(2)*(u+sv*k1b)+lamb(3)*(u-sv*k1b))
	Fp(3)=d/gama2*(gama2_2*lamb(1)*v+lamb(2)*(v+sv*k2b)+lamb(3)*(v-sv*k2b))
	Fp(4)=d/gama2*(gama2_2*lamb(1)*w+lamb(2)*(w+sv*k3b)+lamb(3)*(w-sv*k3b))
	Fp(5)=d/gama2*(gama_1*lamb(1)*(u**2+v**2+w**2)+&
		           0.5*lamb(2)*((u+sv*k1b)**2+(v+sv*k2b)**2+(w+sv*k3b)**2)+&
		           0.5*lamb(3)*((u-sv*k1b)**2+(v-sv*k2b)**2+(w-sv*k3b)**2)+&
				   gama_w*(lamb(2)+lamb(3))*sv**2        )
						   !  +&
						!d*gama2_2*lamb(1)*k1b*(k2b*w-k3b*v)     ) gama2:2*gama ͨ������k1b,k2b,k3b������G,H��ʧͨ������
    do i=1,3
	    lamb(i)=0.5d0*(lambda(i)-(lambda(i)**2+epsilon_2)**0.5)  !��lambnegative
	end do

	Fn(1)=d/gama2*(gama2_2*lamb(1)  +lamb(2)           +lamb(3)           )  !2(r-1) lamb(1)=u lamb(2)=u+a lamb(3)=u-a
	Fn(2)=d/gama2*(gama2_2*lamb(1)*u+lamb(2)*(u+sv*k1b)+lamb(3)*(u-sv*k1b))
	Fn(3)=d/gama2*(gama2_2*lamb(1)*v+lamb(2)*(v+sv*k2b)+lamb(3)*(v-sv*k2b))
	Fn(4)=d/gama2*(gama2_2*lamb(1)*w+lamb(2)*(w+sv*k3b)+lamb(3)*(w-sv*k3b))
	Fn(5)=d/gama2*(gama_1*lamb(1)*(u**2+v**2+w**2)+&
		           0.5*lamb(2)*((u+sv*k1b)**2+(v+sv*k2b)**2+(w+sv*k3b)**2)+&
		           0.5*lamb(3)*((u-sv*k1b)**2+(v-sv*k2b)**2+(w-sv*k3b)**2)+&    !H=0.5*(u**2+v**2+w**2)+a**2/(r-1)
				   gama_w*(lamb(2)+lamb(3))*sv**2        )
						   !  +&
						!d*gama2_2*lamb(1)*k1b*(k2b*w-k3b*v)     )
  return
end
!**********************************************************************************************************************************
!*********************************************************************************************************************************
subroutine dfdx_p(case_lb,case_ub,jj,nn,dx,fp,fpx)
    integer :: case_lb,case_ub
	integer :: jj,nn
	real*8  :: dx,fp(-5:nn+6,5),fpx(1:nn,5),fhp(0:nn,5)
	integer :: i,m

do m=1,5
		do i=jj+3,nn-3
	      call fh_p_weno7(fp(i-3,m),fp(i-2,m),fp(i-1,m),fp(i,m),fp(i+1,m),fp(i+2,m),fp(i+3,m),fhp(i,m))
!          call fh_p_weno5(fp(i-2,m),fp(i-1,m),fp(i,m),fp(i+1,m),fp(i+2,m),fhp(i,m))
        end do
		do i=jj+4,nn-3
		  fpx(i,m)=(fhp(i,m)-fhp(i-1,m))/dx
		end do
        if(case_lb==1)then  !���б߽�
          i=jj+2
		  call fh_p_weno5(fp(i-2,m),fp(i-1,m),fp(i,m),fp(i+1,m),fp(i+2,m),fhp(i,m))
		  fpx(jj+3,m)=(fhp(jj+3,m)-fhp(jj+2,m))/dx
		  i=jj+1
		  call fh_p_weno3(fp(i-1,m),fp(i,m),fp(i+1,m),fhp(i,m))
		  fpx(i+1,m)=(fhp(i+1,m)-fhp(i,m))/dx
          fpx(jj,  m)=(fp(jj+1,m)-fp(jj,m))/dx
		  fpx(jj+1,m)=(fp(jj+1,m)-fp(jj,m))/dx
!		  fpx(jj+2,m)=(fp(jj,m)-6.0*fp(jj+1,m)+3.0*fp(jj+2,m)+2.0*fp(jj+3,m))/(6.0*dx)
		else
		  do i=jj+2,jj-1,-1
		    call fh_p_weno7(fp(i-3,m),fp(i-2,m),fp(i-1,m),fp(i,m),fp(i+1,m),fp(i+2,m),fp(i+3,m),fhp(i,m))
!		    call fh_p_weno5(fp(i-2,m),fp(i-1,m),fp(i,m),fp(i+1,m),fp(i+2,m),fhp(i,m))
		    fpx(i+1,m)=(fhp(i+1,m)-fhp(i,m))/dx
          end do
		endif
		
		if(case_ub==1)then
		  fpx(nn,m)=(fp(nn,m)-fp(nn-1,m))/dx
          i=nn-2
		  call fh_p_weno5(fp(i-2,m),fp(i-1,m),fp(i,m),fp(i+1,m),fp(i+2,m),fhp(i,m))
		  fpx(i,m)=(fhp(i,m)-fhp(i-1,m))/dx
		  i=nn-1
		  call fh_p_weno3(fp(i-1,m),fp(i,m),fp(i+1,m),fhp(i,m))
		  fpx(i,m)=(fhp(i,m)-fhp(i-1,m))/dx
!		  fpx(nn-1,m)=(fp(nn-3,m)-6.0*fp(nn-2,m)+3.0*fp(nn-1,m)+2.0*fp(nn,m))/(6.0*dx)
		else
		  do i=nn-2,nn
		    call fh_p_weno7(fp(i-3,m),fp(i-2,m),fp(i-1,m),fp(i,m),fp(i+1,m),fp(i+2,m),fp(i+3,m),fhp(i,m))
!            call fh_p_weno5(fp(i-2,m),fp(i-1,m),fp(i,m),fp(i+1,m),fp(i+2,m),fhp(i,m))
			fpx(i,m)=(fhp(i,m)-fhp(i-1,m))/dx
		  end do
		  
		endif
enddo		

  return	  
end




subroutine dfdx_n(case_lb,case_ub,jj,nn,dx,fn,fnx)
    integer :: case_lb,case_ub
	integer :: jj,nn
	real*8  :: dx,fn(-5:nn+6,5),fnx(1:nn,5),fhn(0:nn+1,5)
    integer :: i,m
do m=1,5
	    do i=jj+3,nn-3
	      call fh_n_weno7(fn(i-3,m),fn(i-2,m),fn(i-1,m),fn(i,m),fn(i+1,m),fn(i+2,m),fn(i+3,m),fhn(i,m)) !ע��˴�weno�ع���ģ����ʵ��iΪ���ĵ�
!        call fh_n_weno5(fn(i-2,m),fn(i-1,m),fn(i,m),fn(i+1,m),fn(i+2,m),fhn(i,m))
        enddo
        do i=jj+3,nn-4
		  fnx(i,m)=(fhn(i+1,m)-fhn(i,m))/dx
		end do
        if(case_lb==1)then  !���б߽�
		  i=jj+2
		  call fh_n_weno5(fn(i-2,m),fn(i-1,m),fn(i,m),fn(i+1,m),fn(i+2,m),fhn(i,m))
		  fnx(jj+2,m)=(fhn(jj+3,m)-fhn(jj+2,m))/dx
          i=jj+1
		  call fh_n_weno3(fn(i-1,m),fn(i,m),fn(i+1,m),fhn(i,m))
		  fnx(jj+1,m)=(fhn(jj+2,m)-fhn(jj+1,m))/dx
          i=jj
		  fnx(i,m)=(fn(i+1,m)-fn(i,m))/dx
!		  fnx(i,m)=(-2.0*fn(i-1,m)-3.0*fn(i,m)+6.0*fn(i+1,m)-fn(i+2,m))/(6.0*dx)
		else
		  do i=jj,jj+2
		    call fh_n_weno7(fn(i-3,m),fn(i-2,m),fn(i-1,m),fn(i,m),fn(i+1,m),fn(i+2,m),fn(i+3,m),fhn(i,m))
!          call fh_n_weno5(fn(i-2,m),fn(i-1,m),fn(i,m),fn(i+1,m),fn(i+2,m),fhn(i,m))
		  end do
          do i=jj,jj+2
            fnx(i,m)=(fhn(i+1,m)-fhn(i,m))/dx
          end do
        endif

		if(case_ub==1)then
          i=nn-2
	      call fh_n_weno5(fn(i-2,m),fn(i-1,m),fn(i,m),fn(i+1,m),fn(i+2,m),fhn(i,m))
		  fnx(nn-3,m)=(fhn(nn-2,m)-fhn(nn-3,m))/dx
		  i=nn-1
	      call fh_n_weno3(fn(i-1,m),fn(i,m),fn(i+1,m),fhn(i,m))
		  fnx(nn-2,m)=(fhn(nn-1,m)-fhn(nn-2,m))/dx
!		  fnx(i,m)=(-2.0*fn(i-1,m)-3.0*fn(i,m)+6.0*fn(i+1,m)-fn(i+2,m))/(6.0*dx)
		  i=nn-1
		  fnx(i,m)=(fn(nn,m)-fn(nn-1,m))/dx
		  i=nn
		  fnx(i,m)=(fn(nn,m)-fn(nn-1,m))/dx
		else
		  do i=nn-2,nn+1
            call fh_n_weno7(fn(i-3,m),fn(i-2,m),fn(i-1,m),fn(i,m),fn(i+1,m),fn(i+2,m),fn(i+3,m),fhn(i,m))
!		  call fh_n_weno5(fn(i-2,m),fn(i-1,m),fn(i,m),fn(i+1,m),fn(i+2,m),fhn(i,m))
		  end do
          do i=nn-3,nn
          fnx(i,m)=(fhn(i+1,m)-fhn(i,m))/dx
          end do
		endif
end do		
  return
end











   subroutine fh_p_weno3(f_1,f0,f1,fhp)   !positively biased
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
 
 
 
 
 subroutine fh_n_weno3(f0,f1,f2,fhn)	  ! negatively biased
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







subroutine fh_p_weno5(f_2,f_1,f0,f1,f2,fhp)   !positively biased
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


subroutine fh_n_weno5(f_1,f0,f1,f2,f3,fhn)	  ! negatively biased
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

subroutine fh_p_weno7(f_3,f_2,f_1,f0,f1,f2,f3,fhp)   !positively biased
  implicit none
  integer ::  i
  real*8  ::  f_3,f_2,f_1,f0,f1,f2,f3,fhp  
  real*8  ::  w(0:3),beta(0:3),alfa(0:3),sumalfa
  real*8  ::  fh(0:3)
  real*8   ::  d(0:3)
  real*8   ::  TV(0:3),TVR,TV_MAX,TV_MIN
  real*8,parameter :: epsilon=1e-6
!optimal weight
  d(0)=1.0d0/35.0d0
  d(1)=12.0d0/35.0d0
  d(2)=18.0d0/35.0d0
  d(3)=4.0d0/35.0d0

  fh(0)= -1.0d0/4.0d0*f_3 + 13.0d0/12.0d0*f_2 - 23.0d0/12.0d0*f_1 + 25.0d0/12.0d0*f0
  fh(1)= 1.0d0/12.0d0*f_2 - 5.0d0/12.0d0*f_1  + 13.0d0/12.0d0*f0  + 1.0d0/4.0d0*f1
  fh(2)=-1.0d0/12.0d0*f_1 + 7.0d0/12.0d0*f0   + 7.0d0/12.0d0*f1   - 1.0d0/12.0d0*f2
  fh(3)= 1.0d0/4.0d0*f0   + 13.0d0/12.0d0*f1  - 5.0d0/12.0d0*f2   + 1.0d0/12.0d0*f3

  TV(0)=abs(f_3-f_2)+abs(f_2-f_1)+abs(f_1-f0)
  TV(1)=abs(f_2-f_1)+abs(f_1-f0 )+abs(f0 -f1)
  TV(2)=abs(f_1-f0 )+abs(f0 -f1 )+abs(f1 -f2)
  TV(3)=abs(f0 -f1 )+abs(f1 -f2 )+abs(f2 -f3)

  TV_MAX=MaxVal(TV)
  TV_MIN=MinVal(TV)
  TVR=TV_MAX/(TV_MIN+epsilon)

  if( TV_MAX<0.2d0 .and. TVR<5.0d0) then
      do i=0,3
        w(i)=d(i)
      end do
  else
  beta(0)=f_3*(547.0*f_3-3882.0*f_2+4642.0*f_1-1854.0*f0)+f_2*(7043.0*f_2-17246.0*f_1+7042.0*f0)&
  &+f_1*(11003.0*f_1-9402.0*f0)+2107.0*f0**2
  beta(1)=f_2*(267.0*f_2-1642.0*f_1+1602.0*f0 -494.0*f1 )+f_1*(2843.0*f_1-5966.0*f0  +1922.0*f1)&
  &+f0*(3443.0*f0   -2522.0*f1)+547.0*f1**2
  beta(2)=f_1*(547.0*f_1-2522.0*f0 +1922.0*f1 -494.0*f2 )+f0*(3443.0*f0  -5966.0*f1  +1602.0*f2)&
  &+f1*(2843.0*f1   -1642.0*f2)+267.0*f2**2
  beta(3)=f0*(2107.0*f0 -9402.0*f1 +7042.0*f2 -1854.0*f3)+f1*(11003.0*f1 -17246.0*f2 +4642.0*f3)&
  &+f2*(7043.0*f2   -3882.0*f3)+547.0*f3**2

  do i=0,3
	  alfa(i)=d(i)/(epsilon+beta(i))**2
  end do

  sumalfa=alfa(0)+alfa(1)+alfa(2)+alfa(3)

  do i=0,3
	w(i)=alfa(i)/sumalfa
  end do

  endif

  fhp=0.0d0
  do i=0,3
    fhp=fhp+w(i)*fh(i)
  end do
 
 return
end

subroutine fh_n_weno7(f_3,f_2,f_1,f0,f1,f2,f3,fhn)    ! negatively biased
  implicit none
  integer ::  i
  real*8  ::  f_3,f_2,f_1,f0,f1,f2,f3,fhn 
  real*8  ::  w(0:3),beta(0:3),alfa(0:3),sumalfa
  real*8  ::  fh(0:3)
  real*8  ::  d(0:3)
  real*8  ::  TV(0:3),TVR,TV_MAX,TV_MIN
  real*8,parameter :: epsilon=1e-6

  d(0)=1.0/35.0
  d(1)=12.0/35.0
  d(2)=18.0/35.0
  d(3)=4.0/35.0

  fh(0)= -1.0/4.0*f3 + 13.0/12.0*f2 - 23.0/12.0*f1 + 25.0/12.0*f0
  fh(1)= 1.0/12.0*f2 - 5.0/12.0*f1  + 13.0/12.0*f0  + 1.0/4.0*f_1
  fh(2)=-1.0/12.0*f1 + 7.0/12.0*f0   + 7.0/12.0*f_1   - 1.0/12.0*f_2
  fh(3)= 1.0/4.0*f0   + 13.0/12.0*f_1  - 5.0/12.0*f_2   + 1.0/12.0*f_3

  TV(0)=abs(f0 -f1 )+abs(f1 -f2 )+abs(f2 -f3)
  TV(1)=abs(f_1-f0 )+abs(f0 -f1 )+abs(f1 -f2)
  TV(2)=abs(f_2-f_1)+abs(f_1-f0 )+abs(f0 -f1)
  TV(3)=abs(f_3-f_2)+abs(f_2-f_1)+abs(f_1-f0)

  TV_MAX=MaxVal(TV)
  TV_MIN=MinVal(TV)
  TVR=TV_MAX/(TV_MIN+epsilon)

  if( TV_MAX<0.2d0 .and. TVR<5.0d0) then
      do i=0,3
        w(i)=d(i)
      end do
  else
  beta(0)=f3*(547.0*f3-3882.0*f2+4642.0*f1-1854.0*f0)+f2*(7043.0*f2-17246.0*f1+7042.0*f0)&
  &+f1*(11003.0*f1-9402.0*f0)+2107.0*f0**2
  beta(1)=f2*(267.0*f2-1642.0*f1+1602.0*f0 -494.0*f_1 )+f1*(2843.0*f1-5966.0*f0  +1922.0*f_1)&
  &+f0*(3443.0*f0   -2522.0*f_1)+547.0*f_1**2
  beta(2)=f1*(547.0*f1-2522.0*f0 +1922.0*f_1 -494.0*f_2 )+f0*(3443.0*f0  -5966.0*f_1  +1602.0*f_2)&
  &+f_1*(2843.0*f_1   -1642.0*f_2)+267.0*f_2**2
  beta(3)=f0*(2107.0*f0 -9402.0*f_1 +7042.0*f_2 -1854.0*f_3)+f_1*(11003.0*f_1 -17246.0*f_2 +4642.0*f_3)&
  &+f_2*(7043.0*f_2   -3882.0*f_3)+547.0*f_3**2

  do i=0,3
    alfa(i)=d(i)/(epsilon+beta(i))**2
  end do

  sumalfa=alfa(0)+alfa(1)+alfa(2)+alfa(3)

  do i=0,3
  w(i)=alfa(i)/sumalfa
  end do
  endif

  fhn=0.0d0
  do i=0,3
    fhn=fhn+w(i)*fh(i)
  end do
 
 return
end

