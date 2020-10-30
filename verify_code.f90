program verify_code
    implicit none 
    real*8  :: x
	real*8  :: k1=1.1,k2=2.2,k3=3.3
    real*8  :: d=0.1,u=2.3,v=4.5,w=6.7,sv=8.9
    real*8  :: dx=1.6
    integer  :: case_lb,case_ub
    real*8,allocatable :: fp(:,:),fn(:,:), fpx(:,:),fnx(:,:)
    integer :: i, n0=1, nn=64

    write(*,*)"verify function" 
    allocate(fp(1-6:nn+6,5),fn(1-6:nn+6,5),fpx(nn,5),fnx(nn,5))
    do i=1-6,nn+6
        call steger_warming(k1,k2,k3,d,u,v,w,sv,fp(i,:),fn(i,:))
    enddo
    do case_lb=1,2
        do case_ub=1,2
            call dfdx_p(case_lb,case_ub,n0,nn,dx,fp,fpx)
            call dfdx_n(case_lb,case_ub,n0,nn,dx,fn,fnx)
        enddo
    enddo

    deallocate(fp,fn,fpx,fnx)
    write(*,*)"verify over" 
end program verify_code