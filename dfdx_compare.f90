program dfdx_compare
    implicit none

    integer :: case_lb=3,case_ub=3
	integer :: jj=10,nn=17
    real*8  :: dx=1.6,fp(1-6:64+6,5)=0.0
    real*8  :: fpx(1:64,5)=0.0,c_fpx(1:64,5)=0.0
    
    integer :: i,j
    real*8  :: x
    real ::start, finish ,ft, ct    !��ʱ����
    
    do i=1-6,64+6
        do j=1,5
            call random_number(x)
            fp(i,j) = x
        end do
    end do

    write(*,*)"dfdx_p function" 
    ! fortran test
    ! call cpu_time(start)
    ! do i=0,100000
    call dfdx_p(case_lb,case_ub,jj,nn,dx,fp,fpx)
    ! end do
    ! call cpu_time(finish)
    ! ft=finish-start
    ! write(*,*)"fortran running time(s) =",ft  
    
    ! c test
    ! call cpu_time(start)
    ! do i=0,100000
    call dfdx_p_c(case_lb,case_ub,jj,nn,dx,fp,c_fpx)
    ! end do
    ! call cpu_time(finish)
    ! ct=finish-start
    ! write(*,*)"c running time(s) =",ct 
    
    do i=jj,nn
        do j=1,5
            write(*,*)"line = ",i,j,c_fpx(i,j),fpx(i,j)
        end do
    end do

    ! running time comparision 
    ! write(*,*)"(ctime-ftime)/ftime (%) =",(ct-ft)/ft * 100
    write(*,*)"-----------------------------------"     

    
end program dfdx_compare


