program warm_compare
    implicit none 

	real*8  :: k1=1.1,k2=2.2,k3=3.3
    real*8  :: d=0.1,u=2.3,v=4.5,w=6.7,sv=8.9
    real*8,allocatable :: fp(:,:),fn(:,:), c_fp(:,:),c_fn(:,:)
    integer :: i, j, nn=64
    real :: loss, start, finish ,ft, ct

    allocate(fp(1-6:nn+6,5),fn(1-6:nn+6,5),c_fp(1-6:nn+6,5),c_fn(1-6:nn+6,5))

    write(*,*)"dfdx_p function" 
    ! fortran test
    call cpu_time(start)
    do i=0,10000
        do j=1-6,nn+6
            call steger_warming(k1,k2,k3,d,u,v,w,sv,fp(j,:),fn(j,:))
        enddo
    end do
    call cpu_time(finish)
    ft=finish-start
    write(*,*)"fortran running time(s) =",ft  

    ! c test
    call cpu_time(start)
    do i=0,10000
        do j=1-6,nn+6
            call steger_warming_c(k1,k2,k3,d,u,v,w,sv,c_fp(j,:),c_fn(j,:))
        enddo
    end do
    call cpu_time(finish)
    ct=finish-start
    write(*,*)"c running time(s) =",ct 

    ! running time comparision 
    write(*,*)"(ctime-ftime)/ftime (%) =",(ct-ft)/ft * 100

    write(*,*)"-----------------------------------" 

    loss = 0.0
    do i=1-6,nn+6
        do j=1,5
            loss = loss + abs(fp(i,j)-c_fp(i,j))
            loss = loss + abs(fn(i,j)-c_fn(i,j))
        end do
    end do
    write(*,*)"loss=",loss
    deallocate(fp,fn,c_fp,c_fn)

end program warm_compare


