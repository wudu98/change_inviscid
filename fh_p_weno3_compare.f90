program fh_p_weno3_compare
    implicit none
    
    integer ::  i
    real*8  ::  f_1=3.5,f0=4.3,f1=5.0  !weno7�Ĳ���,���⸳ֵ
    real*8  ::  fhp,c_fhp,fhn,c_fhn    
    real ::start, finish ,ft,ct    !��ʱ����
    
    write(*,*)"fh_p_weno3  function" 
    ! fortran test
    call cpu_time(start)
    do i=0,1000000
    call fh_p_weno3(f_1,f0,f1,fhp)
    end do
    call cpu_time(finish)
    ft=finish-start
    write(*,*)"fortran running time(s) =",ft
    write(*,*)"f fhp =",fhp   
    
    ! c test
    call cpu_time(start)
    do i=0,1000000
    call fh_p_weno3_c(f_1,f0,f1,c_fhp)
    end do
    call cpu_time(finish)
    ct=finish-start
    write(*,*)"c running time(s) =",ct 
    write(*,*)"c fhp =",c_fhp 
    
    ! running time comparision 
    write(*,*)"(ctime-ftime)/ftime (%) =",(ct-ft)/ft *100
    write(*,*)"-----------------------------------"     
    
    write(*,*)"fh_n_weno3  function" 
    !fortran test
    call cpu_time(start)
    do i=0,1000000
    call fh_n_weno3(f1,f0,f_1,fhn)
    end do
    call cpu_time(finish)
    ft=finish-start
    write(*,*)"fortran running time(s) =",ft
    write(*,*)"f fhn =",fhn   

    !c test
    call cpu_time(start)
    do i=0,1000000
    call fh_p_weno3_c(f_1,f0,f1,c_fhn)
    end do
    call cpu_time(finish)
    ct=finish-start
    write(*,*)"c running time(s) =",ct 
    write(*,*)"c fhp =",c_fhn 
    
    ! running time comparision 
    write(*,*)"(ctime-ftime)/ftime (%) =",(ct-ft)/ft *100
    write(*,*)"-----------------------------------"     
    
end program fh_p_weno3_compare


