program fh_p_weno7_compare
    implicit none
    
    integer ::  i
    real*8  ::  f_3=1.9,f_2=2.2,f_1=3.5,f0=4.3,f1=5.0,f2=6.8,f3=7.1  !weno7�Ĳ���,���⸳ֵ
    real*8  ::  fhp,c_fhp,fhn,c_fhn    
    real ::start, finish ,ft,ct    !��ʱ����
    
    write(*,*)"fh_p_weno7  function" 
    ! fortran test
    call cpu_time(start)
    do i=0,1000000
    call fh_p_weno7(f_3,f_2,f_1,f0,f1,f2,f3,fhp)
    end do
    call cpu_time(finish)
    ft=finish-start
    write(*,*)"fortran running time(s) =",ft
    write(*,*)"f fhp =",fhp   
    
    ! c test
    call cpu_time(start)
    do i=0,1000000
    call fh_p_weno7_c(f_3,f_2,f_1,f0,f1,f2,f3,c_fhp)
    end do
    call cpu_time(finish)
    ct=finish-start
    write(*,*)"c running time(s) =",ct 
    write(*,*)"c fhp =",c_fhp 
    
    ! running time comparision 
    write(*,*)"(ctime-ftime)/ftime (%) =",(ct-ft)/ft *100
    write(*,*)"-----------------------------------"     
    
    write(*,*)"fh_n_weno7  function" 
    !fortran test
    call cpu_time(start)
    do i=0,1000000
    call fh_n_weno7(f3,f2,f1,f0,f_1,f_2,f_3,fhn)
    end do
    call cpu_time(finish)
    ft=finish-start
    write(*,*)"fortran running time(s) =",ft
    write(*,*)"f fhn =",fhn   

    !c test
    call cpu_time(start)
    do i=0,1000000
    call fh_p_weno7_c(f_3,f_2,f_1,f0,f1,f2,f3,c_fhn)
    end do
    call cpu_time(finish)
    ct=finish-start
    write(*,*)"c running time(s) =",ct 
    write(*,*)"c fhp =",c_fhn 
    
    ! running time comparision 
    write(*,*)"(ctime-ftime)/ftime (%) =",(ct-ft)/ft *100
    write(*,*)"-----------------------------------"     
    
end program fh_p_weno7_compare


