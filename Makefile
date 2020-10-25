all : 
	gcc -c fh_p_weno7_c.c
	gfortran -c fh_p_weno7_f.f90
	gfortran -c fh_n_weno7_f.f90
	gfortran -c fh_p_weno7_compare.f90	
	gfortran -pg -o fh_p_weno7_compare.out fh_p_weno7_c.o fh_p_weno7_f.o fh_n_weno7_f.o fh_p_weno7_compare.o
    
	gcc -c fh_p_weno5_c.c
	gfortran -c fh_p_weno5_f.f90
	gfortran -c fh_n_weno5_f.f90
	gfortran -c fh_p_weno5_compare.f90	
	gfortran -pg -o fh_p_weno5_compare.out fh_p_weno5_c.o fh_p_weno5_f.o fh_n_weno5_f.o fh_p_weno5_compare.o

	gcc -c fh_p_weno3_c.c
	gfortran -c fh_p_weno3_f.f90
	gfortran -c fh_n_weno3_f.f90
	gfortran -c fh_p_weno3_compare.f90	
	gfortran -pg -o fh_p_weno3_compare.out fh_p_weno3_c.o fh_p_weno3_f.o fh_n_weno3_f.o fh_p_weno3_compare.o

clean:
	rm -f *.o *.out