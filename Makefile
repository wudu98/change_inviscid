all : 
	gfortran -c inviscid.f90

	gcc -c fh_p_weno7_c.c
	# gfortran -c fh_p_weno7_compare.f90	
	# gfortran -pg -o fh_p_weno7_compare.out inviscid.o fh_p_weno7_c.o fh_p_weno7_compare.o
    
	gcc -c fh_p_weno5_c.c
	# gfortran -c fh_p_weno5_compare.f90	
	# gfortran -pg -o fh_p_weno5_compare.out inviscid.o fh_p_weno5_c.o fh_p_weno5_compare.o

	gcc -c fh_p_weno3_c.c 
	# gfortran -c fh_p_weno3_compare.f90	
	# gfortran -pg -o fh_p_weno3_compare.out inviscid.o fh_p_weno3_c.o fh_p_weno3_compare.o

	gcc -c dfdx_p_c.c 
	gcc -c dfdx_n_c.c
	gfortran -c dfdx_compare.f90 
	gfortran -pg -o dfdx_compare.out inviscid.o dfdx_p_c.o dfdx_n_c.o dfdx_compare.o fh_p_weno7_c.o fh_p_weno5_c.o fh_p_weno3_c.o

	gcc -c steger_warming_c.c 
	gfortran -c warm_compare.f90 
	gfortran -pg -o warm_compare.out inviscid.o steger_warming_c.o warm_compare.o
clean:
	rm -f *.o *.out