all : 
	gcc -c inviscid_c.c

	gfortran -c inviscid.f90
	gfortran -c inviscid_verify.f90

	# gfortran -c fh_p_weno7_compare.f90
	# gfortran -c fh_p_weno5_compare.f90
	# gfortran -c fh_p_weno3_compare.f90	
	gfortran -c dfdx_compare.f90
	# gfortran -c warm_compare.f90

	# gfortran -pg -o fh_p_weno7_compare.out inviscid.o inviscid_c.o fh_p_weno7_compare.o
	# gfortran -pg -o fh_p_weno5_compare.out inviscid.o inviscid_c.o fh_p_weno5_compare.o
	# gfortran -pg -o fh_p_weno3_compare.out inviscid.o inviscid_c.o fh_p_weno3_compare.o
	gfortran -pg -o dfdx_compare.out inviscid.o inviscid_c.o dfdx_compare.o
	# gfortran -pg -o warm_compare.out inviscid.o inviscid_c.o warm_compare.o

	# # verify
	gfortran -c verify_code.f90
	gfortran -pg -o verify_code.out inviscid_verify.o inviscid_c.o verify_code.o

clean:
	rm -f *.o *.out