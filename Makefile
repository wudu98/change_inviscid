all : 
	gcc -c fh_p_weno7_c.c
	gfortran -c fh_p_weno7_f.f90
	gfortran -c h_p_weno7_compare.f90
	gfortran -pg -o fh_p_weno7_compare.out fh_p_weno7_c.o fh_p_weno7_f.o fh_p_weno7_compare.o
    
clean:
	rm -f *.o *.out