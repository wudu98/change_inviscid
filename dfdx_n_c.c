#include <math.h>
#include <stdio.h>

extern void fh_p_weno7_c_(const double *f_3, const double *f_2, const double *f_1, const double *f0, const double *f1, const double *f2, const double *f3, double* fhp);
extern void fh_p_weno5_c_(const double *f_2, const double *f_1, const double *f0, const double *f1, const double *f2, double* fhp);
extern void fh_p_weno3_c_(const double *f_1, const double *f0, const double *f1, double* fhp);


void dfdx_n_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double* dx, double fn[], double fnx[]){
	int i,m;
	const int BE = 6;
	const int nn0 = *nn;
	const int nn12 = *nn+2*BE;  //数组变化
	double fhn[5][nn12+1];
	for (int m = 1-1; m <= 5-1; m++) {
		for (int i = (*jj)+3+BE-1; i <= (*nn)-3+BE-1; i++)
			fh_p_weno7_c_(&fn[m*nn12+i+3],&fn[m*nn12+i+2],&fn[m*nn12+i+1],&fn[m*nn12+i],&fn[m*nn12+i-1],&fn[m*nn12+i-2],&fn[m*nn12+i-3],&fhn[m][i]);
		for (int i = (*jj)+3+BE-1; i <= (*nn)-4+BE-1; i++)
			fnx[m*nn0+i-BE] = (fhn[m][i+1] - fhn[m][i]) / *dx;
		if(*case_lb == 1){
			i=(*jj)+2+BE-1;
			fh_p_weno5_c_(&fn[m*nn12+i+2],&fn[m*nn12+i+1],&fn[m*nn12+i],&fn[m*nn12+i-1],&fn[m*nn12+i-2],&fhn[m][i]);
			fnx[m*nn0+i-BE] = (fhn[m][i+1] - fhn[m][i]) / *dx;
			i=(*jj)+1+BE-1;
			fh_p_weno3_c_(&fn[m*nn12+i+1],&fn[m*nn12+i],&fn[m*nn12+i-1],&fhn[m][i]);
			fnx[m*nn0+i-BE] = (fhn[m][i+1]  - fhn[m][i]  )     / *dx;
			i=(*jj)+BE-1;
		  	fnx[m*nn0+i  -BE] = (fn[m*nn12+i+1] - fn[m*nn12+i] ) / *dx;
		}
		else{
			for(int i=(*jj)+BE-1; i<=(*jj)+2+BE-1; i++){
				fh_p_weno7_c_(&fn[m*nn12+i+3],&fn[m*nn12+i+2],&fn[m*nn12+i+1],&fn[m*nn12+i],&fn[m*nn12+i-1],&fn[m*nn12+i-2],&fn[m*nn12+i-3],&fhn[m][i]);
			}
			for(int i=(*jj)+BE-1; i<=(*jj)+2+BE-1; i++){
				fnx[m*nn0+i-BE] = (fhn[m][i+1] - fhn[m][i]) / *dx;
			}
		}
		if(*case_ub == 1){
			i=(*nn)-2+BE-1;
			fh_p_weno5_c_(&fn[m*nn12+i+2],&fn[m*nn12+i+1],&fn[m*nn12+i],&fn[m*nn12+i-1],&fn[m*nn12+i-2],&fhn[m][i]);
			fnx[m*nn0+i-1-BE] = (fhn[m][i] - fhn[m][i-1]) / *dx;
			i=(*nn)-1+BE-1;
			fh_p_weno3_c_(&fn[m*nn12+i+1],&fn[m*nn12+i],&fn[m*nn12+i-1],&fhn[m][i]);
			fnx[m*nn0+i-1-BE] = (fhn[m][i] - fhn[m][i-1]) / *dx;
			i=(*nn)+BE-1;
			fnx[m*nn0+i-1-BE] = (fn[m*nn12+nn0+BE-1] - fn[m*nn12+nn0-1+BE-1] ) / *dx;
			fnx[m*nn0+i  -BE] = (fn[m*nn12+nn0+BE-1] - fn[m*nn12+nn0-1+BE-1] ) / *dx;
		}
		else{
			for(int i=(*nn)-2+BE-1; i<=(*nn)+1+BE-1; i++){
				fh_p_weno7_c_(&fn[m*nn12+i+3],&fn[m*nn12+i+2],&fn[m*nn12+i+1],&fn[m*nn12+i],&fn[m*nn12+i-1],&fn[m*nn12+i-2],&fn[m*nn12+i-3],&fhn[m][i]);
				fnx[m*nn0+i-1-BE] = (fhn[m][i] - fhn[m][i-1]) / *dx;
			}
		}
		
		
	}

}



