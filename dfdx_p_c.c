#include <math.h>
#include <stdio.h>

extern void fh_p_weno7_c_(const double *f_3, const double *f_2, const double *f_1, const double *f0, const double *f1, const double *f2, const double *f3, double* fhp);
extern void fh_p_weno5_c_(const double *f_2, const double *f_1, const double *f0, const double *f1, const double *f2, double* fhp);
extern void fh_p_weno3_c_(const double *f_1, const double *f0, const double *f1, double* fhp);


void dfdx_p_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double* dx, double fp[], double fpx[]){
	int i,m;
	const int BE = 6;
	const int nn0 = *nn;
	const int nn12 = *nn+2*BE;  //数组变化
	double fhp[5][nn12+1];
	for (int m = 1-1; m <= 5-1; m++) {
		for (int i = (*jj)+3+BE-1; i <= (*nn)-3+BE-1; i++)
			fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[m][i]);
		for (int i = (*jj)+4+BE-1; i <= (*nn)-3+BE-1; i++)
			fpx[m*nn0+i-BE] = (fhp[m][i] - fhp[m][i-1]) / *dx;
		if(*case_lb == 1){
			i=(*jj)+2+BE-1;
			fh_p_weno5_c_(&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fhp[m][i]);
			fpx[m*nn0+i+1-BE] = (fhp[m][i+1] - fhp[m][i]) / *dx;
			i=(*jj)+1+BE-1;
			fh_p_weno3_c_(&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fhp[m][i]);
			fpx[m*nn0+i+1-BE] = (fhp[m][i+1]  - fhp[m][i]  )     / *dx;
		  	fpx[m*nn0+i  -BE] = (fp[m*nn12+i] - fp[m*nn12+i-1] ) / *dx;
			fpx[m*nn0+i-1-BE] = (fp[m*nn12+i] - fp[m*nn12+i-1] ) / *dx;
		}
		else{
			for(int i=(*jj)+2+BE-1; i>=(*jj)-1+BE-1; i--){
				fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[m][i]);
				fpx[m*nn0+i+1-BE] = (fhp[m][i+1] - fhp[m][i]) / *dx;
			}
		}
		if(*case_ub == 1){
			fpx[m*nn0+nn0+BE-1]=(fp[m*nn12+nn0+BE-1] - fp[m*nn12+nn0-1+BE-1] ) / *dx;
			i=(*nn)-2+BE-1;
			fh_p_weno5_c_(&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fhp[m][i]);
			fpx[m*nn0+i-BE] = (fhp[m][i] - fhp[m][i-1]) / *dx;
			i=(*nn)-1+BE-1;
			fh_p_weno3_c_(&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fhp[m][i]);
			fpx[m*nn0+i-BE] = (fhp[m][i] - fhp[m][i-1]) / *dx;
			i=(*nn)+BE-1;
			fpx[m*nn0+i-BE] = (fp[m*nn12+nn0+BE-1] - fp[m*nn12+nn0-1+BE-1] ) / *dx;
		}
		else{
			for(int i=(*nn)-2+BE-1; i<=(*nn)+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[m][i]);
				fpx[m*nn0+i-BE] = (fhp[m][i] - fhp[m][i-1]) / *dx;
			}
		}
		
		
	}

}



