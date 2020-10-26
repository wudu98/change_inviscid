#include <math.h>
#include <stdio.h>

extern void fh_p_weno7_(double *f_3, double *f_2, double *f_1, double *f0, double *f1, double *f2, double *f3, double* fhp );
extern void fh_p_weno5_(double *f_2, double *f_1, double *f0, double *f1, double *f2, double* fhp );
extern void fh_p_weno3_(double *f_1, double *f0, double *f1, double* fhp );

void dfdx_p_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double* dx, double fp[][64], double fpx[][64]){
	int i,m;
	double fhp[5][*nn];
	for (int m = 0; m <= 4; m++) {
		for (int i = (*jj)+3; i <= (*nn)-3; i++)
			fh_p_weno7_(&fp[m][i-3],&fp[m][i-2],&fp[m][i-1],&fp[m][i],&fp[m][i+1],&fp[m][i+2],&fp[m][i+3],&fhp[m][i]);
		for (int i = (*jj)+4; i <= (*nn)-3; i++)
			fpx[m][i] = (fhp[m][i] - fhp[m][i-1]) / *dx;
		if(*case_lb == 1){
			i=(*jj)+2;
			fh_p_weno5_(&fp[m][i-2],&fp[m][i-1],&fp[m][i],&fp[m][i+1],&fp[m][i+2],&fhp[m][i]);
			fpx[m][i+1] = (fhp[m][i+1] - fhp[m][i]) / *dx;
			i=(*jj)+1;
			fh_p_weno3_(&fp[m][i-1],&fp[m][i],&fp[m][i+1],&fhp[m][i]);
			fpx[m][i+1]=(fhp[m][i+1] - fhp[m][i]  ) / *dx;
          	fpx[m][i-1]=(fp[m][i]    - fp[m][i-1] ) / *dx;
		  	fpx[m][i  ]=(fp[m][i]    - fp[m][i-1] ) / *dx;
		}
		else{
			for(int i=(*jj)+2; i>=(*jj)-1; i--){
				fh_p_weno7_(&fp[m][i-3],&fp[m][i-2],&fp[m][i-1],&fp[m][i],&fp[m][i+1],&fp[m][i+2],&fp[m][i+3],&fhp[m][i]);
				fpx[m][i+1] = (fhp[m][i+1] - fhp[m][i]) / *dx;
			}
		}
		if(*case_ub == 1){
			fpx[m][*nn]=(fp[m][*nn] - fp[m][(*nn)-1] ) / *dx;
			i=*nn-2;
			fh_p_weno5_(&fp[m][i-2],&fp[m][i-1],&fp[m][i],&fp[m][i+1],&fp[m][i+2],&fhp[m][i]);
			fpx[i][m] = (fhp[m][i] - fhp[m][i-1]) / *dx;
			i=(*nn)-1;
			fh_p_weno3_(&fp[m][i-1],&fp[m][i],&fp[m][i+1],&fhp[m][i]);
			fpx[m][i]=(fhp[m][i] - fhp[m][i-1]  ) / *dx;
		}
		else{
			for(int i=(*nn)-2; i<=*nn; i++){
				fh_p_weno7_(&fp[m][i-3],&fp[m][i-2],&fp[m][i-1],&fp[m][i],&fp[m][i+1],&fp[m][i+2],&fp[m][i+3],&fhp[m][i]);
				fpx[m][i] = (fhp[m][i] - fhp[m][i-1]) / *dx;
			}
		}
		
	}
	
	
}



