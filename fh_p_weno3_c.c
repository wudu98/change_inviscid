#include <math.h>
#include <stdio.h>
void fh_p_weno3_c_(const double *f_1, const double *f0, const double *f1, double* fhp){ //positively biased
	int i;
	double w[2], beta[2], alfa[2], sumalfa;
	double fh[2];
	double d[2]={2.0/3.0, 1.0/3.0};
	double epsilon = 1e-6; 

	fh[0] =  0.5 * *f0  + 0.5 * *f1 ;
	fh[1] = -0.5 * *f_1 + 1.5 * *f0 ;

	beta[0] = (*f1 -  *f0 ) * (*f1 -  *f0 );
	beta[1] = (*f0 -  *f_1) * (*f0 -  *f_1);
	
	for (int i = 0; i <= 1; i++) {
		alfa[i]=d[i]/((epsilon+beta[i])*(epsilon+beta[i]));
	}
	sumalfa=alfa[0]+alfa[1];
	for (int i = 0; i <= 1; i++) {
		w[i]=alfa[i]/sumalfa;
	}

	*fhp = 0.0;
	for (int i = 0; i <= 1; i++) {
		*fhp = *fhp + w[i] * fh[i];
	}
	
}



