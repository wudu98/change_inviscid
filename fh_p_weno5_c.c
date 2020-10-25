#include <math.h>
#include <stdio.h>
void fh_p_weno5_c_(const double *f_2, const double *f_1, const double *f0, const double *f1, const double *f2, double* fhp){ //positively biased
	int i;
	double w[3], beta[3], alfa[3], sumalfa;
	double fh[3];
	double d[3]={0.3,0.6,0.1};
	double epsilon = 1e-6; 

	fh[0] =  1.0 / 3.0 * *f0  + 5.0 / 6.0 * *f1  -  1.0 / 6.0 * *f2 ;
	fh[1] = -1.0 / 6.0 * *f_1 + 5.0 / 6.0 * *f0  +  1.0 / 3.0 * *f1;
	fh[2] =  1.0 / 3.0 * *f_2 - 7.0 / 6.0 * *f_1 + 11.0 / 6.0 * *f0;

	beta[0] = 13.0/12.0 * (      *f0 - 2.0 * *f1  +       *f2) * (      *f0 - 2.0 * *f1  +       *f2)
			      +0.25 * (3.0 * *f0 - 4.0 * *f1  +       *f2) * (3.0 * *f0 - 4.0 * *f1  +       *f2);
	beta[1] = 13.0/12.0 * (     *f_1 - 2.0 * *f0  +       *f1) * (     *f_1 - 2.0 * *f0  +       *f1)
			      +0.25 * (     *f_1 -       *f1             ) * (     *f_1 -       *f1             );
	beta[2] = 13.0/12.0 * (     *f_2 - 2.0 * *f_1 +       *f0) * (     *f_2 - 2.0 * *f_1 +       *f0)
			      +0.25 * (     *f_2 - 4.0 * *f_1 + 3.0 * *f0) * (     *f_2 - 4.0 * *f_1 + 3.0 * *f0);
	
	for (int i = 0; i <= 2; i++) {
		alfa[i]=d[i]/((epsilon+beta[i])*(epsilon+beta[i]));
	}
	sumalfa=alfa[0]+alfa[1]+alfa[2];
	for (int i = 0; i <= 2; i++) {
		w[i]=alfa[i]/sumalfa;
	}

	*fhp = 0.0;
	for (int i = 0; i <= 2; i++) {
		*fhp = *fhp + w[i] * fh[i];
	}
	
}



