#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void steger_warming_c_(const double* k1, const double* k2, const double* k3, const double* d, 
                       const double* u , const double* v , const double* w , const double* sv, double* fp, double* fn){
    double k1b,k2b,k3b,kk;
    double lambda[3],lamb[3];
    int i;
    const double gama = 1.4;
    const double gama2_2 = 0.8;
    const double gama2 = 2.8;
    const double gama_1 = 0.4;
    const double gama_w = 2.0;
    const double epsilon_2 = 1e-12;

    kk = sqrt((*k1)*(*k1)+(*k2)*(*k2)+(*k3)*(*k3));
    lambda[0]=(*k1)*(*u)+(*k2)*(*v)+(*k3)*(*w);
	lambda[1]=lambda[0]+(*sv)*kk;
	lambda[2]=lambda[0]-(*sv)*kk;
    k1b=(*k1)/kk;
	k2b=(*k2)/kk;
	k3b=(*k3)/kk;
    double usv  = (*u)+(*sv)*k1b;
    double vsv  = (*v)+(*sv)*k2b;
    double wsv  = (*w)+(*sv)*k3b;
    double usv_ = (*u)-(*sv)*k1b;
    double vsv_ = (*v)-(*sv)*k2b;
    double wsv_ = (*w)-(*sv)*k3b;
    for(int i=0; i<3; i++)
	    lamb[i]=0.5*(lambda[i]+sqrt(lambda[i]*lambda[i]+epsilon_2));
    fp[0]=(*d) / gama2 * (gama2_2 * lamb[0]      + lamb[1]     + lamb[2]     );
	fp[1]=(*d) / gama2 * (gama2_2 * lamb[0]*(*u) + lamb[1]*usv + lamb[2]*usv_);
	fp[2]=(*d) / gama2 * (gama2_2 * lamb[0]*(*v) + lamb[1]*vsv + lamb[2]*vsv_);
	fp[3]=(*d) / gama2 * (gama2_2 * lamb[0]*(*w) + lamb[1]*wsv + lamb[2]*wsv_);
	fp[4]=(*d) / gama2 * (gama_1  * lamb[0]*((*u)*(*u) + (*v)*(*v) + (*w)*(*w))+
		                    0.5   * lamb[1]*( usv*usv  +  vsv*vsv  +  wsv*wsv )+
		                    0.5   * lamb[2]*(usv_*usv_ + vsv_*vsv_ + wsv_*wsv_)+
				          gama_w  * (lamb[1] + lamb[2]) * (*sv)*(*sv)          );
    for(int i=0; i<3; i++)
	    lamb[i]=0.5*(lambda[i]-sqrt(lambda[i]*lambda[i]+epsilon_2));
    fn[0]=(*d) / gama2 * (gama2_2 * lamb[0]      + lamb[1]     + lamb[2]     );
	fn[1]=(*d) / gama2 * (gama2_2 * lamb[0]*(*u) + lamb[1]*usv + lamb[2]*usv_);
	fn[2]=(*d) / gama2 * (gama2_2 * lamb[0]*(*v) + lamb[1]*vsv + lamb[2]*vsv_);
	fn[3]=(*d) / gama2 * (gama2_2 * lamb[0]*(*w) + lamb[1]*wsv + lamb[2]*wsv_);
	fn[4]=(*d) / gama2 * (gama_1  * lamb[0]*((*u)*(*u) + (*v)*(*v) + (*w)*(*w))+
		                    0.5   * lamb[1]*( usv*usv  +  vsv*vsv  +  wsv*wsv )+
		                    0.5   * lamb[2]*(usv_*usv_ + vsv_*vsv_ + wsv_*wsv_)+
				          gama_w  * (lamb[1] + lamb[2]) * (*sv)*(*sv)          );
}



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


void fh_p_weno7_c_(const double *f_3, const double *f_2, const double *f_1, const double *f0, const double *f1, const double *f2, const double *f3, double* fhp){ //positively biased
	int i;
	double w[4], beta[4], alfa[4], sumalfa;
	double fh[4];
	double d[4];
	double TV[4], TVR, TV_MAX, TV_MIN;
	double epsilon = 1e-6; 

	//optimal weight
	d[0] = 1.0 / 35.0;//fortran: d(0)=1.0d0/35.0d0  �������⣿��
	d[1] = 12.0 / 35.0;
	d[2] = 18.0 / 35.0;
	d[3] = 4.0 / 35.0;

	fh[0] = -1.0 / 4.0  * *f_3 + 13.0 / 12.0 * *f_2 - 23.0 / 12.0 * *f_1 + 25.0 / 12.0 * *f0;
	fh[1] =  1.0 / 12.0 * *f_2 -  5.0 / 12.0 * *f_1 + 13.0 / 12.0 * *f0  +  1.0 / 4.0  * *f1;
	fh[2] = -1.0 / 12.0 * *f_1 +  7.0 / 12.0 * *f0  +  7.0 / 12.0 * *f1  -  1.0 / 12.0 * *f2;
	fh[3] =  1.0 / 4.0  * *f0  + 13.0 / 12.0 * *f1  -  5.0 / 12.0 * *f2  +  1.0 / 12.0 * *f3;

	TV[0] = fabs(*f_3 - *f_2) + fabs(*f_2 - *f_1) + fabs(*f_1 - *f0); // TV(0)=abs(*f_3-*f_2)+abs(*f_2-*f_1)+abs(*f_1-*f0)
	TV[1] = fabs(*f_2 - *f_1) + fabs(*f_1 - *f0)  + fabs(*f0 - *f1);
	TV[2] = fabs(*f_1 - *f0)  + fabs(*f0 - *f1)   + fabs(*f1 - *f2);
	TV[3] = fabs(*f0 - *f1)   + fabs(*f1 - *f2)   + fabs(*f2 - *f3);

	//TV_MAX = MaxVal(TV)
	//TV_MIN = MinVal(TV)
	TV_MAX = TV[0];
	TV_MIN = TV[0];
	for (int i = 1; i < 4; i++) {    //4 is the size of TV
		TV_MAX = TV_MAX >= TV[i] ? TV_MAX : TV[i];
		TV_MIN = TV_MIN <= TV[i] ? TV_MIN : TV[i];
	}
	
	
	TVR = TV_MAX / (TV_MIN + epsilon); 
	//for test
	//printf("%lf,%lf,%lf,%lf,", TV[0], TV[1], TV[2], TV[3]);
	//printf("c TV_MAX=%lf,TV_MIN=%lf,TVR=%lf", TV_MAX, TV_MIN, TVR);

	if (TV_MAX < 0.2 && TVR < 5.0) {
		for (int i = 0; i <= 3; i++) {
			w[i] = d[i];
		}
	}

	else {
		beta[0] = *f_3 * (547.0   * *f_3 - 3882.0  * *f_2 + 4642.0 * *f_1 - 1854.0 * *f0) 
				+ *f_2 * (7043.0  * *f_2 - 17246.0 * *f_1 + 7042.0 * *f0)
				+ *f_1 * (11003.0 * *f_1 - 9402.0  * *f0) 
				+ 2107.0 * (*f0**f0);
		beta[1] = *f_2 * (267.0   * *f_2 - 1642.0  * *f_1 + 1602.0 * *f0  - 494.0  * *f1) 
				+ *f_1 * (2843.0  * *f_1 - 5966.0  * *f0  + 1922.0 * *f1)
				+ *f0  * (3443.0  * *f0  - 2522.0  * *f1) 
				+ 547.0 * (*f1**f1);
		beta[2] = *f_1 * (547.0   * *f_1 - 2522.0  * *f0  + 1922.0 * *f1  - 494.0  * *f2) 
				+ *f0  * (3443.0  * *f0  - 5966.0  * *f1  + 1602.0 * *f2) 
				+ *f1  * (2843.0  * *f1  - 1642.0  * *f2) 
				+ 267.0 * (*f2**f2) ;
		beta[3] = *f0  * (2107.0  * *f0  - 9402.0  * *f1  + 7042.0 * *f2  - 1854.0 * *f3) 
				+ *f1  * (11003.0 * *f1  - 17246.0 * *f2  + 4642.0 * *f3) 
				+ *f2  * (7043.0  * *f2  - 3882.0  * *f3) 
				+ 547.0 * (*f3**f3) ;

		for (int i = 0; i <= 3; i++) {
			alfa[i] = d[i] / ((epsilon + beta[i]) * (epsilon + beta[i])); //alfa(i)=d(i)/(epsilon+beta(i))**2
		}

		sumalfa = alfa[0] + alfa[1] + alfa[2] + alfa[3]; 

		//for test
		//printf("sumalfa=%.16f,alfa[0]=%.16f", sumalfa, alfa[0]);

		for (int i = 0; i <= 3; i++) {
			w[i] = alfa[i] / sumalfa;
		}
	}
	*fhp = 0.0;
	for (int i = 0; i <= 3; i++) {
		*fhp = *fhp + w[i] * fh[i];
	}
	
}


void dfdx_p_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double *dx, double fp[], double fpx[]){
	int i,m;
	const int BE_1 = 5;
	const int BE = 6;
	const int nn0 = *nn;
	const int nn12 = *nn+2*BE;  
	double fhp[2];
	int flag = 0;
	if((*case_lb) == 1 && (*case_ub) == 1){
		for (int m = 0; m <= 4; m++) {
			i=(*jj)+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+i+1] - fp[m*nn12+i]) ;
			i=(*jj)+1+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+i  ]   - fp[m*nn12+i-1]) ;
				fh_p_weno3_c_(&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fhp[flag=flag^1]);
			i=(*jj)+2+BE-1;
				fh_p_weno5_c_(&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  )     ;
			for (int i = (*jj)+3+BE-1; i <= (*nn)-3+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			}
			i=(*nn)-2+BE-1;
				fh_p_weno5_c_(&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			i=(*nn)-1+BE-1;
				fh_p_weno3_c_(&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			i=(*nn)+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+nn0+BE-1] - fp[m*nn12+nn0-1+BE-1] ) ;
		}
		
	}
	else if((*case_lb) == 1){
		for (int m = 0; m <= 4; m++) {
			i=(*jj)+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+i+1] - fp[m*nn12+i]) ;
			i=(*jj)+1+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+i  ]   - fp[m*nn12+i-1]) ;
				fh_p_weno3_c_(&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fhp[flag=flag^1]);
			i=(*jj)+2+BE-1;
				fh_p_weno5_c_(&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  )     ;
			for (int i = (*jj)+3+BE-1; i <= (*nn)+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			}
		}
	}
	else if((*case_ub) == 1){
		for (int m = 0; m <= 4; m++) {
			i = (*jj)-1+BE-1;
				fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[flag=flag^1]);
			for (int i = (*jj)+BE-1; i <= (*nn)-3+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			}
			i=(*nn)-2+BE-1;
				fh_p_weno5_c_(&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			i=(*nn)-1-1;
				fh_p_weno3_c_(&fp[m*nn12+i-1+BE],&fp[m*nn12+i+BE],&fp[m*nn12+i+1+BE],&fhp[flag=flag^1]);
				fpx[m*nn0+i] = (fhp[flag]  - fhp[flag^1]  ) ;
			i=(*nn)-1;
				fpx[m*nn0+i] = (fp[m*nn12+nn0+BE-1] - fp[m*nn12+nn0-1+BE-1] ) ;
		}
	}
	else{
		for (int m = 0; m <= 4; m++) {
			i = (*jj)-1+BE-1;
				fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[flag=flag^1]);
			for (int i = (*jj)+BE-1; i <= (*nn)+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			}
		}
	}
	double _dx = 1.0 / *dx;
	for (int m = 0; m <= 4; m++) {
		for (int i = (*jj)-1; i <= (*nn)-1; i++)
			fpx[m*nn0+i] = fpx[m*nn0+i] * _dx;
	}
}


void dfdx_p_v0_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double* dx, double fp[], double fpx[]){
	int i,m;
	const int BE = 6;
	const int nn0 = *nn;
	const int nn12 = *nn+2*BE; 
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

/*  same as dfdx_p_c_  just adjust "fh_p_weno7_c_" input*/
void dfdx_n_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double *dx, double fp[], double fpx[]){
	int i,m;
	const int BE_1 = 5;
	const int BE = 6;
	const int nn0 = *nn;
	const int nn12 = *nn+2*BE;  
	double fhp[2];
	int flag = 0;
	if((*case_lb) == 1 && (*case_ub) == 1){
		for (int m = 0; m <= 4; m++) {
			i=(*jj)+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+i+1] - fp[m*nn12+i]) ;
			i=(*jj)+1+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+i  ]   - fp[m*nn12+i-1]) ;
				fh_p_weno3_c_(&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fhp[flag=flag^1]);
			i=(*jj)+2+BE-1;
				fh_p_weno5_c_(&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  )     ;
			for (int i = (*jj)+3+BE-1; i <= (*nn)-3+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i+3],&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fp[m*nn12+i-3],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			}
			i=(*nn)-2+BE-1;
				fh_p_weno5_c_(&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			i=(*nn)-1+BE-1;
				fh_p_weno3_c_(&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			i=(*nn)+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+nn0+BE-1] - fp[m*nn12+nn0-1+BE-1] ) ;
		}
		
	}
	else if((*case_lb) == 1){
		for (int m = 0; m <= 4; m++) {
			i=(*jj)+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+i+1] - fp[m*nn12+i]) ;
			i=(*jj)+1+BE-1;
				fpx[m*nn0+i-BE] = (fp[m*nn12+i  ]   - fp[m*nn12+i-1]) ;
				fh_p_weno3_c_(&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fhp[flag=flag^1]);
			i=(*jj)+2+BE-1;
				fh_p_weno5_c_(&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  )     ;
			for (int i = (*jj)+3+BE-1; i <= (*nn)+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i+3],&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fp[m*nn12+i-3],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			}
		}
	}
	else if((*case_ub) == 1){
		for (int m = 0; m <= 4; m++) {
			i = (*jj)-1+BE-1;
				fh_p_weno7_c_(&fp[m*nn12+i+3],&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fp[m*nn12+i-3],&fhp[flag=flag^1]);
			for (int i = (*jj)+BE-1; i <= (*nn)-3+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i+3],&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fp[m*nn12+i-3],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			}
			i=(*nn)-2+BE-1;
				fh_p_weno5_c_(&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			i=(*nn)-1-1;
				fh_p_weno3_c_(&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fhp[flag=flag^1]);
				fpx[m*nn0+i] = (fhp[flag]  - fhp[flag^1]  ) ;
			i=(*nn)-1;
				fpx[m*nn0+i] = (fp[m*nn12+nn0+BE-1] - fp[m*nn12+nn0-1+BE-1] ) ;
		}
	}
	else{
		for (int m = 0; m <= 4; m++) {
			i = (*jj)-1+BE-1;
				fh_p_weno7_c_(&fp[m*nn12+i+3],&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fp[m*nn12+i-3],&fhp[flag=flag^1]);
			for (int i = (*jj)+BE-1; i <= (*nn)+BE-1; i++){
				fh_p_weno7_c_(&fp[m*nn12+i+3],&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fp[m*nn12+i-3],&fhp[flag=flag^1]);
				fpx[m*nn0+i-BE] = (fhp[flag]  - fhp[flag^1]  ) ;
			}
		}
	}
	double _dx = 1.0 / *dx;
	for (int m = 0; m <= 4; m++) {
		for (int i = (*jj)-1; i <= (*nn)-1; i++)
			fpx[m*nn0+i] = fpx[m*nn0+i] * _dx;
	}
}

void dfdx_n_v0_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double* dx, double fn[], double fnx[]){
	int i,m;
	const int BE = 6;
	const int nn0 = *nn;
	const int nn12 = *nn+2*BE;  
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









