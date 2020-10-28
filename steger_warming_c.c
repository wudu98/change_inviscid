#include <math.h>
#include <stdio.h>
#include <math.h>

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



