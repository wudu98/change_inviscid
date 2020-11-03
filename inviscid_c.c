#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void steger_warming(const double k1, const double k2, const double k3, const double d,
	const double u, const double v, const double w, const double sv, double* fp, double* fn, const int step);
void dfdx_p(const int case_lb, const int case_ub, const int jj, const int nn, const double dx, double* fp, double* fpx);
void dfdx_n(const int case_lb, const int case_ub, const int jj, const int nn, const double dx, double* fn, double* fnx);
void fh_p_weno3_c_(const double* f_1, const double* f0, const double* f1, double* fhp);
void fh_p_weno5_c_(const double* f_2, const double* f_1, const double* f0, const double* f1, const double* f2, double* fhp);
void fh_p_weno7_c_(const double* f_3, const double* f_2, const double* f_1, const double* f0, const double* f1, const double* f2, const double* f3, double* fhp);


void inviscidc_(double* Lu1, double* Lu2, double* Lu3, double* Lu4, double* Lu5, int* na0, int* nr0,
	double* z, double* da, double* dc, double* dr, double* dadx, double* dcdy, double* drdz,
	double* JJ, double* Ua, double* Uc, double* Ur, double* d, double* p, double* sv) {  //!��ճ��ķ���ֵ

//Declaration�г����ض���
	int nx = 720, ny = 64, nz = 240, BE = 6;
	int dims[3] = { 16, 1, 4 };
	int nx1 = nx / 18;
	int nz1 = nz / 4;
	int nz2 = nz / 2.4;
	int na = 45, nc = 64, nr = 60;
	//�Զ�������
	double df1[na][nc][nr], df2[na][nc][nr], df3[na][nc][nr], df4[na][nc][nr], df5[na][nc][nr]; 
	double* fp, * fn, * fpx, * fnx;  
	//�Զ������
	int i, k, j, nn, n0;
	double k1, k2, k3, dx;
	double ddda, dpda, duada, ducda, durda;// ������߽����
	double dddr, dpdr, duadr, ducdr, durdr; //�������
	double phi1, phi2, phi3, phi4, phi5, d1, d2, d3, d4, d5;
	int case_lb, case_ub;
	int len_na, len_nc, len_nr;//used for malloc and calculate idx
	int idx, Lu_idx;

	len_na = na + BE - (1 - BE) + 1;
	len_nc = nc + BE - (1 - BE) + 1;
	len_nr = nr + BE - (1 - BE) + 1;

	//!������Ϊ��:��ȡ��һ��������ƽ�е�������, ȷ�������ұ߽��λ��n0��NN,
	//!���жϳ���߽�λ���Ƿ�Ϊ�����߽�, ����ǵĻ���ȡΪcase_LB����case_UB����1
	//!Ȼ�����dfdx.f90�Ƚ���ʸͨ���ķ���, �ٵ���WENO��ʽ�����ɢ����
	dx = *da;
	k2 = 0;
	k3 = 0;
	nn = na;
	////printf("41 inviscidc initial success\n");
	fp = (double*)malloc(5 * len_na * sizeof(double));
	fn = (double*)malloc(5 * len_na * sizeof(double));
	fpx = (double*)malloc(5 * na * sizeof(double));
	fnx = (double*)malloc(5 * na * sizeof(double));
	for (k = 1; k <= nr; k++) {
		for (j = 1; j <= nc; j++) {
			if (*nr0 + k < nz2) {
				if (nx1 > * na0 && nx1 <= *na0 + na) {
					for (i = 1; i <= na; i++) {
						if (*na0 + i == nx1) {
							n0 = i;
						}
					}
					case_lb = 1;
					if (*na0 + na == nx) {
						case_ub = 1;
					}
					else {
						case_ub = 0;
					}
				}
				else if (*na0 > nx1) {
					n0 = 1;
					case_lb = 0;
					if (*na0 + na == nx) {
						case_ub = 1;
					}
					else {
						case_ub = 0;  //!1�����������߽��ڴ�
					}
				}
			}
			else {//!(*nr0 + k >= nz2)then
				n0 = 1;
				if (*na0 == 0) {
					case_lb = 1;
					if (*na0 + na == nx) {
						case_ub = 1;
					}
					else {
						case_ub = 0;
					}//!1�����������߽��ڴ�
				}
				else {
					case_lb = 0;
					if (*na0 + na == nx) {
						case_ub = 1;
					}
					else {
						case_ub = 0;
					}//!1�����������߽��ڴ�
				}
			}
			{//include "dfdx.f90" 
				idx = len_na * len_nc * (k + 5) + len_na * (j + 5);//+ (x + 5);
				if (case_lb == 1 && case_ub == 1) {
					idx += n0 + 5;
					for (i = n0; i <= nn; i++) {
						//!k1 = dadx[i-1+BE] * JJ(i, j, k)
						k1 = 1.0 / drdz[k - 1 + BE];
						steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[i + BE - 1], &fn[i + BE - 1], len_na);
						idx++;
					}
					dfdx_p(case_lb, case_ub, n0, nn, dx, fp, fpx);
					dfdx_n(case_lb, case_ub, n0, nn, dx, fn, fnx);
				}
				if (case_lb == 1 && case_ub == 0) {
					idx += n0 + 5;
					for (i = n0; i <= nn + 6; i++) {
						//!k1 = dadx[i-1+BE] * JJ(i, j, k)
						k1 = 1.0 / drdz[k - 1 + BE]; 
						steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[i + BE - 1], &fn[i + BE - 1], len_na);
						idx++;
					}
					dfdx_p(case_lb, case_ub, n0, nn, dx, fp, fpx);
					dfdx_n(case_lb, case_ub, n0, nn, dx, fn, fnx);
				}
				if (case_lb == 0 && case_ub == 1) {
					idx += n0 - 6 + 5;
					for (i = n0 - 6; i <= nn; i++) {
						//!k1 = dadx[i-1+BE] * JJ(i, j, k)
						k1 = 1.0 / drdz[k - 1 + BE];
						steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[i + BE - 1], &fn[i + BE - 1], len_na);
						idx++;
					}
					dfdx_p(case_lb, case_ub, n0, nn, dx, fp, fpx);
					dfdx_n(case_lb, case_ub, n0, nn, dx, fn, fnx);
				}
				if (case_lb == 0 && case_ub == 0) {
					idx += n0 - 6 + 5;
					for (i = n0 - 6; i <= nn + 6; i++) {
						//!k1 = dadx[i-1+BE] * JJ(i, j, k)
						k1 = 1.0 / drdz[k - 1 + BE];
						steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[i + BE - 1], &fn[i + BE - 1], len_na);
						idx++;
					}
					dfdx_p(case_lb, case_ub, n0, nn, dx, fp, fpx);
					dfdx_n(case_lb, case_ub, n0, nn, dx, fn, fnx);
				}
				for (i = n0; i <= nn; i++) { 
					df1[i - 1][j - 1][k - 1] = fpx[i - 1] + fnx[i - 1];
					df2[i - 1][j - 1][k - 1] = fpx[i - 1 + nn] + fnx[i - 1 + nn];
					df3[i - 1][j - 1][k - 1] = fpx[i - 1 + nn * 2] + fnx[i - 1 + nn * 2];
					df4[i - 1][j - 1][k - 1] = fpx[i - 1 + nn * 3] + fnx[i - 1 + nn * 3];
					df5[i - 1][j - 1][k - 1] = fpx[i - 1 + nn * 4] + fnx[i - 1 + nn * 4];
				}
			}//include "dfdx.f90" 
		}//for
	}//for
	free(fp);
	free(fn);
	free(fpx);
	free(fnx);
	if (*na0 == 0) {
		i = 1;
		for (k = 1; k <= nr; k++) {
			if (*nr0 + k >= nz2) {
				for (j = 1; j <= nc; j++) {
					idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
					dpda = (-3.0 * p[idx] + 4.0 * p[idx + 1] - p[idx + 2]) / (2.0 * *da);
					ddda = (-3.0 * d[idx] + 4.0 * d[idx + 1] - d[idx + 2]) / (2.0 * *da);
					duada = (-3.0 * Ua[idx] + 4.0 * Ua[idx + 1] - Ua[idx + 2]) / (2.0 * *da);
					ducda = (-3.0 * Uc[idx] + 4.0 * Uc[idx + 1] - Uc[idx + 2]) / (2.0 * *da);
					durda = (-3.0 * Ur[idx] + 4.0 * Ur[idx + 1] - Ur[idx + 2]) / (2.0 * *da);
					if ((Ua[idx] - sv[idx]) < 0.0) {
						phi1 = (dadx[i - 1 + BE] * Ua[idx] - sv[idx] * dadx[i - 1 + BE]) * (-d[idx] * sv[idx] * duada + dpda);
					}
					else {
						phi1 = 0.0;
					}
					if (Ua[idx] < 0.0) {
						phi2 = dadx[i - 1 + BE] * Ua[idx] * (sv[idx] * sv[idx] * ddda - dpda); //**2
						phi3 = dadx[i - 1 + BE] * Ua[idx] * (dadx[i - 1 + BE] * ducda);
						phi4 = dadx[i - 1 + BE] * Ua[idx] * (dadx[i - 1 + BE] * durda);
					}
					else {
						phi2 = 0.0;
						phi3 = 0.0;
						phi4 = 0.0;
					}
					if ((Ua[idx] + sv[idx]) < 0.0) {
						phi5 = (dadx[i - 1 + BE] * Ua[idx] + sv[idx] * dadx[i - 1 + BE]) * (d[idx] * sv[idx] * duada + dpda);
					}
					else {
						phi5 = 0.0;
					}
					d1 = (0.5 * (phi1 + phi5) + phi2) / (sv[idx] * sv[idx]);//**2
					d2 = (phi5 - phi1) / (2.0 * d[idx] * sv[idx]);
					d3 = phi3 / dadx[i - 1 + BE];
					d4 = phi4 / dadx[i - 1 + BE];
					d5 = 0.5 * (phi1 + phi5);
					df1[i - 1][j - 1][k - 1] = (d1)*JJ[idx];
					df2[i - 1][j - 1][k - 1] = (Ua[idx] * d1 + d[idx] * d2) * JJ[idx];
					df3[i - 1][j - 1][k - 1] = (Uc[idx] * d1 + d[idx] * d3) * JJ[idx];
					df4[i - 1][j - 1][k - 1] = (Ur[idx] * d1 + d[idx] * d4) * JJ[idx];
					df5[i - 1][j - 1][k - 1] = (0.5 * (Ua[idx] * Ua[idx] + Uc[idx] * Uc[idx] + Ur[idx] * Ur[idx]) * d1 + d[idx] * Ua[idx] * d2 + d[idx] * Uc[idx] * d3 + d[idx] * Ur[idx] * d4 + 2.5 * d5) * JJ[idx];//**2
				}//for
			}//if
		}//for
	}//if
	//!Զ���߽�  �����һ������(*na0 + na == nx)�����һ�������(i = na)
	if (*na0 + na == nx) {
		i = na;
		for (k = 1; k <= nr; k++) {
			for (j = 1; j <= nc; j++) {
				idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
				dpda = (3.0 * p[idx] - 4.0 * p[idx - 1] + p[idx - 2]) / (2.0 * *da);
				ddda = (3.0 * d[idx] - 4.0 * d[idx - 1] + d[idx - 2]) / (2.0 * *da);
				duada = (3.0 * Ua[idx] - 4.0 * Ua[idx - 1] + Ua[idx - 2]) / (2.0 * *da);
				ducda = (3.0 * Uc[idx] - 4.0 * Uc[idx - 1] + Uc[idx - 2]) / (2.0 * *da);
				durda = (3.0 * Ur[idx] - 4.0 * Ur[idx - 1] + Ur[idx - 2]) / (2.0 * *da);
				if ((Ua[idx] - sv[idx]) <= 0.0) {
					phi1 = 0.0;
				}
				else {
					phi1 = (dadx[i - 1 + BE] * Ua[idx] - sv[idx] * dadx[i - 1 + BE]) * (-d[idx] * sv[idx] * duada + dpda);
				}
				if (Ua[idx] <= 0.0) {
					phi2 = 0.0;
					phi3 = 0.0;
					phi4 = 0.0;
				}
				else {
					phi2 = dadx[i - 1 + BE] * Ua[idx] * (sv[idx] * sv[idx] * ddda - dpda); //**2
					phi3 = dadx[i - 1 + BE] * Ua[idx] * (dadx[i - 1 + BE] * ducda);
					phi4 = dadx[i - 1 + BE] * Ua[idx] * (dadx[i - 1 + BE] * durda);
				}
				if ((Ua[idx] + sv[idx]) <= 0.0) {
					phi5 = 0.0;
				}
				else {
					phi5 = (dadx[i - 1 + BE] * Ua[idx] + sv[idx] * dadx[i - 1 + BE]) * (d[idx] * sv[idx] * duada + dpda);
				}
				d1 = (0.5 * (phi1 + phi5) + phi2) / (sv[idx] * sv[idx]);//**2
				d2 = (phi5 - phi1) / (2.0 * d[idx] * sv[idx]);
				d3 = phi3 / dadx[i - 1 + BE];
				d4 = phi4 / dadx[i - 1 + BE];
				d5 = 0.5 * (phi1 + phi5);
				df1[i - 1][j - 1][k - 1] = (d1)*JJ[idx];
				df2[i - 1][j - 1][k - 1] = (Ua[idx] * d1 + d[idx] * d2) * JJ[idx];
				df3[i - 1][j - 1][k - 1] = (Uc[idx] * d1 + d[idx] * d3) * JJ[idx];
				df4[i - 1][j - 1][k - 1] = (Ur[idx] * d1 + d[idx] * d4) * JJ[idx];
				df5[i - 1][j - 1][k - 1] = (0.5 * (Ua[idx] * Ua[idx] + Uc[idx] * Uc[idx] + Ur[idx] * Ur[idx]) * d1 + d[idx] * Ua[idx] * d2 + d[idx] * Uc[idx] * d3 + d[idx] * Ur[idx] * d4 + 2.5 * d5) * JJ[idx]; //**2
			}//for
		}//for
	}//if
	idx = 0;
	for (k = 0; k <= nr - 1; k++) {
		for (j = 0; j <= nc - 1; j++) {
			for (i = 0; i <= na - 1; i++) {
				if ((*na0 + i + 1 >= nx1 && *nr0 + k + 1 > 1) || (*nr0 + k + 1 >= nz2)) {  
					Lu1[idx] = df1[i][j][k];
					Lu2[idx] = df2[i][j][k];
					Lu3[idx] = df3[i][j][k];
					Lu4[idx] = df4[i][j][k];
					Lu5[idx] = df5[i][j][k];
				}
				idx++;
			}
		}
	}
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	k1 = 0;
	k3 = 0;
	n0 = 1;
	nn = nc;
	case_lb = 0;
	case_ub = 0;
	fp = (double*)malloc(5 * len_nc * sizeof(double));
	fn = (double*)malloc(5 * len_nc * sizeof(double));
	for (k = 1; k <= nr; k++) {
		dx = z[k - 1 + BE] * *dc;
		for (i = 1; i <= na; i++) {
			if ((*na0 + i >= nx1 && *nr0 + k > 1) || (*nr0 + k >= nz2)) {
				for (j = 1 - BE; j <= nc + BE; j++) {
					idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
					int idx_Lu = na * nc * (k-1) + na * (j-1-3) + (i-1);
					k2 = dcdy[j - 1 + BE] * JJ[idx];
					steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[j - 1 + BE], &fn[j - 1 + BE], len_nc);
					if(j>=n0+3 && j<=nn+3)
					{
						dfdx_p_v2_c(case_lb, case_ub, j-8, nn, dx, fp, len_nc);
						dfdx_n_v2_c(case_lb, case_ub, j-8, nn, dx, fn, len_nc);
						if(j==n0+4) continue;
						Lu1[idx_Lu] = Lu1[idx_Lu] + fp[j-3+BE-1] + fn[j-3+BE-1];
						Lu2[idx_Lu] = Lu2[idx_Lu] + fp[j-3+BE-1+len_nc*1] + fn[j-3+BE-1+len_nc*1];
						Lu3[idx_Lu] = Lu3[idx_Lu] + fp[j-3+BE-1+len_nc*2] + fn[j-3+BE-1+len_nc*2];
						Lu4[idx_Lu] = Lu4[idx_Lu] + fp[j-3+BE-1+len_nc*3] + fn[j-3+BE-1+len_nc*3];
						Lu5[idx_Lu] = Lu5[idx_Lu] + fp[j-3+BE-1+len_nc*4] + fn[j-3+BE-1+len_nc*4];
					}
				}
			}
		}
	}
	free(fp);
	free(fn);

		///----------------------------------------------------------------------------------------------------------------------------------------------------------------
	dx = *dr;
	k1 = 0;
	k2 = 0;
	nn = nr;
	fp = (double*)malloc(5 * len_nr * sizeof(double));
	fn = (double*)malloc(5 * len_nr * sizeof(double));
	fpx = (double*)malloc(5 * nr * sizeof(double));
	fnx = (double*)malloc(5 * nr * sizeof(double));
	for (j = 1; j <= nc; j++) {
		for (i = 1; i <= na; i++) {
			if (*na0 + i < nx1) {
				if (nz2 > * nr0 && nz2 <= *nr0 + nr) {
					for (k = 1; k <= nr; k++) {
						if (*nr0 + k == nz2) {
							n0 = k;
						}
					}
					case_lb = 1;
					if (*nr0 + nr == nz) {
						case_ub = 1;
					}
					else {
						case_ub = 0;
					}
				}
				else if (*nr0 > nz2) {
					n0 = 1;
					case_lb = 0;
					if (*nr0 + nr == nz) {
						case_ub = 1;
					}
					else {
						case_ub = 0;// !1�����������߽��ڴ�
					}
				}
			}
			else {
				n0 = 1;
				if (*nr0 == 0) {
					case_lb = 1;
					if (*nr0 + nr == nz) {
						case_ub = 1;
					}
					else {
						case_ub = 0;        // !1�����������߽��ڴ�
					}
				}
				else {
					case_lb = 0;
					if (*nr0 + nr == nz) {
						case_ub = 1;
					}
					else {
						case_ub = 0;
					}//!1�����������߽��ڴ�
				}
			}
			{//include "dfdz.f90"
				if (case_lb == 1 && case_ub == 1) {
					for (k = n0; k <= nn; k++) {
						//!k3 = dadx[k-1+BE] * JJ(i, j, k)
						k3 = 1.0 / dadx[i - 1 + BE];
						idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
						steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[k + BE - 1], &fn[k + BE - 1], len_nr);
					}
					dfdx_p(case_lb, case_ub, n0, nn, dx, fp, fpx);
					dfdx_n(case_lb, case_ub, n0, nn, dx, fn, fnx);
				}
				if (case_lb == 1 && case_ub == 0) {
					for (k = n0; k <= nn + 6; k++) {
						//!k3 = dadx[k-1+BE] * JJ(i, j, k)
						k3 = 1.0 / dadx[i - 1 + BE];
						idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
						steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[k + BE - 1], &fn[k + BE - 1], len_nr);
					}
					dfdx_p(case_lb, case_ub, n0, nn, dx, fp, fpx);
					dfdx_n(case_lb, case_ub, n0, nn, dx, fn, fnx);
				}
				if (case_lb == 0 && case_ub == 1) {
					for (k = n0 - 6; k <= nn; k++) {
						//!k3 = dadx[k-1+BE] * JJ(i, j, k)
						k3 = 1.0 / dadx[i - 1 + BE];
						idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
						steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[k + BE - 1], &fn[k + BE - 1], len_nr);
					}
					dfdx_p(case_lb, case_ub, n0, nn, dx, fp, fpx);
					dfdx_n(case_lb, case_ub, n0, nn, dx, fn, fnx);
				}
				if (case_lb == 0 && case_ub == 0) {
					for (k = n0 - 6; k <= nn + 6; k++) {
						//!k3 = dadx[k-1+BE] * JJ(i, j, k)
						k3 = 1.0 / dadx[i - 1 + BE];
						idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
						steger_warming(k1, k2, k3, d[idx], Ua[idx], Uc[idx], Ur[idx], sv[idx], &fp[k + BE - 1], &fn[k + BE - 1], len_nr);
					}
					dfdx_p(case_lb, case_ub, n0, nn, dx, fp, fpx);
					dfdx_n(case_lb, case_ub, n0, nn, dx, fn, fnx);
				}
				for (k = n0; k <= nn; k++) { 
					df1[i - 1][j - 1][k - 1] = fpx[k - 1] + fnx[k - 1];
					df2[i - 1][j - 1][k - 1] = fpx[k - 1 + nn] + fnx[k - 1 + nn];
					df3[i - 1][j - 1][k - 1] = fpx[k - 1 + nn * 2] + fnx[k - 1 + nn * 2];
					df4[i - 1][j - 1][k - 1] = fpx[k - 1 + nn * 3] + fnx[k - 1 + nn * 3];
					df5[i - 1][j - 1][k - 1] = fpx[k - 1 + nn * 4] + fnx[k - 1 + nn * 4];
				}
			}//include "dfdz.f90"
		}//for
	}//for
	free(fp);
	free(fn);
	free(fpx);
	free(fnx);
	////!���ڱ߽�
	if (*nr0 + nr == nz) {
		k = nr;
		for (j = 1; j <= nc; j++) {
			for (i = 1; i <= na; i++) {
				idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
				dddr = (3.0 * d[idx] - 4.0 * d[idx - len_na * len_nc] + d[idx - 2 * len_na * len_nc]) / (2.0 * *dr);
				dpdr = (3.0 * p[idx] - 4.0 * p[idx - len_na * len_nc] + p[idx - 2 * len_na * len_nc]) / (2.0 * *dr);
				duadr = (3.0 * Ua[idx] - 4.0 * Ua[idx - len_na * len_nc] + Ua[idx - 2 * len_na * len_nc]) / (2.0 * *dr);
				ducdr = (3.0 * Uc[idx] - 4.0 * Uc[idx - len_na * len_nc] + Uc[idx - 2 * len_na * len_nc]) / (2.0 * *dr);
				durdr = (3.0 * Ur[idx] - 4.0 * Ur[idx - len_na * len_nc] + Ur[idx - 2 * len_na * len_nc]) / (2.0 * *dr);
				if ((Ur[idx] - sv[idx]) <= 0.0) {
					phi1 = 0.0;
				}
				else {
					phi1 = (dadx[k - 1 + BE] * Ur[idx] - sv[idx] * dadx[k - 1 + BE]) * (-d[idx] * sv[idx] * durdr + dpdr);
				}
				if (Ur[idx] <= 0.0) {
					phi2 = 0.0;
					phi3 = 0.0;
					phi4 = 0.0;
				}
				else {
					phi2 = dadx[k - 1 + BE] * Ur[idx] * (sv[idx] * sv[idx] * dddr - dpdr); //**2
					phi3 = dadx[k - 1 + BE] * Ur[idx] * (dadx[k - 1 + BE] * duadr);
					phi4 = dadx[k - 1 + BE] * Ur[idx] * (dadx[k - 1 + BE] * ducdr);
				}
				if ((Ur[idx] + sv[idx]) <= 0.0) {
					phi5 = 0.0;
				}
				else {
					phi5 = (dadx[k - 1 + BE] * Ur[idx] + sv[idx] * dadx[k - 1 + BE]) * (d[idx] * sv[idx] * durdr + dpdr);
				}
				d1 = (0.5 * (phi1 + phi5) + phi2) / (sv[idx] * sv[idx]); //**2
				d2 = phi3 / dadx[k - 1 + BE];
				d3 = phi4 / dadx[k - 1 + BE];
				d4 = (phi5 - phi1) / (2.0 * d[idx] * sv[idx]);
				d5 = 0.5 * (phi1 + phi5);
				df1[i - 1][j - 1][k - 1] = (d1)*JJ[idx]; 
				df2[i - 1][j - 1][k - 1] = (Ua[idx] * d1 + d[idx] * d2) * JJ[idx];
				df3[i - 1][j - 1][k - 1] = (Uc[idx] * d1 + d[idx] * d3) * JJ[idx];
				df4[i - 1][j - 1][k - 1] = (Ur[idx] * d1 + d[idx] * d4) * JJ[idx];
				df5[i - 1][j - 1][k - 1] = (0.5 * (Ua[idx] * Ua[idx] + Uc[idx] * Uc[idx] + Ur[idx] * Ur[idx]) * d1 + d[idx] * Ua[idx] * d2 + d[idx] * Uc[idx] * d3 + d[idx] * Ur[idx] * d4 + 2.5 * d5) * JJ[idx]; //**2
			}//for
		}//for
	}//if
	//!��ճ����ķ���ֵ
	Lu_idx = 0;
	for (k = 1; k <= nr; k++) {
		for (j = 1; j <= nc; j++) {
			for (i = 1; i <= na; i++) {
				if ((*na0 + i >= nx1 && *nr0 + k > 1) || (*nr0 + k >= nz2)) { 
					idx = len_na * len_nc * (k + 5) + len_na * (j + 5) + (i + 5);
					Lu1[Lu_idx] = (Lu1[Lu_idx] + df1[i - 1][j - 1][k - 1]) / JJ[idx];
					Lu2[Lu_idx] = (Lu2[Lu_idx] + df2[i - 1][j - 1][k - 1]) / JJ[idx];
					Lu3[Lu_idx] = (Lu3[Lu_idx] + df3[i - 1][j - 1][k - 1]) / JJ[idx];
					Lu4[Lu_idx] = (Lu4[Lu_idx] + df4[i - 1][j - 1][k - 1]) / JJ[idx];
					Lu5[Lu_idx] = (Lu5[Lu_idx] + df5[i - 1][j - 1][k - 1]) / JJ[idx];
				}
				Lu_idx++;
			}
		}
	}
}

void steger_warming(const double k1, const double k2, const double k3, const double d,
	const double u, const double v, const double w, const double sv, double* fp, double* fn, int step) {

	double k1b, k2b, k3b, kk;
	double lambda[3], lamb[3];
	int i;
	const double gama = 1.4;
	const double gama2_2 = 0.8;
	const double gama2 = 2.8;
	const double gama_1 = 0.4;
	const double gama_w = 2.0;
	const double epsilon_2 = 1e-12;

	kk = sqrt((k1) * (k1)+(k2) * (k2)+(k3) * (k3));
	lambda[0] = (k1) * (u)+(k2) * (v)+(k3) * (w);
	lambda[1] = lambda[0] + (sv)*kk;
	lambda[2] = lambda[0] - (sv)*kk;
	k1b = (k1) / kk;
	k2b = (k2) / kk;
	k3b = (k3) / kk;
	double usv = (u)+(sv)*k1b;
	double vsv = (v)+(sv)*k2b;
	double wsv = (w)+(sv)*k3b;
	double usv_ = (u)-(sv)*k1b;
	double vsv_ = (v)-(sv)*k2b;
	double wsv_ = (w)-(sv)*k3b;
	for (int i = 0; i < 3; i++)
		lamb[i] = 0.5 * (lambda[i] + sqrt(lambda[i] * lambda[i] + epsilon_2));
	fp[0] = (d) / gama2 * (gama2_2 * lamb[0] + lamb[1] + lamb[2]);
	fp[step] = (d) / gama2 * (gama2_2 * lamb[0] * (u)+lamb[1] * usv + lamb[2] * usv_);
	fp[step * 2] = (d) / gama2 * (gama2_2 * lamb[0] * (v)+lamb[1] * vsv + lamb[2] * vsv_);
	fp[step * 3] = (d) / gama2 * (gama2_2 * lamb[0] * (w)+lamb[1] * wsv + lamb[2] * wsv_);
	fp[step * 4] = (d) / gama2 * (gama_1 * lamb[0] * ((u) * (u)+(v) * (v)+(w) * (w)) +
		0.5 * lamb[1] * (usv * usv + vsv * vsv + wsv * wsv) +
		0.5 * lamb[2] * (usv_ * usv_ + vsv_ * vsv_ + wsv_ * wsv_) +
		gama_w * (lamb[1] + lamb[2]) * (sv) * (sv));
	for (int i = 0; i < 3; i++)
		lamb[i] = 0.5 * (lambda[i] - sqrt(lambda[i] * lambda[i] + epsilon_2));
	fn[0] = (d) / gama2 * (gama2_2 * lamb[0] + lamb[1] + lamb[2]);
	fn[step] = (d) / gama2 * (gama2_2 * lamb[0] * (u)+lamb[1] * usv + lamb[2] * usv_);
	fn[step * 2] = (d) / gama2 * (gama2_2 * lamb[0] * (v)+lamb[1] * vsv + lamb[2] * vsv_);
	fn[step * 3] = (d) / gama2 * (gama2_2 * lamb[0] * (w)+lamb[1] * wsv + lamb[2] * wsv_);
	fn[step * 4] = (d) / gama2 * (gama_1 * lamb[0] * ((u) * (u)+(v) * (v)+(w) * (w)) +
		0.5 * lamb[1] * (usv * usv + vsv * vsv + wsv * wsv) +
		0.5 * lamb[2] * (usv_ * usv_ + vsv_ * vsv_ + wsv_ * wsv_) +
		gama_w * (lamb[1] + lamb[2]) * (sv) * (sv));
}


void dfdx_p(const int case_lb, const int case_ub, const int jj, const int nn, const double dx, double* fp, double* fpx) {
	int i, m;
	const int BE = 6;
	const int nn0 = nn;
	const int nn12 = nn + 2 * BE;
	double fhp[5][nn12 + 1];
	for (int m = 1 - 1; m <= 5 - 1; m++) {
		for (int i = (jj)+3 + BE - 1; i <= (nn)-3 + BE - 1; i++)
			fh_p_weno7_c_(&fp[m * nn12 + i - 3], &fp[m * nn12 + i - 2], &fp[m * nn12 + i - 1], &fp[m * nn12 + i], &fp[m * nn12 + i + 1], &fp[m * nn12 + i + 2], &fp[m * nn12 + i + 3], &fhp[m][i]);
		for (int i = (jj)+4 + BE - 1; i <= (nn)-3 + BE - 1; i++)
			fpx[m * nn0 + i - BE] = (fhp[m][i] - fhp[m][i - 1]) / dx;
		if (case_lb == 1) {
			i = (jj)+2 + BE - 1;
			fh_p_weno5_c_(&fp[m * nn12 + i - 2], &fp[m * nn12 + i - 1], &fp[m * nn12 + i], &fp[m * nn12 + i + 1], &fp[m * nn12 + i + 2], &fhp[m][i]);
			fpx[m * nn0 + i + 1 - BE] = (fhp[m][i + 1] - fhp[m][i]) / dx;
			i = (jj)+1 + BE - 1;
			fh_p_weno3_c_(&fp[m * nn12 + i - 1], &fp[m * nn12 + i], &fp[m * nn12 + i + 1], &fhp[m][i]);
			fpx[m * nn0 + i + 1 - BE] = (fhp[m][i + 1] - fhp[m][i]) / dx;
			fpx[m * nn0 + i - BE] = (fp[m * nn12 + i] - fp[m * nn12 + i - 1]) / dx;
			fpx[m * nn0 + i - 1 - BE] = (fp[m * nn12 + i] - fp[m * nn12 + i - 1]) / dx;
		}
		else {
			for (int i = (jj)+2 + BE - 1; i >= (jj)-1 + BE - 1; i--) {
				fh_p_weno7_c_(&fp[m * nn12 + i - 3], &fp[m * nn12 + i - 2], &fp[m * nn12 + i - 1], &fp[m * nn12 + i], &fp[m * nn12 + i + 1], &fp[m * nn12 + i + 2], &fp[m * nn12 + i + 3], &fhp[m][i]);
				fpx[m * nn0 + i + 1 - BE] = (fhp[m][i + 1] - fhp[m][i]) / dx;
			}
		}
		if (case_ub == 1) {
			fpx[m * nn0 + nn0 + BE - 1] = (fp[m * nn12 + nn0 + BE - 1] - fp[m * nn12 + nn0 - 1 + BE - 1]) / dx;
			i = (nn)-2 + BE - 1;
			fh_p_weno5_c_(&fp[m * nn12 + i - 2], &fp[m * nn12 + i - 1], &fp[m * nn12 + i], &fp[m * nn12 + i + 1], &fp[m * nn12 + i + 2], &fhp[m][i]);
			fpx[m * nn0 + i - BE] = (fhp[m][i] - fhp[m][i - 1]) / dx;
			i = (nn)-1 + BE - 1;
			fh_p_weno3_c_(&fp[m * nn12 + i - 1], &fp[m * nn12 + i], &fp[m * nn12 + i + 1], &fhp[m][i]);
			fpx[m * nn0 + i - BE] = (fhp[m][i] - fhp[m][i - 1]) / dx;
			i = (nn)+BE - 1;
			fpx[m * nn0 + i - BE] = (fp[m * nn12 + nn0 + BE - 1] - fp[m * nn12 + nn0 - 1 + BE - 1]) / dx;
		}
		else {
			for (int i = (nn)-2 + BE - 1; i <= (nn)+BE - 1; i++) {
				fh_p_weno7_c_(&fp[m * nn12 + i - 3], &fp[m * nn12 + i - 2], &fp[m * nn12 + i - 1], &fp[m * nn12 + i], &fp[m * nn12 + i + 1], &fp[m * nn12 + i + 2], &fp[m * nn12 + i + 3], &fhp[m][i]);
				fpx[m * nn0 + i - BE] = (fhp[m][i] - fhp[m][i - 1]) / dx;
			}
		}
	}
}


void dfdx_n(const int case_lb, const int case_ub, const int jj, const int nn, const double dx, double* fn, double* fnx) {
	int i, m;
	const int BE = 6;
	const int nn0 = nn;
	const int nn12 = nn + 2 * BE;
	double fhn[5][nn12 + 1];
	for (int m = 1 - 1; m <= 5 - 1; m++) {
		for (int i = (jj)+3 + BE - 1; i <= (nn)-3 + BE - 1; i++)
			fh_p_weno7_c_(&fn[m * nn12 + i + 3], &fn[m * nn12 + i + 2], &fn[m * nn12 + i + 1], &fn[m * nn12 + i], &fn[m * nn12 + i - 1], &fn[m * nn12 + i - 2], &fn[m * nn12 + i - 3], &fhn[m][i]);
		for (int i = (jj)+3 + BE - 1; i <= (nn)-4 + BE - 1; i++)
			fnx[m * nn0 + i - BE] = (fhn[m][i + 1] - fhn[m][i]) / dx;
		if (case_lb == 1) {
			i = (jj)+2 + BE - 1;
			fh_p_weno5_c_(&fn[m * nn12 + i + 2], &fn[m * nn12 + i + 1], &fn[m * nn12 + i], &fn[m * nn12 + i - 1], &fn[m * nn12 + i - 2], &fhn[m][i]);
			fnx[m * nn0 + i - BE] = (fhn[m][i + 1] - fhn[m][i]) / dx;
			i = (jj)+1 + BE - 1;
			fh_p_weno3_c_(&fn[m * nn12 + i + 1], &fn[m * nn12 + i], &fn[m * nn12 + i - 1], &fhn[m][i]);
			fnx[m * nn0 + i - BE] = (fhn[m][i + 1] - fhn[m][i]) / dx;
			i = (jj)+BE - 1;
			fnx[m * nn0 + i - BE] = (fn[m * nn12 + i + 1] - fn[m * nn12 + i]) / dx;
		}
		else {
			for (int i = (jj)+BE - 1; i <= (jj)+2 + BE - 1; i++) {
				fh_p_weno7_c_(&fn[m * nn12 + i + 3], &fn[m * nn12 + i + 2], &fn[m * nn12 + i + 1], &fn[m * nn12 + i], &fn[m * nn12 + i - 1], &fn[m * nn12 + i - 2], &fn[m * nn12 + i - 3], &fhn[m][i]);
			}
			for (int i = (jj)+BE - 1; i <= (jj)+2 + BE - 1; i++) {
				fnx[m * nn0 + i - BE] = (fhn[m][i + 1] - fhn[m][i]) / dx;
			}
		}
		if (case_ub == 1) {
			i = (nn)-2 + BE - 1;
			fh_p_weno5_c_(&fn[m * nn12 + i + 2], &fn[m * nn12 + i + 1], &fn[m * nn12 + i], &fn[m * nn12 + i - 1], &fn[m * nn12 + i - 2], &fhn[m][i]);
			fnx[m * nn0 + i - 1 - BE] = (fhn[m][i] - fhn[m][i - 1]) / dx;
			i = (nn)-1 + BE - 1;
			fh_p_weno3_c_(&fn[m * nn12 + i + 1], &fn[m * nn12 + i], &fn[m * nn12 + i - 1], &fhn[m][i]);
			fnx[m * nn0 + i - 1 - BE] = (fhn[m][i] - fhn[m][i - 1]) / dx;
			i = (nn)+BE - 1;  
			fnx[m * nn0 + i - 1 - BE] = (fn[m * nn12 + nn0 + BE - 1] - fn[m * nn12 + nn0 - 1 + BE - 1]) / dx;
			fnx[m * nn0 + i - BE] = (fn[m * nn12 + nn0 + BE - 1] - fn[m * nn12 + nn0 - 1 + BE - 1]) / dx;
		}
		else {
			for (int i = (nn)-2 + BE - 1; i <= (nn)+1 + BE - 1; i++) { 
				fh_p_weno7_c_(&fn[m * nn12 + i + 3], &fn[m * nn12 + i + 2], &fn[m * nn12 + i + 1], &fn[m * nn12 + i], &fn[m * nn12 + i - 1], &fn[m * nn12 + i - 2], &fn[m * nn12 + i - 3], &fhn[m][i]);
			}
			for (int i = (nn)-2 + BE - 1; i <= (nn)+1 + BE - 1; i++) {
				fnx[m * nn0 + i - 1 - BE] = (fhn[m][i] - fhn[m][i - 1]) / dx;
			}
		}
	}
}


void fh_p_weno3_c_(const double* f0, const double* f1, const double* f2, double* fhp) { //positively biased
	int i;
	double w[2], beta[2], alfa[2], sumalfa;
	double fh[2];
	double d[2] = { 2.0 / 3.0, 1.0 / 3.0 };
	double epsilon = 1e-6;

	fh[0] = 0.5 * *f1 + 0.5 * *f2;
	fh[1] = 1.5 * *f1 - 0.5 * *f0;

	beta[0] = (*f2 - *f1) * (*f2 - *f1);
	beta[1] = (*f1 - *f0) * (*f1 - *f0);

	for (int i = 0; i <= 1; i++) {
		alfa[i] = d[i] / ((epsilon + beta[i]) * (epsilon + beta[i]));
	}
	sumalfa = alfa[0] + alfa[1];
	for (int i = 0; i <= 1; i++) {
		w[i] = alfa[i] / sumalfa;
	}

	*fhp = 0.0;
	for (int i = 0; i <= 1; i++) {
		*fhp = *fhp + w[i] * fh[i];
	}

}

void fh_p_weno5_c_(const double* f_2, const double* f_1, const double* f0, const double* f1, const double* f2, double* fhp) { //positively biased
	int i;
	double w[3], beta[3], alfa[3], sumalfa;
	double fh[3];
	double d[3] = { 0.3,0.6,0.1 };
	double epsilon = 1e-6;

	fh[0] = 1.0 / 3.0 * *f0 + 5.0 / 6.0 * *f1 - 1.0 / 6.0 * *f2;
	fh[1] = -1.0 / 6.0 * *f_1 + 5.0 / 6.0 * *f0 + 1.0 / 3.0 * *f1;
	fh[2] = 1.0 / 3.0 * *f_2 - 7.0 / 6.0 * *f_1 + 11.0 / 6.0 * *f0;

	beta[0] = 13.0 / 12.0 * (*f0 - 2.0 * *f1 + *f2) * (*f0 - 2.0 * *f1 + *f2) + 0.25 * (3.0 * *f0 - 4.0 * *f1 + *f2) * (3.0 * *f0 - 4.0 * *f1 + *f2);
	beta[1] = 13.0 / 12.0 * (*f_1 - 2.0 * *f0 + *f1) * (*f_1 - 2.0 * *f0 + *f1) + 0.25 * (*f_1 - *f1) * (*f_1 - *f1);
	beta[2] = 13.0 / 12.0 * (*f_2 - 2.0 * *f_1 + *f0) * (*f_2 - 2.0 * *f_1 + *f0) + 0.25 * (*f_2 - 4.0 * *f_1 + 3.0 * *f0) * (*f_2 - 4.0 * *f_1 + 3.0 * *f0);

	for (int i = 0; i <= 2; i++) {
		alfa[i] = d[i] / ((epsilon + beta[i]) * (epsilon + beta[i]));
	}
	sumalfa = alfa[0] + alfa[1] + alfa[2];
	for (int i = 0; i <= 2; i++) {
		w[i] = alfa[i] / sumalfa;
	}

	*fhp = 0.0;
	for (int i = 0; i <= 2; i++) {
		*fhp = *fhp + w[i] * fh[i];
	}

}

void fh_p_weno7_c_(const double* f_3, const double* f_2, const double* f_1, const double* f0, const double* f1, const double* f2, const double* f3, double* fhp) { //positively biased
	int i;
	double w[4], beta[4], alfa[4], sumalfa;
	double fh[4];
	double d[4];
	double TV[4], TVR, TV_MAX, TV_MIN;
	double epsilon = 1e-6;

	d[0] = 1.0 / 35.0;//fortran: d(0)=1.0d0/35.0d0
	d[1] = 12.0 / 35.0;
	d[2] = 18.0 / 35.0;
	d[3] = 4.0 / 35.0;

	fh[0] = -1.0 / 4.0 * *f_3 + 13.0 / 12.0 * *f_2 - 23.0 / 12.0 * *f_1 + 25.0 / 12.0 * *f0;
	fh[1] = 1.0 / 12.0 * *f_2 - 5.0 / 12.0 * *f_1 + 13.0 / 12.0 * *f0 + 1.0 / 4.0 * *f1;
	fh[2] = -1.0 / 12.0 * *f_1 + 7.0 / 12.0 * *f0 + 7.0 / 12.0 * *f1 - 1.0 / 12.0 * *f2;
	fh[3] = 1.0 / 4.0 * *f0 + 13.0 / 12.0 * *f1 - 5.0 / 12.0 * *f2 + 1.0 / 12.0 * *f3;

	TV[0] = fabs(*f_3 - *f_2) + fabs(*f_2 - *f_1) + fabs(*f_1 - *f0);
	TV[1] = fabs(*f_2 - *f_1) + fabs(*f_1 - *f0) + fabs(*f0 - *f1);
	TV[2] = fabs(*f_1 - *f0) + fabs(*f0 - *f1) + fabs(*f1 - *f2);
	TV[3] = fabs(*f0 - *f1) + fabs(*f1 - *f2) + fabs(*f2 - *f3);

	TV_MAX = TV[0];
	TV_MIN = TV[0];
	for (int i = 1; i < 4; i++) {    //4 is the size of TV
		TV_MAX = TV_MAX >= TV[i] ? TV_MAX : TV[i];
		TV_MIN = TV_MIN <= TV[i] ? TV_MIN : TV[i];
	}


	TVR = TV_MAX / (TV_MIN + epsilon);

	if (TV_MAX < 0.2 && TVR < 5.0) {
		for (int i = 0; i <= 3; i++) {
			w[i] = d[i];
		}
	}

	else {
		beta[0] = *f_3 * (547.0 * *f_3 - 3882.0 * *f_2 + 4642.0 * *f_1 - 1854.0 * *f0)
			+ *f_2 * (7043.0 * *f_2 - 17246.0 * *f_1 + 7042.0 * *f0)
			+ *f_1 * (11003.0 * *f_1 - 9402.0 * *f0)
			+ 2107.0 * (*f0 * *f0);
		beta[1] = *f_2 * (267.0 * *f_2 - 1642.0 * *f_1 + 1602.0 * *f0 - 494.0 * *f1)
			+ *f_1 * (2843.0 * *f_1 - 5966.0 * *f0 + 1922.0 * *f1)
			+ *f0 * (3443.0 * *f0 - 2522.0 * *f1)
			+ 547.0 * (*f1 * *f1);
		beta[2] = *f_1 * (547.0 * *f_1 - 2522.0 * *f0 + 1922.0 * *f1 - 494.0 * *f2)
			+ *f0 * (3443.0 * *f0 - 5966.0 * *f1 + 1602.0 * *f2)
			+ *f1 * (2843.0 * *f1 - 1642.0 * *f2)
			+ 267.0 * (*f2 * *f2);
		beta[3] = *f0 * (2107.0 * *f0 - 9402.0 * *f1 + 7042.0 * *f2 - 1854.0 * *f3)
			+ *f1 * (11003.0 * *f1 - 17246.0 * *f2 + 4642.0 * *f3)
			+ *f2 * (7043.0 * *f2 - 3882.0 * *f3)
			+ 547.0 * (*f3 * *f3);

		for (int i = 0; i <= 3; i++) {
			alfa[i] = d[i] / ((epsilon + beta[i]) * (epsilon + beta[i])); //alfa(i)=d(i)/(epsilon+beta(i))**2
		}

		sumalfa = alfa[0] + alfa[1] + alfa[2] + alfa[3];

		for (int i = 0; i <= 3; i++) {
			w[i] = alfa[i] / sumalfa;
		}
	}
	*fhp = 0.0;
	for (int i = 0; i <= 3; i++) {
		*fhp = *fhp + w[i] * fh[i];
	}
}


void dfdx_p_v2_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double *dx, double fp[], int *pstep){
	int i;
	const int BE = 6;
	const int nn12 = *nn+2*BE;  
	double fhp[5][2];
	int step = *pstep;
	int flag = 0;
	double _dx = 1.0 / *dx;

		for (int m = 0; m <= 4; m++) {
			i = (*jj)-1+BE-1;
			fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[m][flag=flag^1]);
			i = (*jj)+BE-1;
			fh_p_weno7_c_(&fp[m*nn12+i-3],&fp[m*nn12+i-2],&fp[m*nn12+i-1],&fp[m*nn12+i],&fp[m*nn12+i+1],&fp[m*nn12+i+2],&fp[m*nn12+i+3],&fhp[flag=flag^1]);
			fp[m*nn12+i-4] = (fhp[m][flag]  - fhp[m][flag^1]  ) ;
		}
}

void dfdx_n_v2_c_(const int *case_lb, const int *case_ub, const int *jj, const int *nn, const double *dx, double fp[], int *pstep){
	int i;
	const int BE = 6;
	const int nn12 = *nn+2*BE;  
	double fhp[5][2];
	int step = *pstep;
	int flag = 0;
	double _dx = 1.0 / *dx;

		for (int m = 0; m <= 4; m++) {
			i = (*jj)-1+BE-1;
			fh_p_weno7_c_(&fp[m*nn12+i+3],&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fp[m*nn12+i-3],&fhp[m][flag=flag^1]);
			i = (*jj)+BE-1;
			fh_p_weno7_c_(&fp[m*nn12+i+3],&fp[m*nn12+i+2],&fp[m*nn12+i+1],&fp[m*nn12+i],&fp[m*nn12+i-1],&fp[m*nn12+i-2],&fp[m*nn12+i-3],&fhp[m][flag=flag^1]);
			fp[m*nn12+i-4] = (fhp[m][flag]  - fhp[m][flag^1]  ) ;
		}
}