#ifndef _DIAGNOSTICS_EDGB_H_
#define _DIAGNOSTICS_EDGB_H_

void compute_diagnostics_massless_scalar_EdGB(
	int Nx, int exc_jC, 
	double s_L,double c_gbs,
	double dt, double dx,
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* P_n, double* P_nm1, double* P_nm2,
	double* Q_n, double* Q_nm1, double* Q_nm2,
	double* eom_TR,
	double* eom_ThTh,
	double* eom_scalar)
;

#endif /* _DIAGNOSTICS_EDGB_H_ */

