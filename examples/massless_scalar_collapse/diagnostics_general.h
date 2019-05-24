#ifndef _DIAGNOSTICS_GENERAL_H_
#define _DIAGNOSTICS_GENERAL_H_

void compute_diagnostics_general(
	int Nx, int exc_jC, 
	double s_L,
	double dt, double dx,
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* mass_aspect,
	double* ingoing_null_characteristic,
	double* outgoing_null_characteristic,
	double* Ricci_scalar,
	double* Gauss_Bonnet_scalar)
;

#endif /* _DIAGNOSTICS_GENERAL_H_ */
