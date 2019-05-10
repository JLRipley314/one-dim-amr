#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static inline double max_fabs(double var_1, double var_2) 
{
	return (fabs(var_1)>fabs(var_2)) ? fabs(var_1) : fabs(var_2) ;
}
/*****************************************************************************
 * Dirichlet boundary conditions. 
 ****************************************************************************/
void advance_tStep_wave(
	int Nx,
	double dt, 	double dx,
	double* P_n, 	double* P_nm1,
	double* Q_n,	double* Q_nm1)
{
	double t_der_P, x_der_P, t_der_Q, x_der_Q ;
	double res_P, jac_P, res_Q, jac_Q ;

	double res_infty_norm = 0 ;
	do {
		for (int iC=1; iC<Nx-1; iC++) {
			t_der_P = (P_n[iC]-P_nm1[iC])/dt ;
			t_der_Q = (Q_n[iC]-Q_nm1[iC])/dt ;

			r_der_P  = (P_n[iC+1]  -P_n[iC-1]  )/(2.*dr) ;
			r_der_P += (P_nm1[iC+1]-P_nm1[iC-1])/(2.*dr) ;
			r_der_P /= 2. ;

			r_der_Q  = (Q_n[iC+1]  -Q_n[iC-1]  )/(2.*dr) ;
			r_der_Q += (Q_nm1[iC+1]-Q_nm1[iC-1])/(2.*dr) ;
			r_der_Q /= 2. ;

			res_P = t_der_P - r_der_Q ;
			res_Q = t_der_Q - r_der_P ;

			jac_P = (1/dt) ;
			jac_Q = (1/dt) ;

			P_n[iC] -= res_P/jac_P ;
			Q_n[iC] -= res_Q/jac_Q ;

			res_infty_norm = max_fabs(res_infty_norm,res_Q) ;
			res_infty_norm = max_fabs(res_infty_norm,res_P) ;
		}
	while (res_infty_norm > 1e-10) ;

	return ;
}
/*****************************************************************************
 ****************************************************************************/
void Kreiss_Oliger_Dissipation(
	int Nx,
	double* field)
{
	for (int iC=2; iC<Nx-2; iC++) {
		field[iC] -= 
	}

	return ;
}
