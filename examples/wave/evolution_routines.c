#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define ERR_TOLERANCE ((double)1e-10)

/****************************************************************************/
static inline double max_fabs(double var_1, double var_2) 
{
	return (fabs(var_1)>fabs(var_2)) ? fabs(var_1) : fabs(var_2) ;
}
/****************************************************************************/
void copy_to_2nd_array(int Nx, double* array_1, double* array_2) 
{
	for (int iC=0; iC<Nx; iC++) {
		array_2[iC] = array_1[iC] ;
	}
	return ;
}
/*****************************************************************************
 * leftgoing Gaussian pulse 
 ****************************************************************************/
void initial_Data(
	int Nx, 	double dx,
	double left_point,
	double* P_n, 	double* Q_n)
{
	double amp = 1 ;
	double width = 5 ;
	double x_0 = 0 ;
	double x = 0 ;

	for (int iC=0; iC<Nx; iC++) {
		x = (iC * dx) - left_point ;
		Q_n[iC] = amp * (-(x-x_0)/pow(width,2)) * exp(-pow((x-x_0)/width,2)) ;
		P_n[iC] = - Q_n[iC] ;
	}
	return ;
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

			x_der_P  = (P_n[iC+1]  -P_n[iC-1]  )/(2.*dx) ;
			x_der_P += (P_nm1[iC+1]-P_nm1[iC-1])/(2.*dx) ;
			x_der_P /= 2. ;

			x_der_Q  = (Q_n[iC+1]  -Q_n[iC-1]  )/(2.*dx) ;
			x_der_Q += (Q_nm1[iC+1]-Q_nm1[iC-1])/(2.*dx) ;
			x_der_Q /= 2. ;

			res_P = t_der_P - x_der_Q ;
			res_Q = t_der_Q - x_der_P ;

			jac_P = (1/dt) ;
			jac_Q = (1/dt) ;

			P_n[iC] -= res_P/jac_P ;
			Q_n[iC] -= res_Q/jac_Q ;

			res_infty_norm = max_fabs(res_infty_norm,res_Q) ;
			res_infty_norm = max_fabs(res_infty_norm,res_P) ;
		}
	} while (res_infty_norm > ERR_TOLERANCE) ;

	return ;
}
/*****************************************************************************
 ****************************************************************************/
void Kreiss_Oliger_Filter(
	int Nx,
	double* field)
{
	double epsilon_ko = 0.5 ;
	for (int iC=2; iC<Nx-2; iC++) {
		field[iC] -= epsilon_ko * (
			field[iC+2] + (-4.*field[iC+1]) + (6.*field[iC]) + (-4.*field[iC-1]) + field[iC-2] 
		)
		;
	}

	return ;
}
/*****************************************************************************
 ****************************************************************************/
void save_to_txt_file(
	int Nx,
	FILE* output_file,
	double* field)
{
	for (int iC=0; iC<Nx; iC++) {
		fprintf(output_file, "%.10f;", field[iC]) ;
	}
	fprintf(output_file, "\n") ;
	return ;
}
