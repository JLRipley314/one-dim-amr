#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "evolution_routines.h"

#define ERR_TOLERANCE ((double)1e-10)

#define S_L 100

/*==========================================================================*/
inline double D1_CrankNicolson_2ndOrder(
	double f_np1, double f_n, double dt)
{
	return (f_np1-f_n)/dt ;
}
inline double D1_center_2ndOrder(
	double f_np1, double f_nm1, double dx)
{
	return (f_np1-f_nm1)/(2*dx) ;
}
inline double D1_forward_2ndOrder(
	double f_np2, double f_np1, double f_n, double dx)
{
	return ((-0.5*f_np2)+(2*f_np1)+(-1.5*f_n))/dx ;
}
/*==========================================================================*/
inline double stereographic_r(double s_L, double x_j)
{
        return pow(1. - (x_j/s_L), -1) * x_j ;
}
inline double stereographic_dr(double s_L, double x_j, double dx)
{
        return pow(1. - (x_j/s_L), -2) * dx ;
}
/*==========================================================================*/
inline double compute_weighted_infty_norm(double weight, double val_1, double val_2)
{
        return (fabs(val_2)>fabs(val_1)) ? fabs(weight*val_2) : fabs(weight*val_1) ;
}
/*==========================================================================*/
static inline double max_fabs(double var_1, double var_2) 
{
	return (fabs(var_1)>fabs(var_2)) ? fabs(var_1) : fabs(var_2) ;
}
/*==========================================================================*/
static void set_array_val(int Nx, double val, double* array) 
{
	for (int iC=0; iC<Nx; iC++) {
		array[iC] = val ;
	}
	return ;
}/*==========================================================================*/
static void copy_to_2nd_array(int Nx, double* array_1, double* array_2) 
{
	for (int iC=0; iC<Nx; iC++) {
		array_2[iC] = array_1[iC] ;
	}
	return ;
}
/*==========================================================================*/
/* leftgoing Gaussian pulse */
/*==========================================================================*/
void initial_data_Gaussian(
	int Nx, 	double dx,
	double bbox[2],
	double* Al_n, double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double*  P_n, double*  P_nm1, double*  Q_n, double*  Q_nm1)
{
	double left_point = bbox[0] ;

	double amp = 1 ;
	double width = 5 ;
	double r_0 = 10 ;
	double x = 0 ;
	double r = 0 ;

	set_array_val(Nx, 1., Al_n  ) ;
	set_array_val(Nx, 0., Ze_n  ) ;

	for (int iC=0; iC<Nx; iC++) {
		x = (iC * dx) + left_point ;
		r = stereographic_r(S_L, x) ;
		Q_n[iC] = amp * exp(-pow((r-r_0)/width,2)) * (
			(-(r-r_0)/pow(width,2)) * pow(r,2) 
		+	2*r
		) ;
		P_n[iC] = Q_n[iC] ;
	}
	copy_to_2nd_array(Nx, Al_n, Al_nm1) ;
	copy_to_2nd_array(Nx, Ze_n, Ze_nm1) ;
	copy_to_2nd_array(Nx,  P_n,  P_nm1) ;
	copy_to_2nd_array(Nx,  Q_n,  Q_nm1) ;
	return ;
}
/*==========================================================================*/
void Kreiss_Oliger_Filter(
	int Nx,
	double* field)
{
	double epsilon_ko = 0.5 ;
	for (int iC=2; iC<Nx-2; iC++) {
		field[iC] -= (epsilon_ko/16.) * (
			field[iC+2] 
		+ 	(-4.*field[iC+1]) 
		+ 	(6.*field[iC]) 
		+ 	(-4.*field[iC-1]) 
		+ 	field[iC-2] 
		)
		;
	}

	return ;
}
/*****************************************************************************
 * Dirichlet boundary conditions. 
 ****************************************************************************/
static double compute_iteration_GR_Crank_Nicolson_PQ(
	int Nx,
	double dt, 	double dx,
	double bbox[2],
	bool perim_interior[2],
	double* Al_n, 	double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double*  P_n, 	double*  P_nm1, double*  Q_n, double*  Q_nm1)
{
	int exc_jC = 0 ;

	double lower_x = perim_interior[0] * dx ;

	double  
		x_jm1, x_j, x_jp1, x_jp2, 
		r_jm1, r_j, r_jp1, r_jp2, dr 
	;
	double res_Q = 0 ;
	double res_P = 0 ;
	double jac_Q = 1 ;
	double jac_P = 1 ;

	double res_infty_norm = 0 ; /* returning this */
/****************************************************************************/
/* interior */
/****************************************************************************/
	for (int jC=exc_jC+1;jC<Nx-1;jC++) {
		x_j   = lower_x + (dx * (jC)  ) ;
		x_jp1 = lower_x + (dx * (jC+1)) ;
		x_jm1 = lower_x + (dx * (jC-1)) ;

		r_j   = stereographic_r(S_L, x_j  ) ;
		r_jp1 = stereographic_r(S_L, x_jp1) ;
		r_jm1 = stereographic_r(S_L, x_jm1) ;

		dr = stereographic_dr(S_L, x_j, dx) ;
	/* Q field */
		res_Q = D1_CrankNicolson_2ndOrder(
			Q_n[jC], 
			Q_nm1[jC], 
			dt)
		;
		res_Q -= (1./2.)*D1_center_2ndOrder(
			Al_n[jC+1]*(P_n[jC+1] + Ze_n[jC+1]*Q_n[jC+1]),
			Al_n[jC-1]*(P_n[jC-1] + Ze_n[jC-1]*Q_n[jC-1]),
			dr)
		;
		res_Q -= (1./2.)*D1_center_2ndOrder(
			Al_nm1[jC+1]*(P_nm1[jC+1] + Ze_nm1[jC+1]*Q_nm1[jC+1]),
			Al_nm1[jC-1]*(P_nm1[jC-1] + Ze_nm1[jC-1]*Q_nm1[jC-1]),
			dr)
		;
		jac_Q = (1./dt) ;
	/* P field */
		res_P = 
			D1_CrankNicolson_2ndOrder(
			P_n[jC], 
			P_nm1[jC], 
			dt)
		;
		res_P -= (1./2.)*pow(r_j,-2)*D1_center_2ndOrder(
			pow(r_jp1,2)*Al_n[jC+1]*(Q_n[jC+1] + Ze_n[jC+1]*P_n[jC+1]),
			pow(r_jm1,2)*Al_n[jC-1]*(Q_n[jC-1] + Ze_n[jC-1]*P_n[jC-1]),
			dr)
		;
		res_P -= (1./2.)*pow(r_j,-2)*D1_center_2ndOrder(
			pow(r_jp1,2)*Al_nm1[jC+1]*(Q_nm1[jC+1] + Ze_nm1[jC+1]*P_nm1[jC+1]),
			pow(r_jm1,2)*Al_nm1[jC-1]*(Q_nm1[jC-1] + Ze_nm1[jC-1]*P_nm1[jC-1]),
			dr)
		;
		jac_P = (1./dt) ;
	/****/
		Q_n[jC] -= res_Q / jac_Q ;
		P_n[jC] -= res_P / jac_P ;

		res_infty_norm = compute_weighted_infty_norm(1-x_j/S_L, res_Q, res_infty_norm) ;
		res_infty_norm = compute_weighted_infty_norm(1-x_j/S_L, res_P, res_infty_norm) ;
	}
/****************************************************************************/
/* lower */
/* Q[0] = 0, so do not change anything. r_Der_P=0 at r=0 */	
/****************************************************************************/
	if ((exc_jC == 0) 
	) {
		x_j = 0 ;
		dr  = pow(1. - (x_j/S_L), -2) * dx;

		res_P = D1_forward_2ndOrder(
			P_n[2], P_n[1], P_n[0], 
		dr) ;
		jac_P = - (3./2.) / dr ;

		P_n[0] -= res_P / jac_P ;
	}
	if (exc_jC > 0) {
		x_j   = lower_x + (dx * (exc_jC))   ;
		x_jp1 = lower_x + (dx * (exc_jC+1)) ;
		x_jp2 = lower_x + (dx * (exc_jC+2)) ;

		r_j   = stereographic_r(S_L, x_j  ) ;
		r_jp1 = stereographic_r(S_L, x_jp1) ;
		r_jp2 = stereographic_r(S_L, x_jp2) ;

		dr = stereographic_dr(S_L, x_j, dx) ;
	/* Q field */
		res_Q = D1_CrankNicolson_2ndOrder(
			Q_n[exc_jC], 
			Q_nm1[exc_jC], 
			dt)
		;
		res_Q -= (1./2.)*D1_forward_2ndOrder(
			Al_n[exc_jC+2]*(P_n[exc_jC+2] + Ze_n[exc_jC+2]*Q_n[exc_jC+2]),
			Al_n[exc_jC+1]*(P_n[exc_jC+1] + Ze_n[exc_jC+1]*Q_n[exc_jC+1]),
			Al_n[exc_jC+0]*(P_n[exc_jC+0] + Ze_n[exc_jC+0]*Q_n[exc_jC+0]),
			dr)
		;
		res_Q -= (1./2.)*D1_forward_2ndOrder(
			Al_nm1[exc_jC+2]*(P_nm1[exc_jC+2] + Ze_nm1[exc_jC+2]*Q_nm1[exc_jC+2]),
			Al_nm1[exc_jC+1]*(P_nm1[exc_jC+1] + Ze_nm1[exc_jC+1]*Q_nm1[exc_jC+1]),
			Al_nm1[exc_jC+0]*(P_nm1[exc_jC+0] + Ze_nm1[exc_jC+0]*Q_nm1[exc_jC+0]),
			dr)
		;
		jac_Q =	1/dt 
		- 	(1./2.) * (-3./(2*dr)) * Al_n[exc_jC]*Ze_n[exc_jC] 
		;
	/* P field */
		res_P = D1_CrankNicolson_2ndOrder(
			P_n[exc_jC], 
			P_nm1[exc_jC], 
			dt)
		;
		res_P -= (1./2.)*D1_forward_2ndOrder(
			pow(r_jp2,2)*Al_n[exc_jC+2]*(Q_n[exc_jC+2] + Ze_n[exc_jC+2]*P_n[exc_jC+2]),
			pow(r_jp1,2)*Al_n[exc_jC+1]*(Q_n[exc_jC+1] + Ze_n[exc_jC+1]*P_n[exc_jC+1]),
			pow(r_j,  2)*Al_n[exc_jC+0]*(Q_n[exc_jC+0] + Ze_n[exc_jC+0]*P_n[exc_jC+0]),
			dr)
		;
		res_P -= (1./2.)*D1_forward_2ndOrder(
			pow(r_jp2,2)*Al_nm1[exc_jC+2]*(Q_nm1[exc_jC+2] + Ze_nm1[exc_jC+2]*P_nm1[exc_jC+2]),
			pow(r_jp1,2)*Al_nm1[exc_jC+1]*(Q_nm1[exc_jC+1] + Ze_nm1[exc_jC+1]*P_nm1[exc_jC+1]),
			pow(r_j  ,2)*Al_nm1[exc_jC+0]*(Q_nm1[exc_jC+0] + Ze_nm1[exc_jC+0]*P_nm1[exc_jC+0]),
			dr)
		;
		jac_P = 1./dt
		-	(1./2.) * (-3./(2*dr)) * pow(r_j,2) * Al_n[exc_jC] * Ze_n[exc_jC]
		;
	/****/
		Q_n[exc_jC] -= res_Q / jac_Q ;
		P_n[exc_jC] -= res_P / jac_P ;

		res_infty_norm = compute_weighted_infty_norm(1-x_j/S_L, res_Q, res_infty_norm) ;
		res_infty_norm = compute_weighted_infty_norm(1-x_j/S_L, res_P, res_infty_norm) ;
	}
/***************************************************************************/
/* dirichlet outer boundary conditions if outermost level;
 * otherwise do not evolve outer boundary anyways (it is interpolated) */
/***************************************************************************/
/***************************************************************************/
	return res_infty_norm ;
} 
void advance_tStep_massless_scalar(
	int Nx, 
	double dt, double dx, double bbox[2], 
	bool perim_interior[2],
	double* Al_n, double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double*  P_n, double*  P_nm1, double*  Q_n, double*  Q_nm1)
{ 
	copy_to_2nd_array(Nx, Al_n, Al_nm1) ;
	copy_to_2nd_array(Nx, Ze_n, Ze_nm1) ;
	copy_to_2nd_array(Nx,  P_n,  P_nm1) ;
	copy_to_2nd_array(Nx,  Q_n,  Q_nm1) ;

	compute_iteration_GR_Crank_Nicolson_PQ(
		Nx,
		dt, dx,
		bbox,
		perim_interior,
		Al_n, Al_nm1, Ze_n, Ze_nm1,
		 P_n,  P_nm1,  Q_n,  Q_nm1)
	;
	Kreiss_Oliger_Filter(Nx, P_n) ;
	Kreiss_Oliger_Filter(Nx, Q_n) ;

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
		fprintf(output_file, "%.10f\t", field[iC]) ;
	}
	fprintf(output_file, "\n") ;
	return ;
}
