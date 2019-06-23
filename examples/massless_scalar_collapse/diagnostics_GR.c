#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "stencils.h"
#include "diagnostics_GR.h"
#include "diagnostics_general.h"

/*==========================================================================*/
static void set_array_val(int start, int end, double val, double* array) 
{
	if (end<start) return ;
	for (int iC=start; iC<end; iC++) {
		array[iC] = val ;
	}
	return ;
}
/*==========================================================================*/
static void compute_SE_LL_components_massless_scalar(
	int Nx, double dx,
	double s_L,   int exc_jC,
	double* Al_n, double* Ze_n,
	double* P_n,  double* Q_n,
	double* SE_LL_TR,
	double* SE_LL_ThTh)
{
	for (int jC=exc_jC; jC<Nx; jC++) {
		double x_j = dx*jC ;
		double r_j = stereographic_r(s_L, x_j) ;

		double Al = Al_n[jC] ;
		double Ze = Ze_n[jC] ;

		double P = P_n[jC] ;
		double Q = Q_n[jC] ;

		SE_LL_TR[jC] =  (Al*(2*P*Q + (pow(P,2) + pow(Q,2))*Ze))/2. ;

		SE_LL_ThTh[jC] = (pow(r_j,2)*(pow(P,2) - pow(Q,2)))/2. ;
	}
	return ;
}
/*==========================================================================*/
static void compute_eom_TR(
	int Nx, 
	double dt, double dx,
	double s_L, int exc_jC, 
	double* Al_n, 
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* SE_LL_TR,
	double* eom_TR)
{
	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		double x_j = jC * dx ;
                double r_j = stereographic_r(s_L, x_j) ; 

		double Al = Al_n[jC] ;
		double Ze = Ze_n[jC] ;

		double r_Der_Ze = D1_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1], Ze_n[jC-1], dx) ;

		double t_Der_Ze = D1_backward_2ndOrder(Ze_n[jC], Ze_nm1[jC], Ze_nm2[jC], dt) ;

		eom_TR[jC] = 
		+ 	(2*t_Der_Ze*Ze)/r_j 
		- 	(2*r_Der_Ze*Al*pow(Ze,2))/r_j 
		- 	(Al*pow(Ze,3))/pow(r_j,2) 
		-	SE_LL_TR[jC] 
		;
	}
	return ;
}/*==========================================================================*/
static void compute_eom_ThTh(
	int Nx, 
	double dt, double dx,
	double s_L, int exc_jC, 
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* SE_LL_ThTh,
	double* eom_ThTh)
{
	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		double x_j = jC * dx ;
                double r_j = stereographic_r(s_L, x_j) ; 

		double Al = Al_n[jC] ;
		double Ze = Ze_n[jC] ;

		double r_Der_Al = D1_stereographic_center_2ndOrder(s_L, x_j, Al_n[jC+1], Al_n[jC-1], dx) ; 
		double r_Der_Ze = D1_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1], Ze_n[jC-1], dx) ;

		double rr_Der_Al = D2_stereographic_center_2ndOrder(s_L, x_j, Al_n[jC+1], Al_n[jC], Al_n[jC-1], dx) ;
		double rr_Der_Ze = D2_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1], Ze_n[jC], Ze_n[jC-1], dx) ;

		double t_Der_Al = D1_backward_2ndOrder(Al_n[jC], Al_nm1[jC], Al_nm2[jC], dt) ;
		double t_Der_Ze = D1_backward_2ndOrder(Ze_n[jC], Ze_nm1[jC], Ze_nm2[jC], dt) ;

		double tr_Der_Al = D1_backward_2ndOrder(
				D1_stereographic_center_2ndOrder(s_L, x_j, Al_n[jC+1],   Al_n[jC-1],   dx),
				D1_stereographic_center_2ndOrder(s_L, x_j, Al_nm1[jC+1], Al_nm1[jC-1], dx), 
				D1_stereographic_center_2ndOrder(s_L, x_j, Al_nm2[jC+1], Al_nm2[jC-1], dx),
			dt)
		;
		double tr_Der_Ze = D1_backward_2ndOrder(
				D1_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1],   Ze_n[jC-1],   dx),
				D1_stereographic_center_2ndOrder(s_L, x_j, Ze_nm1[jC+1], Ze_nm1[jC-1], dx), 
				D1_stereographic_center_2ndOrder(s_L, x_j, Ze_nm2[jC+1], Ze_nm2[jC-1], dx),
			dt)
		;
		eom_ThTh[jC] = 
		-	(pow(r_j,2)*pow(r_Der_Ze,2)) 
		+ 	t_Der_Ze*((pow(r_j,2)*r_Der_Al)/pow(Al,2) + r_j/Al) 
		+ 	(pow(r_j,2)*tr_Der_Ze)/Al 
		- 	pow(r_j,2)*rr_Der_Ze*Ze 
		- 	2*r_j*r_Der_Ze*Ze 
		- 	(pow(r_j,2)*r_Der_Al*t_Der_Al*Ze)/pow(Al,3) 
		+ 	(pow(r_j,2)*tr_Der_Al*Ze)/pow(Al,2) 
		- 	(pow(r_j,2)*rr_Der_Al*(-1 + pow(Ze,2)))/Al 
		+ 	r_Der_Al*((-3*pow(r_j,2)*r_Der_Ze*Ze)/Al + (r_j - r_j*pow(Ze,2))/Al) 
		- 	SE_LL_ThTh[jC] 
		;
		eom_ThTh[jC] /= pow(0.1+r_j,2) ; 
	}
	return ;
}
/*==========================================================================*/
static void compute_eom_scalar(
	int Nx, 
	double dt, double dx,
	double s_L, int exc_jC, 
	double* Al_n, 
	double* Ze_n,
	double* P_n, double* P_nm1, double* P_nm2,
	double* Q_n, 
	double* eom_scalar)
{
	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		double x_j = jC * dx ;
                double r_j = stereographic_r(s_L, x_j) ; 

		double Al = Al_n[jC] ;
		double Ze = Ze_n[jC] ;

		double P = P_n[jC] ;
		double Q = Q_n[jC] ;

		double r_Der_Al = D1_stereographic_center_2ndOrder(s_L, x_j, Al_n[jC+1], Al_n[jC-1], dx) ; 
		double r_Der_Ze = D1_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1], Ze_n[jC-1], dx) ;

		double r_Der_P = D1_stereographic_center_2ndOrder(s_L, x_j, P_n[jC+1], P_n[jC-1], dx) ; 
		double r_Der_Q = D1_stereographic_center_2ndOrder(s_L, x_j, Q_n[jC+1], Q_n[jC-1], dx) ;

		double t_Der_P = D1_backward_2ndOrder(P_n[jC], P_nm1[jC], P_nm2[jC], dt) ;

		eom_scalar[jC] = 
		-	t_Der_P/Al
		+	P*r_Der_Ze
		+	r_Der_Q
		+	Ze*r_Der_P
		+	(Q+P*Ze)*r_Der_Al/Al
		+	2*(Q/r_j+P*(Ze/r_j))
		;
		eom_scalar[jC] *= r_j ; 
	}
	return ;
}
/*==========================================================================*/
void compute_diagnostics_massless_scalar_GR(
	int Nx, int exc_jC, 
	double s_L,
	double dt, double dx,
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* P_n, double* P_nm1, double* P_nm2,
	double* Q_n, double* Q_nm1, double* Q_nm2,
	double* eom_TR,
	double* eom_ThTh,
	double* eom_scalar)
{
	double* SE_LL_TR   = calloc(Nx, sizeof(double)) ;
	double* SE_LL_ThTh = calloc(Nx, sizeof(double)) ;

	compute_SE_LL_components_massless_scalar(
		Nx, dx,
		s_L, exc_jC,
		Al_n, Ze_n,
		P_n,  Q_n,
		SE_LL_TR,
		SE_LL_ThTh)
	;
	compute_eom_TR(
		Nx, 
		dt, dx,
		s_L, exc_jC, 
		Al_n, 
		Ze_n, Ze_nm1, Ze_nm2,
		SE_LL_TR,
		eom_TR)
	;
	compute_eom_ThTh(
		Nx, 
		dt, dx,
		s_L, exc_jC, 
		Al_n, Al_nm1, Al_nm2,
		Ze_n, Ze_nm1, Ze_nm2,
		SE_LL_ThTh,
		eom_ThTh)
	;
	compute_eom_scalar(
		Nx, 
		dt, dx,
		s_L, exc_jC, 
		Al_n, 
		Ze_n,
		P_n, P_nm1, P_nm2,
		Q_n, 
		eom_scalar)
	;
	set_array_val(0, exc_jC, 0, P_n) ;
	set_array_val(0, exc_jC, 0, Q_n) ;
	set_array_val(0, exc_jC, 0, eom_TR) ;
	set_array_val(0, exc_jC, 0, eom_ThTh) ;
	set_array_val(0, exc_jC, 0, eom_scalar) ;

	free(SE_LL_TR) ;
	free(SE_LL_ThTh) ;

	return ;
}
