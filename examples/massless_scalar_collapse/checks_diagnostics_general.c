#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "checks_diagnostics_general.h"
#include "stereographic_routines.h"

/*==========================================================================*/
/* general checks and diagnostics */
/*==========================================================================*/
int compute_ingoing_outgoing_null_characteristics(
	int Nx, double* Al, double* Ze,
	double* ingoing_characteristic,
	double* outgoing_characteristic)
{
	for (int jC=0; jC<Nx; jC++) {
		ingoing_characteristic[jC]  = - Al[jC] * (1 + Ze[jC]) ;
		outgoing_characteristic[jC] = + Al[jC] * (1 - Ze[jC]) ;
	}
	return 0 ;
} 
/*==========================================================================*/
int compute_excision_point(
	int Nx, int exc_jC, int buffer_size, double* outgoing_null_characteristic)
{
	int inner_trapped = 0 ;
	int outer_trapped = 0 ;
	for (int jC=exc_jC; jC<Nx; jC++) {
		if ((outgoing_null_characteristic[jC] < 0)
		&&  (inner_trapped = 0)
		) {
			inner_trapped = jC ;
		} 
		if ((outgoing_null_characteristic[jC] < 0)
		&&  (outer_trapped > 0)
		) {
			outer_trapped = jC ;
		} 
	}
	if ((outer_trapped-inner_trapped) > buffer_size) {
		return outer_trapped-buffer_size ;
	}
	return 0 ;
}
/*==========================================================================*/
int compute_Ricci_scalar(
	int Nx, int exc_jC, 
	double s_L,
	double dx,
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* Ricci_scalar)
{
	for (int jC=exc_jC; jC<Nx; jC++) {
/*		x_j = jC * dx ;

		Al = Al_n[jC] ;
		Ze = Ze_n[jC] ;

		x_Der_Al = D1_center_2ndOrder(Al_n[jC+1], Al_n[jC-1], dx) ;
		x_Der_Ze = D1_center_2ndOrder(Ze_n[jC+1], Ze_n[jC-1], dx) ;

		xx_Der_Al = D2_center_2ndOrder(Al_n[jC+1], Al_n[jC], Al_n[jC-1], dx) ;
		xx_Der_Ze = D2_center_2ndOrder(Ze_n[jC+1], Ze_n[jC], Ze_n[jC-1], dx) ;

		r_Der_Al = compute_r_Der(s_L, x_j, x_Der_Al) ;
		r_Der_Ze = compute_r_Der(s_L, x_j, x_Der_Ze) ;
		

		Ricci_scalar[jC] = 2*(
			pow(r_j,2)*r_Der_Al*t_Der_Al*Ze 
		- 	pow(r_j,2)*Al*(r_Der_Al*t_Der_Ze + tr_Der_Al*Ze) 
		+ 	pow(Al,3)*(pow(r_j,2)*pow(r_Der_Ze,2) + r_j*(r_j*rr_Der_Ze + 4*r_Der_Ze)*Ze + pow(Ze,2)) 
		+ 	r_j*pow(Al,2)*(
			-	(r_j*tr_Der_Ze) 
			- 	2*t_Der_Ze 
			+ 	r_j*rr_Der_Al*(-1 + pow(Ze,2)) 
			+ 	r_Der_Al*(-2 + 3*r_j*r_Der_Ze*Ze + 2*pow(Ze,2))
			)
		)/(pow(r_j,2)*pow(Al,3))
		;
*/	}
	return 0 ;
}
/*==========================================================================*/
 
