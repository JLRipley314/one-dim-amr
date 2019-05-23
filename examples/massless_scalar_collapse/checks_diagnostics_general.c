#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "checks_diagnostics_general.h"
#include "stereographic_routines.h"

/*==========================================================================*/
/* general checks and diagnostics */
/*==========================================================================*/
static int compute_mass_aspect(
	int Nx,
	double s_L, double dx,
	double* Ze,
	double* mass_aspect)
{
	double r_j ;
	for (int jC=0; jC<Nx; jC++) {
		r_j = stereographic_r(s_L, dx*jC) ;
		mass_aspect[jC] = (r_j/2) * pow(Ze[jC],2) ;
	}
	mass_aspect[Nx-1] = mass_aspect[Nx-2] ;
	return 0 ;
} 
/*==========================================================================*/
/*==========================================================================*/
static int compute_ingoing_outgoing_null_characteristics(
	int Nx, double* Al, double* Ze,
	double* ingoing_null_characteristic,
	double* outgoing_null_characteristic)
{
	for (int jC=0; jC<Nx; jC++) {
		ingoing_null_characteristic[jC]  = - Al[jC] * (1 + Ze[jC]) ;
		outgoing_null_characteristic[jC] = + Al[jC] * (1 - Ze[jC]) ;
	}
	return 0 ;
} 
/*==========================================================================*/
static int compute_excision_point(
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
		&&  (inner_trapped > 0)
		) {
			outer_trapped = jC ;
		} 
	}
	if (inner_trapped>0 && outer_trapped>0) printf("inner trapped %d\touter trapped %d\n", inner_trapped, outer_trapped) ;
	if ((outer_trapped-inner_trapped) > buffer_size) {
		return outer_trapped-buffer_size ;
	}
	return 0 ;
}
/*==========================================================================*/
static int compute_Ricci_scalar(
	int Nx, int exc_jC, 
	double s_L,
	double dt, double dx,
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* Ricci_scalar)
{
	double x_j, r_j,
		Al, Ze, r_Der_Al, r_Der_Ze, rr_Der_Al, rr_Der_Ze,
		t_Der_Al, t_Der_Ze, tr_Der_Al, tr_Der_Ze 
	;
	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		x_j = jC * dx ;
                r_j = stereographic_r(s_L, x_j) ; 

		Al = Al_n[jC] ;
		Ze = Ze_n[jC] ;

		r_Der_Al = D1_stereographic_center_2ndOrder(x_j, s_L, Al_n[jC+1], Al_n[jC-1], dx) ; 
		r_Der_Ze = D1_stereographic_center_2ndOrder(x_j, s_L, Ze_n[jC+1], Ze_n[jC-1], dx) ;

		rr_Der_Al = D2_stereographic_center_2ndOrder(x_j, s_L, Al_n[jC+1], Al_n[jC], Al_n[jC-1], dx) ;
		rr_Der_Ze = D2_stereographic_center_2ndOrder(x_j, s_L, Ze_n[jC+1], Ze_n[jC], Ze_n[jC-1], dx) ;

		t_Der_Al = D1_backward_2ndOrder(Ze_n[jC], Al_nm1[jC], Al_nm2[jC], dt) ;
		t_Der_Ze = D1_backward_2ndOrder(Ze_n[jC], Ze_nm1[jC], Ze_nm2[jC], dt) ;

		tr_Der_Al = D1_backward_2ndOrder(
				D1_stereographic_center_2ndOrder(x_j, s_L, Al_n[jC+1],   Al_n[jC-1],   dx),
				D1_stereographic_center_2ndOrder(x_j, s_L, Al_nm1[jC+1], Al_nm1[jC-1], dx), 
				D1_stereographic_center_2ndOrder(x_j, s_L, Al_nm2[jC+1], Al_nm2[jC-1], dx),
			dt)
		;
		tr_Der_Ze = D1_backward_2ndOrder(
				D1_stereographic_center_2ndOrder(x_j, s_L, Ze_n[jC+1],   Ze_n[jC-1],   dx),
				D1_stereographic_center_2ndOrder(x_j, s_L, Ze_nm1[jC+1], Ze_nm1[jC-1], dx), 
				D1_stereographic_center_2ndOrder(x_j, s_L, Ze_nm2[jC+1], Ze_nm2[jC-1], dx),
			dt)
		;
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
	}
	Ricci_scalar[exc_jC] = Ricci_scalar[exc_jC+1] ; 
	Ricci_scalar[Nx-1]   = Ricci_scalar[Nx-2] ; 

	return 0 ;
}
/*==========================================================================*/
static int compute_Gauss_Bonnet_scalar(
	int Nx, int exc_jC, 
	double s_L,
	double dt, double dx,
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* Gauss_Bonnet_scalar)
{
	double x_j, r_j,
		Al, Ze, r_Der_Al, r_Der_Ze, rr_Der_Al, rr_Der_Ze,
		t_Der_Al, t_Der_Ze, tr_Der_Al, tr_Der_Ze 
	;
	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		x_j = jC * dx ;
                r_j = stereographic_r(s_L, x_j) ; 

		Al = Al_n[jC] ;
		Ze = Ze_n[jC] ;

		r_Der_Al = D1_stereographic_center_2ndOrder(x_j, s_L, Al_n[jC+1], Al_n[jC-1], dx) ; 
		r_Der_Ze = D1_stereographic_center_2ndOrder(x_j, s_L, Ze_n[jC+1], Ze_n[jC-1], dx) ;

		rr_Der_Al = D2_stereographic_center_2ndOrder(x_j, s_L, Al_n[jC+1], Al_n[jC], Al_n[jC-1], dx) ;
		rr_Der_Ze = D2_stereographic_center_2ndOrder(x_j, s_L, Ze_n[jC+1], Ze_n[jC], Ze_n[jC-1], dx) ;

		t_Der_Al = D1_backward_2ndOrder(Ze_n[jC], Al_nm1[jC], Al_nm2[jC], dt) ;
		t_Der_Ze = D1_backward_2ndOrder(Ze_n[jC], Ze_nm1[jC], Ze_nm2[jC], dt) ;

		tr_Der_Al = D1_backward_2ndOrder(
				D1_stereographic_center_2ndOrder(x_j, s_L, Al_n[jC+1],   Al_n[jC-1],   dx),
				D1_stereographic_center_2ndOrder(x_j, s_L, Al_nm1[jC+1], Al_nm1[jC-1], dx), 
				D1_stereographic_center_2ndOrder(x_j, s_L, Al_nm2[jC+1], Al_nm2[jC-1], dx),
			dt)
		;
		tr_Der_Ze = D1_backward_2ndOrder(
				D1_stereographic_center_2ndOrder(x_j, s_L, Ze_n[jC+1],   Ze_n[jC-1],   dx),
				D1_stereographic_center_2ndOrder(x_j, s_L, Ze_nm1[jC+1], Ze_nm1[jC-1], dx), 
				D1_stereographic_center_2ndOrder(x_j, s_L, Ze_nm2[jC+1], Ze_nm2[jC-1], dx),
			dt)
		;
		Gauss_Bonnet_scalar[jC] = 8*Ze*(
			r_Der_Al*t_Der_Al*pow(Ze,2) 
		+ 	pow(Al,3)*Ze*(3*pow(r_Der_Ze,2) + rr_Der_Ze*Ze) 
		- 	Al*Ze*(3*r_Der_Al*t_Der_Ze + tr_Der_Al*Ze) 
		+ 	pow(Al,2)*(
			-	2*r_Der_Ze*t_Der_Ze 
			- 	(rr_Der_Al + tr_Der_Ze)*Ze 
			+ 	rr_Der_Al*pow(Ze,3) 
			+ 	r_Der_Al*r_Der_Ze*(-2 + 5*pow(Ze,2))
			)
		)/(pow(r_j,2)*pow(Al,3))
		;
	}
	Gauss_Bonnet_scalar[exc_jC] = Gauss_Bonnet_scalar[exc_jC+1] ; 
	Gauss_Bonnet_scalar[Nx-1]   = Gauss_Bonnet_scalar[Nx-2] ; 
	return 0 ;
}
/*==========================================================================*/
void compute_checks_diagnostics_general(
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
{
	int buffer_size = 10 ; /* TO DO: set buffer size to some percentage of ADM mass */
	compute_excision_point(
		Nx, exc_jC, buffer_size, outgoing_null_characteristic)
	;
	compute_mass_aspect(
		Nx,
		s_L, dx,
		Ze_n,
		mass_aspect)
	;
	compute_ingoing_outgoing_null_characteristics(
		Nx, Al_n, Ze_n,
		ingoing_null_characteristic,
		outgoing_null_characteristic)
	;
	compute_Ricci_scalar(
		Nx, exc_jC, 
		s_L,
		dt, dx,
		Al_n, Al_nm1, Al_nm2,
		Ze_n, Ze_nm1, Ze_nm2,
		Ricci_scalar)
	;
	compute_Gauss_Bonnet_scalar(
		Nx, exc_jC, 
		s_L,
		dt, dx,
		Al_n, Al_nm1, Al_nm2,
		Ze_n, Ze_nm1, Ze_nm2,
		Gauss_Bonnet_scalar)
	;
	return ;
}
