#include <math.h>
#include <string.h>
#include <assert.h>

#include "stencils.h"

/*
 * Dn: n = order of derivative (1: first, 2, second, etc derivative)
 * center: center difference about grid point n
 * forward: forward ' '
 * backward: backward ' '
 * CrankNicolson: about grid point n+1/2
 */
/*==========================================================================*/
extern double stereographic_r(double s_L, double x) ;

extern double stereographic_dr(double s_L, double x, double dx) ;
/*==========================================================================*/
/* 2nd order accurate stencils */
/*==========================================================================*/
extern double D1_center_2ndOrder(double vp1, double vm1, double dx) ;

extern double D2_center_2ndOrder(double vp1, double v0 , double vm1, double dx) ;

extern double D4_center_2ndOrder(double vp2, double vp1, double v0 , double vm1, double vm2, double dx) ;

extern double D1_forward_2ndOrder(double vp2, double vp1, double v0,  double dx) ;

extern double D1_backward_2ndOrder(double v0 , double vm1, double vm2, double dx) ;

extern double D2_forward_2ndOrder(double vp3, double vp2, double vp1, double v0,  double dx) ;

extern double D2_backward_2ndOrder(double v0, double vm1, double vm2, double vm3,  double dx) ;

extern double D1_CrankNicolson_2ndOrder(double vtp1, double vt, double dt) ;

extern double AVG_CrankNicolson_2ndOrder(double v1, double v2) ;
/*==========================================================================*/
/* 2nd order stencils
 * stereographic projection: r = x / (1-x/s_L) */
/*==========================================================================*/
extern double D1_stereographic_center_2ndOrder(double s_L, double x, double vp1, double vm1, double dx) ;

extern double D1_stereographic_forward_2ndOrder(double s_L, double x, double vp2, double vp1, double v0,  double dx) ;

extern double D1_stereographic_backward_2ndOrder(double s_L, double x, double v0 , double vm1, double vm2, double dx) ;

extern double D2_stereographic_center_2ndOrder(double s_L, double x, double vp1 , double v0, double vm1, double dx) ;

/*==========================================================================*/
/* see gr-qc/0302072 */
/*==========================================================================*/
void Kreiss_Oliger_filter(
	int Nx,
	int exc_jC,
	double* field)
{
	double epsilon_ko = 0.8 ;
	for (int iC=exc_jC+2; iC<Nx-2; iC++) {
		field[iC] -= (epsilon_ko/16.) * (
			field[iC+2] 
		+ 	(-4.*field[iC+1]) 
		+ 	(6.*field[iC]) 
		+ 	(-4.*field[iC-1]) 
		+ 	field[iC-2] 
		)
		;
	}
/* for outer excision boundary */
	epsilon_ko = 1.0 ;
	field[Nx-2] += (epsilon_ko/16.) * (
		field[Nx-1] 
	+ 	(-4.*field[Nx-2]) 
	+ 	(6.*field[Nx-3]) 
	+ 	(-4.*field[Nx-4]) 
	+ 	field[Nx-5] 
	) ;

	return ;
}
/*==========================================================================*/
void Kreiss_Oliger_left_inner_grid(
	double* field)
{
	double epsilon_ko = 1.0 ;
	field[1] += (epsilon_ko/16.) * (
			field[0] 
		+ 	(-4.*field[1]) 
		+ 	(6.*field[2]) 
		+ 	(-4.*field[3]) 
		+ 	field[4] 
	) ;

	return ;
}
/*==========================================================================*/
void Kreiss_Oliger_filter_origin(
	double *field,
	char *parity)
{
	double epsilon_ko = 0.9 ;
	if (strcmp(parity,"even")==0) {
		field[1] -= (epsilon_ko/16.) * (
			field[3] 
		+ 	(-4.*field[2]) 
		+ 	(6.*field[1]) 
		+ 	(-4.*field[0]) 
		+ 	field[1] 
		)
		;
		field[0] -= (epsilon_ko/16.) * (
			field[2] 
		+ 	(-4.*field[1]) 
		+ 	(6.*field[0]) 
		+ 	(-4.*field[1]) 
		+ 	field[2] 
		)
		;
	} else if (strcmp(parity,"odd")==0) {
		field[1] -= (epsilon_ko/16.) * (
			field[3] 
		+ 	(-4.*field[2]) 
		+ 	(6.*field[1]) 
		+ 	(-4.*field[0]) 
		+ 	(-field[1]) 
		)
		;	
		field[0] -= (epsilon_ko/16.) * (
			field[2] 
		+ 	(-4.*field[1]) 
		+ 	(6.*field[0]) 
		+ 	(-4.*(-field[1])) 
		+ 	(-field[2]) 
		)
		;
	} else if (strcmp(parity,"ignore_origin")==0) { 
		/* do nothing */
	} else {
		assert(strcmp(parity,"even")==0
		||	strcmp(parity,"odd")==0
		||	strcmp(parity,"ignore_origin")==0
		) ;
	}
	return ;
}
