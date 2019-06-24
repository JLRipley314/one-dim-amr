#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include "evolution_routines_GR.h"
#include "stencils.h"

static const double machine_epsilon = 1e-14 ;

/*==========================================================================*/
/* 	Notations:
	s_L: stereographic length
	exc_jC: excision grid point (index jC)
*/
/*==========================================================================*/

/*==========================================================================*/
static inline double weighted_infty_norm(double weight, double val_1, double val_2)
{
        return (fabs(val_2)>fabs(val_1)) ? fabs(weight*val_2) : fabs(weight*val_1) ;
}
/*==========================================================================*/
/* ODE solvers for lapse and shift */
/*==========================================================================*/
/* there is Jr/Ze term that is zero in empty
 * space but undefined if below machine precision.
 * The slope is zero though in empty spacce so
 * to that precision we do this. */
/*==========================================================================*/
static double compute_iteration_Al(
	double s_L,
	int size,
	double dx,
	int start_jC,
	double x_lower,
	double* Al, 	double* Ze, 
	double*  P, 	double*  Q)
{
        double res_infty_norm = 0 ; /* returning this */
 /* scalar field functions */
        for (int jC=start_jC; jC<size-1; jC++) {
                double x_joh = x_lower + ((jC+1) + jC) * dx / 2 ;

                double r_joh = stereographic_r(s_L, x_joh) ;

                double dr = stereographic_dr(s_L, x_joh, dx) ;
                double Al_joh = (Al[jC+1] + Al[jC]) / 2 ;
                double Ze_joh = (Ze[jC+1] + Ze[jC]) / 2 ;

                double Q_joh = (Q[jC+1] + Q[jC]) / 2 ;
                double P_joh = (P[jC+1] + P[jC]) / 2 ;

                double r_Der_Al_joh = (Al[jC+1] - Al[jC]) / dr ;

                double Jr_joh = - Q_joh * P_joh ;

                if ((fabs(Jr_joh) < 10*machine_epsilon)
                &&  (fabs(Ze_joh) < 10*machine_epsilon)
                ) {
                        Al[jC+1] = Al[jC] ;
                } else {
                        double res_Al =
                        +       r_Der_Al_joh*Ze_joh
                        -       (r_joh*Al_joh*Jr_joh)/2.
                        ;
                        double jac_Al =
                                Ze_joh/dr
                        -       r_joh*Jr_joh/4.
                        ;
                        Al[jC+1] -= res_Al / jac_Al ;
			if ((isnan(res_Al) != 0)
			||  (isnan(jac_Al) != 0)
			) {
				printf("jC:%d\tres_Al:%.e\tjac_Al:%.e\n", jC, res_Al, jac_Al) ;
				exit(EXIT_FAILURE) ;
			}
			res_infty_norm = weighted_infty_norm(1-x_joh/s_L, res_Al, res_infty_norm) ;
                }
        }
        return res_infty_norm ;
}
/*==========================================================================*/
static double compute_iteration_Ze(
	double s_L,
	int size,
	double dx,
	int start_jC,
	double x_lower,
	double* Al, 	double* Ze, 
	double*  P, 	double*  Q)
{
        double res_infty_norm = 0 ; /* returning this */

        for (int jC=start_jC; jC<size-1; jC++) {
                double x_joh = x_lower + ((jC+1) + jC) * dx / 2 ; 

                double x_jp1 = (jC+1) * dx ;
                double x_j   = (jC+0) * dx ;

                double r_joh = stereographic_r(s_L, x_joh) ;
                double r_jp1 = stereographic_r(s_L, x_jp1) ;
                double r_j   = stereographic_r(s_L, x_j  ) ; 

                double dr = stereographic_dr(s_L, x_joh, dx) ;

                double Al_joh = (Al[jC+1] + Al[jC]) / 2. ;

                double Al_sqrd_jp1 = pow(Al[jC+1], 2) ;
                double Al_sqrd_j   = pow(Al[jC+0], 2) ;

                double Ze_sqrd_jp1 = pow(Ze[jC+1], 2) ;
                double Ze_sqrd_j   = pow(Ze[jC+0], 2) ;

                double Q_joh = (Q[jC+1] + Q[jC]) / 2. ;
                double P_joh = (P[jC+1] + P[jC]) / 2. ;

                double rho_joh = (1./2) * (pow(Q_joh,2) + pow(P_joh,2)) ;
/*---------------------------------------------------------------------------*/
                double res_Ze_sqrd = 
                        (       (r_jp1)*Al_sqrd_jp1*Ze_sqrd_jp1
                        -       (r_j  )*Al_sqrd_j  *Ze_sqrd_j
                        )/dr
                -       pow(r_joh,2)*pow(Al_joh,2)*rho_joh
                ;
                double jac_Ze_sqrd = (r_jp1)*Al_sqrd_jp1/dr
                ;
                Ze_sqrd_jp1 -= res_Ze_sqrd/jac_Ze_sqrd ;
                Ze[jC+1]  = sqrt(Ze_sqrd_jp1) ;
/*---------------------------------------------------------------------------*/
                if ((isnan(res_Ze_sqrd) != 0)    
                ||  (isnan(jac_Ze_sqrd) != 0)
                ) {    
                        printf("jC:%d\tZe:%.6e\tres_Ze_sqrd:%.e\tjac_Ze_sqrd:%.e\n", jC, Ze[jC+1], res_Ze_sqrd, jac_Ze_sqrd) ;
                        exit(EXIT_FAILURE) ;
                }
                res_infty_norm = weighted_infty_norm(1-x_joh/s_L, res_Ze_sqrd, res_infty_norm) ;
        }
        return res_infty_norm ;
}
/*==========================================================================*/
static double compute_iteration_excision_boundary_condition_Ze(
	double s_L,
	double dt, 	double dx,
	int exc_jC,
	double x_lower,
	double* Al_n,  double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double* P_n,   double* P_nm1,  double* Q_n,  double* Q_nm1)
{
        double x_j = x_lower + (dx * exc_jC) ;
        double r_j = stereographic_r( s_L, x_j) ;
        double dr  = stereographic_dr(s_L, x_j, dx) ;

        double Al = (Al_n[exc_jC] + Al_nm1[exc_jC]) / 2. ;
        double Ze = (Ze_n[exc_jC] + Ze_nm1[exc_jC]) / 2. ;
        double P  = (P_n[exc_jC]  + P_nm1[exc_jC])  / 2. ;
        double Q  = (Q_n[exc_jC]  + Q_nm1[exc_jC])  / 2. ;

        double t_Der_Ze = D1_CrankNicolson_2ndOrder(Ze_n[exc_jC], Ze_nm1[exc_jC], dt) ;

        double	r_Der_Ze  = D1_forward_2ndOrder(Ze_n[exc_jC+2],   Ze_n[exc_jC+1],   Ze_n[exc_jC],   dr) ;
		r_Der_Ze += D1_forward_2ndOrder(Ze_nm1[exc_jC+2], Ze_nm1[exc_jC+1], Ze_nm1[exc_jC], dr) ;
		r_Der_Ze /= 2 ;

        double SE_LL_TR        = (Al*(2*P*Q + (pow(P,2) + pow(Q,2))*Ze))/2. ;
        double Ze_Der_SE_LL_TR = (Al*(0     + (pow(P,2) + pow(Q,2))*1 ))/2. ;

        double res_Ze = 
		t_Der_Ze
	-       (r_j*SE_LL_TR)/(2.*Ze)
	-       r_Der_Ze*Al*Ze
	-       (Al*pow(Ze,2))/(2.*r_j)

        ;
        double jac_Ze = 
		1/dt
	-       (r_Der_Ze*Al)/2.
	+       (
			pow(r_j,2)*SE_LL_TR
		-       pow(r_j,2)*Ze_Der_SE_LL_TR*Ze
		-       2*Al*pow(Ze,3)
	)/(4.*r_j*pow(Ze,2))
        ;

        Ze_n[exc_jC] -= res_Ze/jac_Ze ;

        return fabs(res_Ze) ;
}
/*==========================================================================*/
void solve_Al_Ze_GR(
	double s_L,
	int Nx,
	double dt, 	double dx,
	double err_tolerance,
	bool excision_on,
	int start_jC,
	double bbox[2],
	double* Al_n, 	double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double*  P_n, 	double*  P_nm1, double*  Q_n, double*  Q_nm1)
{
	if (start_jC==Nx-1) { /* if interior to the excision point */
		return ;
	}
	double res = 0 ;
	double x_lower = bbox[0] ;
/* to avoid problems with r=infty when x=s_L */	
	int size = Nx ;
	if (fabs(bbox[1]-s_L)<machine_epsilon) size = Nx-1 ; 
	do {
		res = 0 ;
		if ((excision_on==true)
		&&  (start_jC>0)
		) {	
			res += compute_iteration_excision_boundary_condition_Ze(
				s_L,
				dt, 	dx,
				start_jC,
				x_lower,
				Al_n,  Al_nm1, Ze_n, Ze_nm1,
				 P_n,   P_nm1,  Q_n,  Q_nm1)
			;
		}
		res += compute_iteration_Ze(
			s_L,
			size,
			dx,
			start_jC,
			x_lower,
			Al_n, 	Ze_n, 
			P_n, 	Q_n)
		;
		res += compute_iteration_Al(
			s_L,
			size,
			dx,
			start_jC,
			x_lower,
			Al_n, 	Ze_n, 
			P_n, 	Q_n)
		;
		if (size==Nx-1) {
			Al_n[Nx-1] = Al_n[Nx-2] ;
			Ze_n[Nx-1] = 0 ;
		}
	} while (res>err_tolerance) ;
	return ;
}
/*==========================================================================*/
static double compute_iteration_Crank_Nicolson_PQ(
	double s_L,
	int Nx,
	double dt, 	double dx,
	bool excision_on,
	int exc_jC,
	double bbox[2],
	bool perim_interior[2],
	double* Al_n, 	double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double*  P_n, 	double*  P_nm1, double*  Q_n, double*  Q_nm1)
{
	double lower_x = bbox[0] ;

	int size = Nx ;
	if (fabs(bbox[1]-s_L)<machine_epsilon) size = Nx-1 ; /* to avoid problems with r=infty when x=s_L */
	double res_infty_norm = 0 ; /* returning this */
/*--------------------------------------------------------------------------*/
/* interior: we go to Nx-2 as we do not want to actually include the point
   at infinity in our computational domain */
/*--------------------------------------------------------------------------*/
	for (int jC=exc_jC+1;jC<size-1;jC++) {
		double x_j   = lower_x + (dx * (jC)  ) ;
		double x_jp1 = lower_x + (dx * (jC+1)) ;
		double x_jm1 = lower_x + (dx * (jC-1)) ;

		double r_jp1 = stereographic_r(s_L, x_jp1) ;
		double r_jm1 = stereographic_r(s_L, x_jm1) ;

		double dr = stereographic_dr(s_L, x_j, dx) ;

		double Al_jp1 = (Al_n[jC+1]+Al_nm1[jC+1])/2 ;
		double Al_jm1 = (Al_n[jC-1]+Al_nm1[jC-1])/2 ;

		double Ze_jp1 = (Ze_n[jC+1]+Ze_nm1[jC+1])/2 ;
		double Ze_jm1 = (Ze_n[jC-1]+Ze_nm1[jC-1])/2 ;

		double P_jp1 = (P_n[jC+1]+P_nm1[jC+1])/2 ;
		double P_jm1 = (P_n[jC-1]+P_nm1[jC-1])/2 ;

		double Q_jp1 = (Q_n[jC+1]+Q_nm1[jC+1])/2 ;
		double Q_jm1 = (Q_n[jC-1]+Q_nm1[jC-1])/2 ;
	/* Q field 
	*/
		double res_Q = D1_CrankNicolson_2ndOrder(
				Q_n[jC], 
				Q_nm1[jC], 
				dt)
		-	D1_center_2ndOrder(
				Al_jp1*(P_jp1 + Ze_jp1*Q_jp1),
				Al_jm1*(P_jm1 + Ze_jm1*Q_jm1),
			dr)
		;
		double jac_Q = (1./dt) ;
	/* P field 
	*/
		double res_P = D1_CrankNicolson_2ndOrder(
				P_n[jC], 
				P_nm1[jC], 
				dt)
		-	3*(	pow(r_jp1,2)*Al_jp1*(Q_jp1 + Ze_jp1*P_jp1)
			-	pow(r_jm1,2)*Al_jm1*(Q_jm1 + Ze_jm1*P_jm1)
			)/(pow(r_jp1,3)-pow(r_jm1,3))
		;
		double jac_P = (1./dt) ;
	/* one iteration 
	*/
		Q_n[jC] -= res_Q / jac_Q ;
		P_n[jC] -= res_P / jac_P ;

		res_infty_norm = weighted_infty_norm(1-x_j/s_L, res_Q, res_infty_norm) ;
		res_infty_norm = weighted_infty_norm(1-x_j/s_L, res_P, res_infty_norm) ;
	}
/*--------------------------------------------------------------------------*/
/* lower */
/* Q[0] = 0, so do not change anything. r_Der_P=0 at r=0 */	
/*--------------------------------------------------------------------------*/
	if ((exc_jC == 0) 
	&&  (perim_interior[0] == false)
	) {
		double x_j = 0 ;
		double dr  = pow(1. - (x_j/s_L), -2) * dx;

		double res_P = D1_forward_2ndOrder(
				P_n[2], P_n[1], P_n[0], 
				dr) ;
		double jac_P = - (3./2.) / dr ;

		P_n[0] -= res_P / jac_P ;

		res_infty_norm = weighted_infty_norm(1-x_j/s_L, res_P, res_infty_norm) ;
	}
	if ((exc_jC > 0) 
	&&  (excision_on==true)
	) {
		double x_j   = lower_x + (dx * (exc_jC))   ;
		double x_jp1 = lower_x + (dx * (exc_jC+1)) ;
		double x_jp2 = lower_x + (dx * (exc_jC+2)) ;

		double r_j   = stereographic_r(s_L, x_j  ) ;
		double r_jp1 = stereographic_r(s_L, x_jp1) ;
		double r_jp2 = stereographic_r(s_L, x_jp2) ;

		double dr = stereographic_dr(s_L, x_j, dx) ;

		double Al_j   = (Al_n[exc_jC+0]+Al_nm1[exc_jC+0])/2 ;
		double Al_jp1 = (Al_n[exc_jC+1]+Al_nm1[exc_jC+1])/2 ;
		double Al_jp2 = (Al_n[exc_jC+2]+Al_nm1[exc_jC+2])/2 ;

		double Ze_j   = (Ze_n[exc_jC+0]+Ze_nm1[exc_jC+0])/2 ;
		double Ze_jp1 = (Ze_n[exc_jC+1]+Ze_nm1[exc_jC+1])/2 ;
		double Ze_jp2 = (Ze_n[exc_jC+2]+Ze_nm1[exc_jC+2])/2 ;

		double P_j   = (P_n[exc_jC+0]+P_nm1[exc_jC+0])/2 ;
		double P_jp1 = (P_n[exc_jC+1]+P_nm1[exc_jC+1])/2 ;
		double P_jp2 = (P_n[exc_jC+2]+P_nm1[exc_jC+2])/2 ;

		double Q_j   = (Q_n[exc_jC+0]+Q_nm1[exc_jC+0])/2 ;
		double Q_jp1 = (Q_n[exc_jC+1]+Q_nm1[exc_jC+1])/2 ;
		double Q_jp2 = (Q_n[exc_jC+2]+Q_nm1[exc_jC+2])/2 ;
	/* Q field 
	*/
		double res_Q = D1_CrankNicolson_2ndOrder(
				Q_n[exc_jC], 
				Q_nm1[exc_jC], 
				dt)
		-	D1_forward_2ndOrder(
				Al_jp2*(P_jp2 + Ze_jp2*Q_jp2),
				Al_jp1*(P_jp1 + Ze_jp1*Q_jp1),
				Al_j  *(P_j   + Ze_j*  Q_j  ),
				dr)
		;
		double jac_Q =	1/dt 
		- 	(1./2.)*(-3./(2*dr))*Al_j*Ze_j 
		;
	/* P field 
	*/
		double res_P = D1_CrankNicolson_2ndOrder(
				P_n[exc_jC], 
				P_nm1[exc_jC], 
				dt)
		-	 pow(r_j,-2)*D1_forward_2ndOrder(
				pow(r_jp2,2)*Al_jp2*(Q_jp2 + Ze_jp2*P_jp2),
				pow(r_jp1,2)*Al_jp1*(Q_jp1 + Ze_jp1*P_jp1),
				pow(r_j,  2)*Al_j  *(Q_j   + Ze_j  *P_j  ),
				dr)
		;
		double jac_P = 1./dt
		-	(1./2.)*(-3./(2*dr))*Al_j*Ze_j
		;
	/****/
		Q_n[exc_jC] -= res_Q / jac_Q ;
		P_n[exc_jC] -= res_P / jac_P ;

		res_infty_norm = weighted_infty_norm(1-x_j/s_L, res_Q, res_infty_norm) ;
		res_infty_norm = weighted_infty_norm(1-x_j/s_L, res_P, res_infty_norm) ;
	}
/*--------------------------------------------------------------------------*/
/* dirichlet outer boundary conditions if outermost level;
 * otherwise do not evolve outer boundary anyways (it is interpolated) */
/*--------------------------------------------------------------------------*/
	return res_infty_norm ;
} 
/*===========================================================================*/
void advance_tStep_PQ_massless_scalar_GR(
	double s_L,
	int Nx, 
	double dt, double dx, 
	double err_tolerance,
	bool excision_on,
	int exc_jC,
	double bbox[2], 
	bool perim_interior[2],
	double* Al_n, double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double*  P_n, double*  P_nm1, double*  Q_n, double*  Q_nm1)
{ 
	if (exc_jC==Nx-1) {
		return ;
	}
	double res = 0 ;
	do {
		res = compute_iteration_Crank_Nicolson_PQ(
			s_L,
			Nx,
			dt, 	dx,
			excision_on,
			exc_jC,
			bbox,
			perim_interior,
			Al_n, 	Al_nm1, Ze_n, Ze_nm1,
			 P_n, 	 P_nm1,  Q_n,  Q_nm1)
		;
	} while (res>err_tolerance) ;

	Kreiss_Oliger_filter(Nx, exc_jC, P_n) ;
	Kreiss_Oliger_filter(Nx, exc_jC, Q_n) ;

	if ((fabs(bbox[0])<machine_epsilon)
	&&  (exc_jC>0)
	) {
		Kreiss_Oliger_filter_origin(P_n, "even") ;
		Kreiss_Oliger_filter_origin(Q_n, "odd") ;
	}
	return ;
}	
