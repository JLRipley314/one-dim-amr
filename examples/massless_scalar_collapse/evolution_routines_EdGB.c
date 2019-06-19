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
	double s_L,	double c_gbs,
	int Nx,
	double dt, 	double dx,
	int start_jC,
	double bbox[2],
	double *Al, 	double *Ze, 
	double * P, 	double * Q)
{
	double phiphi_Der_f_joh = 1 ;
	double phi_Der_f_joh = 1 ;

	int size = 0 ;
	if (fabs(bbox[1]-s_L)<machine_epsilon) size = Nx-1 ; /* to avoid problems with r=infty when x=s_L */
	else size = Nx ;
        double res_infty_norm = 0 ; /* returning this */
 /* scalar field functions */
        for (int jC=start_jC; jC<Nx-1; jC++) {
                double x_joh = ((jC+1) + jC) * dx / 2 ;

                double r_joh = stereographic_r(s_L, x_joh) ;

                double dr = stereographic_dr(s_L, x_joh, dx) ;
                double Al_joh = (Al[jC+1] + Al[jC]) / 2 ;
                double Ze_joh = (Ze[jC+1] + Ze[jC]) / 2 ;

                double Q_joh = (Q[jC+1] + Q[jC]) / 2 ;
                double P_joh = (P[jC+1] + P[jC]) / 2 ;

                double r_Der_Al_joh = (Al[jC+1] - Al[jC]) / dr ;
                double r_Der_Ze_joh = (Ze[jC+1] - Ze[jC]) / dr ;

                double r_Der_P_joh = (P[jC+1] - P[jC]) / dr ;

                double Jr_joh = - Q_joh * P_joh ;
		double Al_Der_Jr_joh = 0 ;

                if ((fabs(Jr_joh) < 10*machine_epsilon)
                &&  (fabs(Ze_joh) < 10*machine_epsilon)
                ) {
                        Al[jC+1] = Al[jC] ;
                } else {
                        double res_Al = 			
			-	(r_joh*Al_joh*Jr_joh)/(2.*Ze_joh) 
			+ 	(4*phi_Der_f_joh*r_Der_P_joh*c_gbs*Al_joh*Ze_joh)/r_joh 
			+ 	(4*phi_Der_f_joh*r_Der_Ze_joh*c_gbs*Al_joh*Q_joh*Ze_joh)/r_joh 
			+ 	(4*phiphi_Der_f_joh*c_gbs*Al_joh*P_joh*Q_joh*Ze_joh)/r_joh 
			+ 	(
					1 
				- 	(8*phi_Der_f_joh*c_gbs*Q_joh)/r_joh 
				- 	(8*phi_Der_f_joh*c_gbs*P_joh*Ze_joh)/r_joh 
				+ 	(4*phi_Der_f_joh*c_gbs*Q_joh*pow(Ze_joh,2))/r_joh
			)*r_Der_Al_joh
			;
                        double jac_Al = 
			-	(Al_Der_Jr_joh*r_joh*Al_joh)/(4.*Ze_joh) 
			- 	(r_joh*Jr_joh)/(4.*Ze_joh) 
			+ 	(2*phi_Der_f_joh*r_Der_P_joh*c_gbs*Ze_joh)/r_joh 
			+ 	(2*phi_Der_f_joh*r_Der_Ze_joh*c_gbs*Q_joh*Ze_joh)/r_joh 
			+ 	(2*phiphi_Der_f_joh*c_gbs*P_joh*Q_joh*Ze_joh)/r_joh 
			+ 	(
					1 
				- 	(8*phi_Der_f_joh*c_gbs*Q_joh)/r_joh 
				- 	(8*phi_Der_f_joh*c_gbs*P_joh*Ze_joh)/r_joh 
				+ 	(4*phi_Der_f_joh*c_gbs*Q_joh*pow(Ze_joh,2))/r_joh
			)/dr
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
	if (size==Nx-1) Al[Nx-1] = Al[Nx-2] ;

        return res_infty_norm ;
}
/*==========================================================================*/
static double compute_iteration_Ze(
	double s_L,	double c_gbs,
	int Nx,
	double dt, 	double dx,
	int start_jC,
	double bbox[2],
	double *Al, 	double *Ze, 
	double * P, 	double * Q)
{
	double phi_Der_f_jp1 = 1 ;
	double phi_Der_f_joh = 1 ;
	double phi_Der_f_j   = 1 ;

	int size = 0 ;
	if (fabs(bbox[1]-s_L)<machine_epsilon) size = Nx-1 ; /* to avoid problems with r=infty when x=s_L */
	else size = Nx ;
 /* scalar field functions */   

        double res_infty_norm = 0 ; /* returning this */

        for (int jC=start_jC; jC<size-1; jC++) {
                double x_joh = ((jC+1) + jC) * dx / 2 ; 

                double x_jp1 = (jC+1) * dx ;
                double x_j   = (jC+0) * dx ;

                double r_joh = stereographic_r(s_L, x_joh) ;
                double r_jp1 = stereographic_r(s_L, x_jp1) ;
                double r_j   = stereographic_r(s_L, x_j  ) ; 

                double dr = stereographic_dr(s_L, x_joh, dx) ;

                double Al_joh = (Al[jC+1] + Al[jC]) / 2. ;

                double Q_jp1 = Q[jC+1] ;
                double Q_j   = Q[jC] ;

                double Al_sqrd_jp1 = pow(Al[jC+1], 2) ;
                double Al_sqrd_j   = pow(Al[jC+0], 2) ;
                double Ze_sqrd_jp1 = pow(Ze[jC+1], 2) ;
                double Ze_sqrd_j   = pow(Ze[jC+0], 2) ;

                double Q_joh = (Q[jC+1] + Q[jC]) / 2. ;
                double P_joh = (P[jC+1] + P[jC]) / 2. ;

                double rho_joh = (1./2) * (pow(Q_joh,2) + pow(P_joh,2)) ;
/*---------------------------------------------------------------------------*/
		double res_Ze_sqrd = 
		( 
			(	(r_jp1 - 8*c_gbs*phi_Der_f_jp1*Q_jp1)*Al_sqrd_jp1*Ze_sqrd_jp1
			-	(r_j   - 8*c_gbs*phi_Der_f_j  *Q_j  )*Al_sqrd_j  *Ze_sqrd_j
			)/dr
		)
		- (	
			(8*c_gbs*phi_Der_f_joh*P_joh/Al_joh)*(
				pow(Al_sqrd_jp1*Ze_sqrd_jp1,3./2.) 
			- 	pow(Al_sqrd_j  *Ze_sqrd_j  ,3./2.)
			)/dr
		)	
		-	pow(r_joh,2)*pow(Al_joh,2)*rho_joh
		;
		double jac_Ze_sqrd = 
		(
			(
				(r_jp1 - 8*c_gbs*phi_Der_f_jp1*Q_jp1)*Al_sqrd_jp1
			)/dr
		)
		-(	
			(8*c_gbs*phi_Der_f_joh*P_joh/Al_joh)*(
				(3./2.)*Al_sqrd_jp1*pow(Al_sqrd_jp1*Ze_sqrd_jp1,1./2.)
			)/dr
		)	
		;
                Ze_sqrd_jp1 -= res_Ze_sqrd/jac_Ze_sqrd ;
                Ze[jC+1]  = sqrt(Ze_sqrd_jp1) ;
/*---------------------------------------------------------------------------*/	
                if ((isnan(res_Ze_sqrd) != 0)    
                ||  (isnan(jac_Ze_sqrd) != 0)
                ) {    
                        printf("jC:%d\tZe_j:%.6e\tZe_jp1:%.6e\tres_Ze_sqrd:%.e\tjac_Ze_sqrd:%.e\n", jC, Ze[jC], Ze[jC+1], res_Ze_sqrd, jac_Ze_sqrd) ;
                        exit(EXIT_FAILURE) ;
                }
                res_infty_norm = weighted_infty_norm(1-x_joh/s_L, res_Ze_sqrd, res_infty_norm) ;
        }
	if (size==Nx-1) Ze[Nx-1] = 0 ;

        return res_infty_norm ;
}
/*==========================================================================*/
static double free_evolution_Ze_res(
	double c_gbs, 
	double r_j, 
	double Al, double Ze, double Q, double P,
	double phi_Der_f, double phiphi_Der_f,
	double SE_LL_TR, double SE_LL_ThTh,
	double r_Der_Al, double r_Der_Ze, double r_Der_Q, double r_Der_P,
	double t_Der_Ze)
{
	return 
		4*phi_Der_f*c_gbs*P*SE_LL_TR - (r_j*SE_LL_TR)/(2.*Ze) + (4*phi_Der_f*c_gbs*Q*SE_LL_TR)/Ze - (4*phiphi_Der_f*c_gbs*Al*P*Q*Ze)/r_j + (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*P*pow(Q,2)*Ze)/pow(r_j,2) - (Al*pow(Ze,2))/(2.*r_j) - (4*phiphi_Der_f*c_gbs*Al*pow(P,2)*pow(Ze,2))/r_j - (4*phi_Der_f*c_gbs*Al*Q*pow(Ze,2))/pow(r_j,2) + (64*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*pow(P,2)*Q*pow(Ze,2))/pow(r_j,2) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,2)*pow(Ze,2))/pow(r_j,3) - (4*phi_Der_f*c_gbs*Al*P*pow(Ze,3))/pow(r_j,2) + (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*pow(P,3)*pow(Ze,3))/pow(r_j,2) + (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*Q*pow(Ze,3))/pow(r_j,3) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,4))/pow(r_j,3) + r_Der_P*((-4*phi_Der_f*c_gbs*Al*Ze)/r_j + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*Ze)/pow(r_j,2) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Ze,2))/pow(r_j,2)) + r_Der_Q*((-4*phi_Der_f*c_gbs*Al*pow(Ze,2))/r_j + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,2) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,2)) + pow(r_Der_Ze,2)*((-128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,4))/pow(r_j,3) + (768*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*Q*pow(Ze,4))/pow(r_j,4) + (768*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*P*pow(Ze,5))/pow(r_j,4)) + r_Der_Ze*(-(Al*Ze) + (12*phi_Der_f*c_gbs*Al*Q*Ze)/r_j - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,2)*Ze)/pow(r_j,2) + (12*phi_Der_f*c_gbs*Al*P*pow(Ze,2))/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*Q*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*SE_LL_TR*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,3))/pow(r_j,2) - (256*pow(phi_Der_f,3)*r_Der_P*pow(c_gbs,3)*Al*pow(Ze,4))/pow(r_j,4) - (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*P*Q*pow(Ze,4))/pow(r_j,4) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,5))/pow(r_j,4) - (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*Al*pow(Ze,5))/pow(r_j,4) - (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,5))/pow(r_j,4)) + pow(r_Der_Al,2)*((256*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,5))/(pow(r_j,4)*Al) + (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,6))/(pow(r_j,4)*Al)) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*SE_LL_ThTh*pow(Ze,4))/(pow(r_j,4)*(r_j - 8*phi_Der_f*c_gbs*(Q + P*Ze))) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*Q*SE_LL_ThTh*pow(Ze,4))/(pow(r_j,5)*(r_j - 8*phi_Der_f*c_gbs*(Q + P*Ze))) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*P*SE_LL_ThTh*pow(Ze,5))/(pow(r_j,5)*(r_j - 8*phi_Der_f*c_gbs*(Q + P*Ze))) + t_Der_Ze*(1 - (16*phi_Der_f*c_gbs*Q)/r_j + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2))/pow(r_j,2) - (16*phi_Der_f*c_gbs*P*Ze)/r_j + (128*pow(phi_Der_f,2)*pow(c_gbs,2)*P*Q*Ze)/pow(r_j,2) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/pow(r_j,4) + (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + r_Der_Ze*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,3) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,3))/pow(r_j,4) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,4))/pow(r_j,4)) + r_Der_Al*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/(pow(r_j,3)*Al) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,4))/(pow(r_j,4)*Al) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,5))/(pow(r_j,4)*Al))) + r_Der_Al*((-4*phi_Der_f*c_gbs*Q*pow(Ze,2))/r_j + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*pow(Ze,2))/pow(r_j,2) - (4*phi_Der_f*c_gbs*P*pow(Ze,3))/r_j + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*P*Q*pow(Ze,3))/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*SE_LL_TR*pow(Ze,3))/(pow(r_j,2)*Al) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/pow(r_j,4) + (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,4))/pow(r_j,2) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + (256*pow(phi_Der_f,3)*r_Der_P*pow(c_gbs,3)*pow(Ze,5))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*P*Q*pow(Ze,5))/pow(r_j,4) + r_Der_Ze*((64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,3) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,3))/pow(r_j,4) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,4))/pow(r_j,4) - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,5))/pow(r_j,3) + (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,5))/pow(r_j,4) + (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,6))/pow(r_j,4)))
	;
}
/*==========================================================================*/
static double free_evolution_Ze_forward_jac(
	double c_gbs, 
	double dt, double dr, double r_j, 
	double Al, double Ze, double Q,   double P,
	double phi_Der_f, double phiphi_Der_f,
	double SE_LL_TR,        double SE_LL_ThTh,
	double Ze_Der_SE_LL_TR, double Ze_Der_SE_LL_ThTh,
	double r_Der_Al, double r_Der_Ze, double r_Der_Q, double r_Der_P,
	double t_Der_Ze)
{
	return 
		(-2*phi_Der_f*r_Der_P*c_gbs*Al*(r_j - 8*phi_Der_f*c_gbs*Q - 16*phi_Der_f*c_gbs*P*Ze))/pow(r_j,2) - (4*phi_Der_f*r_Der_Q*c_gbs*Al*Ze*(r_j - 8*phi_Der_f*c_gbs*Q - 12*phi_Der_f*c_gbs*P*Ze))/pow(r_j,2) + (128*pow(phi_Der_f,2)*pow(r_Der_Ze,2)*pow(c_gbs,2)*Al*pow(Ze,3)*(-2*r_j + 12*phi_Der_f*c_gbs*Q + 15*phi_Der_f*c_gbs*P*Ze))/pow(r_j,4) + (128*pow(phi_Der_f,3)*pow(r_Der_Al,2)*pow(c_gbs,3)*pow(Ze,4)*(5*P + 6*Q*Ze))/(pow(r_j,4)*Al) + t_Der_Ze*((512*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,3))/pow(r_j,4) + (128*pow(phi_Der_f,2)*r_Der_Al*pow(c_gbs,2)*pow(Ze,3)*(2*r_j - 16*phi_Der_f*c_gbs*Q - 15*phi_Der_f*c_gbs*P*Ze))/(pow(r_j,4)*Al) + (192*pow(phi_Der_f,2)*r_Der_Ze*pow(c_gbs,2)*pow(Ze,2)*(r_j - 8*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))/pow(r_j,4) - (8*phi_Der_f*c_gbs*(pow(r_j,2)*P*(r_j - 8*phi_Der_f*c_gbs*Q) - 8*phi_Der_f*pow(r_j,2)*c_gbs*pow(P,2)*Ze + 8*phi_Der_f*c_gbs*(1 - 8*phiphi_Der_f*c_gbs*pow(Q,2))*pow(Ze,3)))/pow(r_j,4)) + ((256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (128*pow(phi_Der_f,2)*r_Der_Ze*pow(c_gbs,2)*pow(Ze,3)*(r_j - 8*phi_Der_f*c_gbs*Q - 6*phi_Der_f*c_gbs*P*Ze))/pow(r_j,4) + (128*pow(phi_Der_f,2)*r_Der_Al*pow(c_gbs,2)*pow(Ze,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 6*phi_Der_f*c_gbs*P*Ze))/(pow(r_j,4)*Al) + (pow(r_j,4) - 16*phi_Der_f*pow(r_j,3)*c_gbs*P*Ze + 64*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2) - 32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4) - 16*phi_Der_f*pow(r_j,2)*c_gbs*Q*(r_j - 8*phi_Der_f*c_gbs*P*Ze) + 64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*(pow(r_j,2) + 4*phiphi_Der_f*c_gbs*pow(Ze,4)))/pow(r_j,4))/dt + (pow(r_j,5)*(r_j - 8*phi_Der_f*c_gbs*Q)*SE_LL_TR + Ze*(pow(r_j,5)*Ze_Der_SE_LL_TR*(-r_j + 8*phi_Der_f*c_gbs*Q + 8*phi_Der_f*c_gbs*P*Ze) + 2*Al*Ze*(-4*phiphi_Der_f*pow(r_j,3)*c_gbs*P*Q*(r_j - 8*phi_Der_f*c_gbs*Q) - pow(r_j,2)*(pow(r_j,2)*(1 + 8*phiphi_Der_f*c_gbs*pow(P,2)) + 8*phi_Der_f*r_j*c_gbs*(1 - 16*phiphi_Der_f*c_gbs*pow(P,2))*Q - 128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2))*Ze + 12*phi_Der_f*pow(r_j,2)*c_gbs*P*(r_j*(-1 + 8*phiphi_Der_f*c_gbs*pow(P,2)) + 32*phi_Der_f*c_gbs*Q)*pow(Ze,2) + 128*pow(phi_Der_f,2)*pow(c_gbs,2)*(2*pow(r_j,2)*pow(P,2) + SE_LL_ThTh)*pow(Ze,3) + 32*pow(phi_Der_f,2)*Ze_Der_SE_LL_ThTh*pow(c_gbs,2)*pow(Ze,4))))/(4.*pow(r_j,5)*pow(Ze,2)) + ((192*pow(phi_Der_f,3)*r_Der_P*pow(c_gbs,3)*Al*pow(Ze,4))/pow(r_j,4) + (192*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*Al*pow(Ze,5))/pow(r_j,4) - (96*pow(phi_Der_f,2)*t_Der_Ze*pow(c_gbs,2)*pow(Ze,3)*(r_j - 8*phi_Der_f*c_gbs*Q - 6*phi_Der_f*c_gbs*P*Ze))/pow(r_j,4) + (192*pow(phi_Der_f,2)*r_Der_Ze*pow(c_gbs,2)*Al*pow(Ze,4)*(r_j - 6*phi_Der_f*c_gbs*Q - 6*phi_Der_f*c_gbs*P*Ze))/pow(r_j,4) - (48*pow(phi_Der_f,2)*r_Der_Al*pow(c_gbs,2)*pow(Ze,3)*(r_j - 8*phi_Der_f*c_gbs*Q - 4*phi_Der_f*c_gbs*P*Ze - 2*(r_j - 8*phi_Der_f*c_gbs*Q)*pow(Ze,2) + 12*phi_Der_f*c_gbs*P*pow(Ze,3)))/pow(r_j,4) + (3*Ze*(32*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*SE_LL_TR*Ze + Al*(pow(r_j,4) - 12*phi_Der_f*pow(r_j,3)*c_gbs*P*Ze + 32*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2) - 32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4) + 32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*(pow(r_j,2) + 8*phiphi_Der_f*c_gbs*pow(Ze,4)) + 4*phi_Der_f*c_gbs*Q*(-3*pow(r_j,3) + 16*phi_Der_f*c_gbs*P*Ze*(pow(r_j,2) + 4*phiphi_Der_f*c_gbs*pow(Ze,2))))))/(4.*pow(r_j,4)))/dr + r_Der_Ze*((-512*pow(phi_Der_f,3)*r_Der_P*pow(c_gbs,3)*Al*pow(Ze,3))/pow(r_j,4) - (640*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*Al*pow(Ze,4))/pow(r_j,4) - (32*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*Ze*(2*SE_LL_TR + Ze_Der_SE_LL_TR*Ze) + Al*(pow(r_j,4) - 24*phi_Der_f*pow(r_j,3)*c_gbs*P*Ze + 96*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2) - 160*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4) + 32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*(pow(r_j,2) + 40*phiphi_Der_f*c_gbs*pow(Ze,4)) + 4*phi_Der_f*c_gbs*Q*(-3*pow(r_j,3) + 32*phi_Der_f*c_gbs*P*Ze*(pow(r_j,2) + 8*phiphi_Der_f*c_gbs*pow(Ze,2)))))/(2.*pow(r_j,4))) + r_Der_Al*((512*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,3))/pow(r_j,4) + (640*pow(phi_Der_f,3)*r_Der_P*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (32*pow(phi_Der_f,2)*r_Der_Ze*pow(c_gbs,2)*pow(Ze,2)*(3*(r_j - 8*phi_Der_f*c_gbs*Q) - 16*phi_Der_f*c_gbs*P*Ze - 10*(r_j - 8*phi_Der_f*c_gbs*Q)*pow(Ze,2) + 72*phi_Der_f*c_gbs*P*pow(Ze,3)))/pow(r_j,4) + (2*phi_Der_f*c_gbs*Ze*(-8*phi_Der_f*pow(r_j,2)*c_gbs*Ze*(3*SE_LL_TR + Ze_Der_SE_LL_TR*Ze) + Al*(Ze*(-3*pow(r_j,3)*P - 32*phi_Der_f*c_gbs*Ze + 32*phi_Der_f*pow(r_j,2)*c_gbs*pow(P,2)*Ze) + 16*phi_Der_f*c_gbs*pow(Q,2)*(pow(r_j,2) + 16*phiphi_Der_f*c_gbs*pow(Ze,2)) + Q*(-2*pow(r_j,3) + 16*phi_Der_f*c_gbs*P*Ze*(3*pow(r_j,2) + 20*phiphi_Der_f*c_gbs*pow(Ze,2))))))/(pow(r_j,4)*Al))
	;
}
/*==========================================================================*/
static double compute_iteration_excision_boundary_condition_Ze(
	double s_L,	double c_gbs,
	int Nx,
	double dt, 	double dx,
	int exc_jC,
	double bbox[2],
	double *Al_n,  double *Al_nm1, double *Ze_n, double *Ze_nm1,
	double *P_n,   double *P_nm1,  double *Q_n,  double *Q_nm1)
{
        double x_j = dx * exc_jC ;
        double r_j = stereographic_r( s_L, x_j) ;
        double dr  = stereographic_dr(s_L, x_j, dx) ;

        double Al     = (Al_n[exc_jC  ] + Al_nm1[exc_jC  ]) / 2. ;
        double Al_jp1 = (Al_n[exc_jC+1] + Al_nm1[exc_jC+1]) / 2. ;
        double Al_jp2 = (Al_n[exc_jC+2] + Al_nm1[exc_jC+2]) / 2. ;

        double Ze     = (Ze_n[exc_jC  ] + Ze_nm1[exc_jC  ]) / 2. ;
        double Ze_jp1 = (Ze_n[exc_jC+1] + Ze_nm1[exc_jC+1]) / 2. ;
        double Ze_jp2 = (Ze_n[exc_jC+2] + Ze_nm1[exc_jC+2]) / 2. ;

        double P     = (P_n[exc_jC  ] + P_nm1[exc_jC  ]) / 2. ;
        double P_jp1 = (P_n[exc_jC+1] + P_nm1[exc_jC+1]) / 2. ;
        double P_jp2 = (P_n[exc_jC+2] + P_nm1[exc_jC+2]) / 2. ;

        double Q     = (Q_n[exc_jC  ] + Q_nm1[exc_jC  ]) / 2. ;
        double Q_jp1 = (Q_n[exc_jC+1] + Q_nm1[exc_jC+1]) / 2. ;
        double Q_jp2 = (Q_n[exc_jC+2] + Q_nm1[exc_jC+2]) / 2. ;

        double t_Der_Ze = D1_CrankNicolson_2ndOrder(Ze_n[exc_jC], Ze_nm1[exc_jC], dt) ;

        double	r_Der_Al = D1_forward_2ndOrder(Al_jp2, Al_jp1, Al, dr) ;
        double	r_Der_Ze = D1_forward_2ndOrder(Ze_jp2, Ze_jp1, Ze, dr) ;
        double	r_Der_P  = D1_forward_2ndOrder( P_jp2,  P_jp1, P,  dr) ;
        double	r_Der_Q  = D1_forward_2ndOrder( Q_jp2,  Q_jp1, Q,  dr) ;

	double SE_LL_ThTh = (1./2) * pow(r_j,2) * (pow(P,2) - pow(Q,2)) ;
	double Ze_Der_SE_LL_ThTh = 0. ;

	double SE_LL_TR        = (Al*(2*P*Q + (pow(P,2) + pow(Q,2))*Ze  ))/2. ;
	double Ze_Der_SE_LL_TR = (Al*(0     + (pow(P,2) + pow(Q,2))*0.5 ))/2. ;

	double phi_Der_f    = 1. ;
	double phiphi_Der_f = 0. ;
	
	double res_Ze = free_evolution_Ze_res(
		c_gbs, 
		r_j, 
		Al, Ze, Q, P,
		phi_Der_f, phiphi_Der_f,
		SE_LL_TR, SE_LL_ThTh,
		r_Der_Al, r_Der_Ze, r_Der_Q, r_Der_P,
		t_Der_Ze)
	;
	double jac_Ze = free_evolution_Ze_forward_jac(
		c_gbs, 
		dt, dr, r_j, 
		Al, Ze, Q,   P,
		phi_Der_f, phiphi_Der_f,
		SE_LL_TR,        SE_LL_ThTh,
		Ze_Der_SE_LL_TR, Ze_Der_SE_LL_ThTh,
		r_Der_Al, r_Der_Ze, r_Der_Q, r_Der_P,
		t_Der_Ze)
	;
	Ze_n[exc_jC] -= res_Ze/jac_Ze ;

	return fabs(res_Ze) ; 
}
/*==========================================================================*/
void solve_Al_Ze_EdGB(
	double s_L,	double c_gbs,
	int Nx,
	double dt, 	double dx,
	double err_tolerance,
	bool excision_on,
	int start_jC,
	double bbox[2],
	double *Al_n, 	double *Al_nm1, double *Ze_n, double *Ze_nm1,
	double * P_n, 	double * P_nm1, double * Q_n, double * Q_nm1)
{
	if (start_jC==Nx-1) { /* if interior to the excision point */
		return ;
	}
	double res= 0 ;
	do {
		res = 0 ;
		if ((excision_on==true)
		&&  (start_jC>0)
		) {	
			res += compute_iteration_excision_boundary_condition_Ze(
				s_L,	c_gbs,
				Nx,
				dt, 	dx,
				start_jC,
				bbox,
				Al_n,  Al_nm1, Ze_n, Ze_nm1,
				 P_n,   P_nm1,  Q_n,  Q_nm1)
			;
		}
		res += compute_iteration_Ze(
			s_L,	c_gbs,
			Nx,
			dt, 	dx,
			start_jC,
			bbox,
			Al_n, 	Ze_n, 
			P_n, 	Q_n)
		;
		res += compute_iteration_Al(
			s_L,	c_gbs,
			Nx,
			dt, 	dx,
			start_jC,
			bbox,
			Al_n, 	Ze_n, 
			P_n, 	Q_n)
		;
	} while (res>err_tolerance) ;
	return ;
}
/*###########################################################################*/
static double compute_res_P(
	double c_gbs, double r_j,
	double phi_Der_f, double phiphi_Der_f, 
	double SE_LL_TR, double SE_LL_ThTh, 
	double Al,       double Ze,       double Q,       double P,
	double r_Der_Al, double r_Der_Ze, double r_Der_Q, double r_Der_P,
                                                      double t_Der_P)
{
	return 
		(-2*Al*Q)/r_j + (32*phi_Der_f*c_gbs*Al*pow(Q,2))/pow(r_j,2) - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,3))/pow(r_j,3) - (2*Al*P*Ze)/r_j + (64*phi_Der_f*c_gbs*Al*P*Q*Ze)/pow(r_j,2) - (384*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Q,2)*Ze)/pow(r_j,3) + (32*phi_Der_f*c_gbs*Al*pow(P,2)*pow(Ze,2))/pow(r_j,2) - (384*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*Q*pow(Ze,2))/pow(r_j,3) - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,3)*pow(Ze,3))/pow(r_j,3) - (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*P*Q*pow(Ze,3))/pow(r_j,4) + (256*pow(phiphi_Der_f,2)*phi_Der_f*pow(c_gbs,3)*Al*P*pow(Q,3)*pow(Ze,3))/pow(r_j,4) - (4*phi_Der_f*c_gbs*Al*pow(Ze,4))/pow(r_j,4) - (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,4))/pow(r_j,4) + (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + (256*pow(phiphi_Der_f,2)*phi_Der_f*pow(c_gbs,3)*Al*pow(P,2)*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + SE_LL_TR*((-4*phi_Der_f*c_gbs*Ze)/pow(r_j,2) + (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*pow(Q,2)*Ze)/pow(r_j,2)) + SE_LL_ThTh*((8*phi_Der_f*c_gbs*Al*pow(Ze,2))/pow(r_j,4) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,5) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,5)) + pow(r_Der_Ze,2)*((64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,3) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,2))/pow(r_j,4) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*P*Q*pow(Ze,3))/pow(r_j,4)) + r_Der_Q*(-Al + (16*phi_Der_f*c_gbs*Al*Q)/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,2))/pow(r_j,2) + (16*phi_Der_f*c_gbs*Al*P*Ze)/r_j - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*Q*Ze)/pow(r_j,2) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*SE_LL_TR*Ze)/pow(r_j,2) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,2))/pow(r_j,2) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*P*Q*pow(Ze,3))/pow(r_j,4) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,4))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(P,2)*pow(Ze,4))/pow(r_j,4)) + pow(r_Der_Al,2)*((64*pow(phi_Der_f,2)*pow(c_gbs,2)*P*pow(Ze,3))/(pow(r_j,3)*Al) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*P*Q*pow(Ze,3))/(pow(r_j,4)*Al) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(P,2)*pow(Ze,4))/(pow(r_j,4)*Al) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Q*pow(Ze,4))/(pow(r_j,3)*Al) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,4))/(pow(r_j,4)*Al) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*P*Q*pow(Ze,5))/(pow(r_j,4)*Al)) + r_Der_P*(-(Al*Ze) + (16*phi_Der_f*c_gbs*Al*Q*Ze)/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,2)*Ze)/pow(r_j,2) + (16*phi_Der_f*c_gbs*Al*P*pow(Ze,2))/r_j - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*Q*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,3))/pow(r_j,4) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,3))/pow(r_j,2) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,3))/pow(r_j,4) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,5))/pow(r_j,4) - (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,5))/pow(r_j,4) + r_Der_Q*((256*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*pow(Ze,3))/pow(r_j,4) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*pow(Ze,5))/pow(r_j,4))) + t_Der_P*(1 - (16*phi_Der_f*c_gbs*Q)/r_j + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2))/pow(r_j,2) - (16*phi_Der_f*c_gbs*P*Ze)/r_j + (128*pow(phi_Der_f,2)*pow(c_gbs,2)*P*Q*Ze)/pow(r_j,2) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/pow(r_j,4) + (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + r_Der_Ze*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,3) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,3))/pow(r_j,4) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,4))/pow(r_j,4)) + r_Der_Al*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/(pow(r_j,3)*Al) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,4))/(pow(r_j,4)*Al) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,5))/(pow(r_j,4)*Al))) + r_Der_Ze*(-(Al*P) + (16*phi_Der_f*c_gbs*Al*P*Q)/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Q,2))/pow(r_j,2) + (16*phi_Der_f*c_gbs*Al*pow(P,2)*Ze)/r_j - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*Q*Ze)/pow(r_j,2) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,3)*pow(Ze,2))/pow(r_j,2) + (64*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*P*Q*pow(Ze,2))/pow(r_j,3) - (512*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*P*pow(Q,2)*pow(Ze,2))/pow(r_j,4) + (16*phi_Der_f*c_gbs*Al*pow(Ze,3))/pow(r_j,3) + (128*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,3))/pow(r_j,3) - (160*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*pow(Ze,3))/pow(r_j,4) + (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*Al*Q*pow(Ze,3))/pow(r_j,4) - (1280*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(P,2)*Q*pow(Ze,3))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(Q,3)*pow(Ze,3))/pow(r_j,4) - (96*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Ze,4))/pow(r_j,4) - (768*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(P,3)*pow(Ze,4))/pow(r_j,4) + SE_LL_TR*((8*phi_Der_f*c_gbs)/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Q)/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*P*Ze)/pow(r_j,2)) + r_Der_P*((64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,2))/pow(r_j,3) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*Q*pow(Ze,2))/pow(r_j,4) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*P*pow(Ze,3))/pow(r_j,4) - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,4))/pow(r_j,3) + (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*Q*pow(Ze,4))/pow(r_j,4) + (768*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*P*pow(Ze,5))/pow(r_j,4))) + r_Der_Al*(-Q + (16*phi_Der_f*c_gbs*pow(Q,2))/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,3))/pow(r_j,2) - P*Ze + (32*phi_Der_f*c_gbs*P*Q*Ze)/r_j - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*P*pow(Q,2)*Ze)/pow(r_j,2) - (8*phi_Der_f*c_gbs*pow(Ze,2))/pow(r_j,3) + (16*phi_Der_f*c_gbs*pow(P,2)*pow(Ze,2))/r_j + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Q*pow(Ze,2))/pow(r_j,4) - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*Q*pow(Ze,2))/pow(r_j,2) + (64*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*pow(Q,2)*pow(Ze,2))/pow(r_j,3) - (512*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,3)*pow(Ze,2))/pow(r_j,4) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*P*pow(Ze,3))/pow(r_j,4) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,3)*pow(Ze,3))/pow(r_j,2) + (192*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*P*Q*pow(Ze,3))/pow(r_j,3) - (2048*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*P*pow(Q,2)*pow(Ze,3))/pow(r_j,4) + (16*phi_Der_f*c_gbs*pow(Ze,4))/pow(r_j,3) + (128*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*pow(P,2)*pow(Ze,4))/pow(r_j,3) - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Q*pow(Ze,4))/pow(r_j,4) - (2304*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(P,2)*Q*pow(Ze,4))/pow(r_j,4) - (96*pow(phi_Der_f,2)*pow(c_gbs,2)*P*pow(Ze,5))/pow(r_j,4) - (768*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(P,3)*pow(Ze,5))/pow(r_j,4) + SE_LL_TR*((8*phi_Der_f*c_gbs*Ze)/(r_j*Al) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Q*Ze)/(pow(r_j,2)*Al) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*P*pow(Ze,2))/(pow(r_j,2)*Al)) + r_Der_Q*((64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,2))/pow(r_j,3) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,2))/pow(r_j,4) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,3))/pow(r_j,4)) + r_Der_Ze*((16*phi_Der_f*c_gbs*Ze)/pow(r_j,2) - (256*pow(phi_Der_f,2)*pow(c_gbs,2)*Q*Ze)/pow(r_j,3) + (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Q,2)*Ze)/pow(r_j,4) - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*P*pow(Ze,2))/pow(r_j,3) + (1536*pow(phi_Der_f,3)*pow(c_gbs,3)*P*Q*pow(Ze,2))/pow(r_j,4) + (512*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(P,2)*pow(Ze,3))/pow(r_j,4) + (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Q*pow(Ze,3))/pow(r_j,3) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,3))/pow(r_j,4) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*Q*pow(Ze,4))/pow(r_j,4)) + r_Der_P*((192*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,3) - (1536*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,3))/pow(r_j,4) - (1280*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,4))/pow(r_j,4) - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,5))/pow(r_j,3) + (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,5))/pow(r_j,4) + (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,6))/pow(r_j,4)))
	;
}
/*###########################################################################*/
static double compute_jac_center_P(
	double c_gbs, double r_j, double dr, double dt,
	double phi_Der_f, double phiphi_Der_f, 
	double SE_LL_TR,       double SE_LL_ThTh, 
	double P_Der_SE_LL_TR, double P_Der_SE_LL_ThTh, 
	double Al,       double Ze,       double Q,       double P,
	double r_Der_Al, double r_Der_Ze, double r_Der_Q, double r_Der_P,
                                                      double t_Der_P)
{
	return 
		(-2*phi_Der_f*P_Der_SE_LL_TR*c_gbs*Ze)/pow(r_j,2) - (Al*Ze)/r_j + (32*phi_Der_f*c_gbs*Al*Q*Ze)/pow(r_j,2) + (16*phiphi_Der_f*phi_Der_f*P_Der_SE_LL_TR*pow(c_gbs,2)*pow(Q,2)*Ze)/pow(r_j,2) - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,2)*Ze)/pow(r_j,3) + (4*phi_Der_f*P_Der_SE_LL_ThTh*c_gbs*Al*pow(Ze,2))/pow(r_j,4) + (32*phi_Der_f*c_gbs*Al*P*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*P_Der_SE_LL_ThTh*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,5) - (384*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*Q*pow(Ze,2))/pow(r_j,3) - (32*pow(phi_Der_f,2)*P_Der_SE_LL_ThTh*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,5) - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,3))/pow(r_j,3) - (16*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*Q*pow(Ze,3))/pow(r_j,4) - (128*pow(phi_Der_f,3)*pow(r_Der_Ze,2)*pow(c_gbs,3)*Al*Q*pow(Ze,3))/pow(r_j,4) + (128*pow(phiphi_Der_f,2)*phi_Der_f*pow(c_gbs,3)*Al*pow(Q,3)*pow(Ze,3))/pow(r_j,4) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*SE_LL_ThTh*pow(Ze,3))/pow(r_j,5) - (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*P*pow(Ze,4))/pow(r_j,4) + (256*pow(phiphi_Der_f,2)*phi_Der_f*pow(c_gbs,3)*Al*P*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + r_Der_P*((8*phi_Der_f*c_gbs*Al*pow(Ze,2))/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,2) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,2)) + r_Der_Q*((16*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*Ze)/pow(r_j,2) + (8*phi_Der_f*c_gbs*Al*Ze)/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*Ze)/pow(r_j,2) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Ze,2))/pow(r_j,2) + (128*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*Q*pow(Ze,3))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*P*pow(Ze,4))/pow(r_j,4)) + t_Der_P*((-8*phi_Der_f*c_gbs*Ze)/r_j + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Q*Ze)/pow(r_j,2) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*P*pow(Ze,2))/pow(r_j,2) - (384*pow(phi_Der_f,3)*r_Der_Ze*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) - (384*pow(phi_Der_f,3)*r_Der_Al*pow(c_gbs,3)*pow(Ze,5))/(pow(r_j,4)*Al)) + pow(r_Der_Al,2)*((32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/(pow(r_j,3)*Al) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,3))/(pow(r_j,4)*Al) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,4))/(pow(r_j,4)*Al) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,5))/(pow(r_j,4)*Al)) + r_Der_Ze*((4*phi_Der_f*P_Der_SE_LL_TR*c_gbs)/r_j - Al/2. - (32*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*Q)/pow(r_j,2) + (8*phi_Der_f*c_gbs*Al*Q)/r_j - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,2))/pow(r_j,2) - (16*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*P*Ze)/pow(r_j,2) + (16*phi_Der_f*c_gbs*Al*P*Ze)/r_j - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*Q*Ze)/pow(r_j,2) - (16*pow(phi_Der_f,2)*pow(c_gbs,2)*SE_LL_TR*Ze)/pow(r_j,2) - (96*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,2))/pow(r_j,2) + (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,3) - (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,2))/pow(r_j,4) + (128*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,3) - (1280*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*P*Q*pow(Ze,3))/pow(r_j,4) - (48*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,4))/pow(r_j,4) - (1152*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(P,2)*pow(Ze,4))/pow(r_j,4) + r_Der_P*((-128*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*pow(Ze,3))/pow(r_j,4) + (384*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*pow(Ze,5))/pow(r_j,4))) + (1 - (16*phi_Der_f*c_gbs*Q)/r_j + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2))/pow(r_j,2) - (16*phi_Der_f*c_gbs*P*Ze)/r_j + (128*pow(phi_Der_f,2)*pow(c_gbs,2)*P*Q*Ze)/pow(r_j,2) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/pow(r_j,4) + (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + r_Der_Ze*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,3) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,3))/pow(r_j,4) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,4))/pow(r_j,4)) + r_Der_Al*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/(pow(r_j,3)*Al) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,4))/(pow(r_j,4)*Al) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,5))/(pow(r_j,4)*Al)))/dt + r_Der_Al*(-Ze/2. + (4*phi_Der_f*P_Der_SE_LL_TR*c_gbs*Ze)/(r_j*Al) + (16*phi_Der_f*c_gbs*Q*Ze)/r_j - (32*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*Q*Ze)/(pow(r_j,2)*Al) - (96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*Ze)/pow(r_j,2) + (16*phi_Der_f*c_gbs*P*pow(Ze,2))/r_j - (16*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*P*pow(Ze,2))/(pow(r_j,2)*Al) - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*P*Q*pow(Ze,2))/pow(r_j,2) - (16*pow(phi_Der_f,2)*pow(c_gbs,2)*SE_LL_TR*pow(Ze,2))/(pow(r_j,2)*Al) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,4) - (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,3))/pow(r_j,4) - (96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,3))/pow(r_j,2) + (96*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Q*pow(Ze,3))/pow(r_j,3) - (1024*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,3))/pow(r_j,4) + (128*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*P*pow(Ze,4))/pow(r_j,3) - (2304*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*P*Q*pow(Ze,4))/pow(r_j,4) - (48*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,5))/pow(r_j,4) - (1152*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(P,2)*pow(Ze,5))/pow(r_j,4) + r_Der_Ze*((-96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,2))/pow(r_j,3) + (768*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,2))/pow(r_j,4) + (512*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,3))/pow(r_j,4) - (384*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,4))/pow(r_j,4)) + r_Der_P*((-640*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (384*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Ze,6))/pow(r_j,4)))
	;
}
/*==========================================================================*/
static double compute_jac_forward_P(
	double c_gbs, double r_j, double dr, double dt,
	double phi_Der_f, double phiphi_Der_f, 
	double SE_LL_TR,       double SE_LL_ThTh, 
	double P_Der_SE_LL_TR, double P_Der_SE_LL_ThTh, 
	double Al,       double Ze,       double Q,       double P,
	double r_Der_Al, double r_Der_Ze, double r_Der_Q, double r_Der_P,
                                                      double t_Der_P)
{
	return 
		(-2*phi_Der_f*P_Der_SE_LL_TR*c_gbs*Ze)/pow(r_j,2) - (Al*Ze)/r_j + (32*phi_Der_f*c_gbs*Al*Q*Ze)/pow(r_j,2) + (16*phiphi_Der_f*phi_Der_f*P_Der_SE_LL_TR*pow(c_gbs,2)*pow(Q,2)*Ze)/pow(r_j,2) - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,2)*Ze)/pow(r_j,3) + (4*phi_Der_f*P_Der_SE_LL_ThTh*c_gbs*Al*pow(Ze,2))/pow(r_j,4) + (32*phi_Der_f*c_gbs*Al*P*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*P_Der_SE_LL_ThTh*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,5) - (384*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*Q*pow(Ze,2))/pow(r_j,3) - (32*pow(phi_Der_f,2)*P_Der_SE_LL_ThTh*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,5) - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,3))/pow(r_j,3) - (16*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*Q*pow(Ze,3))/pow(r_j,4) - (128*pow(phi_Der_f,3)*pow(r_Der_Ze,2)*pow(c_gbs,3)*Al*Q*pow(Ze,3))/pow(r_j,4) + (128*pow(phiphi_Der_f,2)*phi_Der_f*pow(c_gbs,3)*Al*pow(Q,3)*pow(Ze,3))/pow(r_j,4) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*SE_LL_ThTh*pow(Ze,3))/pow(r_j,5) - (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*P*pow(Ze,4))/pow(r_j,4) + (256*pow(phiphi_Der_f,2)*phi_Der_f*pow(c_gbs,3)*Al*P*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + r_Der_P*((8*phi_Der_f*c_gbs*Al*pow(Ze,2))/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,2) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,2)) + r_Der_Q*((16*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*Ze)/pow(r_j,2) + (8*phi_Der_f*c_gbs*Al*Ze)/r_j - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*Ze)/pow(r_j,2) - (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*pow(Ze,2))/pow(r_j,2) + (128*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*Q*pow(Ze,3))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*P*pow(Ze,4))/pow(r_j,4)) + t_Der_P*((-8*phi_Der_f*c_gbs*Ze)/r_j + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*Q*Ze)/pow(r_j,2) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*P*pow(Ze,2))/pow(r_j,2) - (384*pow(phi_Der_f,3)*r_Der_Ze*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) - (384*pow(phi_Der_f,3)*r_Der_Al*pow(c_gbs,3)*pow(Ze,5))/(pow(r_j,4)*Al)) + pow(r_Der_Al,2)*((32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/(pow(r_j,3)*Al) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,3))/(pow(r_j,4)*Al) - (512*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,4))/(pow(r_j,4)*Al) - (256*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,5))/(pow(r_j,4)*Al)) + r_Der_Ze*((4*phi_Der_f*P_Der_SE_LL_TR*c_gbs)/r_j - Al/2. - (32*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*Q)/pow(r_j,2) + (8*phi_Der_f*c_gbs*Al*Q)/r_j - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Q,2))/pow(r_j,2) - (16*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*P*Ze)/pow(r_j,2) + (16*phi_Der_f*c_gbs*Al*P*Ze)/r_j - (128*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*P*Q*Ze)/pow(r_j,2) - (16*pow(phi_Der_f,2)*pow(c_gbs,2)*SE_LL_TR*Ze)/pow(r_j,2) - (96*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,2))/pow(r_j,2) + (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,3) - (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,2))/pow(r_j,4) + (128*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,3) - (1280*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*P*Q*pow(Ze,3))/pow(r_j,4) - (48*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*pow(Ze,4))/pow(r_j,4) - (1152*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*Al*pow(P,2)*pow(Ze,4))/pow(r_j,4) + r_Der_P*((-128*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*pow(Ze,3))/pow(r_j,4) + (384*pow(phi_Der_f,3)*pow(c_gbs,3)*Al*pow(Ze,5))/pow(r_j,4))) + (1 - (16*phi_Der_f*c_gbs*Q)/r_j + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2))/pow(r_j,2) - (16*phi_Der_f*c_gbs*P*Ze)/r_j + (128*pow(phi_Der_f,2)*pow(c_gbs,2)*P*Q*Ze)/pow(r_j,2) + (64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2))/pow(r_j,2) - (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/pow(r_j,4) + (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + r_Der_Ze*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,3) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,3))/pow(r_j,4) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,4))/pow(r_j,4)) + r_Der_Al*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/(pow(r_j,3)*Al) - (1024*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,4))/(pow(r_j,4)*Al) - (768*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,5))/(pow(r_j,4)*Al)))/dt + r_Der_Al*(-Ze/2. + (4*phi_Der_f*P_Der_SE_LL_TR*c_gbs*Ze)/(r_j*Al) + (16*phi_Der_f*c_gbs*Q*Ze)/r_j - (32*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*Q*Ze)/(pow(r_j,2)*Al) - (96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*Ze)/pow(r_j,2) + (16*phi_Der_f*c_gbs*P*pow(Ze,2))/r_j - (16*pow(phi_Der_f,2)*P_Der_SE_LL_TR*pow(c_gbs,2)*P*pow(Ze,2))/(pow(r_j,2)*Al) - (192*pow(phi_Der_f,2)*pow(c_gbs,2)*P*Q*pow(Ze,2))/pow(r_j,2) - (16*pow(phi_Der_f,2)*pow(c_gbs,2)*SE_LL_TR*pow(Ze,2))/(pow(r_j,2)*Al) + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,4) - (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,3))/pow(r_j,4) - (96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,3))/pow(r_j,2) + (96*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*Q*pow(Ze,3))/pow(r_j,3) - (1024*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,3))/pow(r_j,4) + (128*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*P*pow(Ze,4))/pow(r_j,3) - (2304*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*P*Q*pow(Ze,4))/pow(r_j,4) - (48*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,5))/pow(r_j,4) - (1152*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(P,2)*pow(Ze,5))/pow(r_j,4) + r_Der_Ze*((-96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,2))/pow(r_j,3) + (768*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,2))/pow(r_j,4) + (512*pow(phi_Der_f,3)*pow(c_gbs,3)*P*pow(Ze,3))/pow(r_j,4) - (384*pow(phi_Der_f,3)*pow(c_gbs,3)*Q*pow(Ze,4))/pow(r_j,4)) + r_Der_P*((-640*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (384*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Ze,6))/pow(r_j,4)))
	;
}
/*==========================================================================*/
static double compute_iteration_Crank_Nicolson_PQ(
	double s_L,	double c_gbs,
	int Nx,
	double dt, 	double dx,
	bool excision_on,
	int exc_jC,
	double bbox[2],
	bool perim_interior[2],
	double *Al_n, 	double *Al_nm1, double *Ze_n, double *Ze_nm1,
	double *P_n, 	double  *P_nm1, double  *Q_n, double  *Q_nm1)
{
	double lower_x = bbox[0] ;

	int size = 0 ;
	if (fabs(bbox[1]-s_L)<machine_epsilon) size = Nx-1 ; /* to avoid problems with r=infty when x=s_L */
	else size = Nx ;
	double res_infty_norm = 0 ; /* returning this */

	for (int jC=exc_jC+1;jC<size-1;jC++) {
		double x_j = lower_x + (dx * (jC)  ) ;

		double r_j = stereographic_r(s_L, x_j  ) ;
		double dr = stereographic_dr(s_L, x_j, dx) ;

		double Al     = (Al_n[jC  ]+Al_nm1[jC  ])/2 ;
		double Al_jp1 = (Al_n[jC+1]+Al_nm1[jC+1])/2 ;
		double Al_jm1 = (Al_n[jC-1]+Al_nm1[jC-1])/2 ;

		double Ze     = (Ze_n[jC  ]+Ze_nm1[jC  ])/2 ;
		double Ze_jp1 = (Ze_n[jC+1]+Ze_nm1[jC+1])/2 ;
		double Ze_jm1 = (Ze_n[jC-1]+Ze_nm1[jC-1])/2 ;

		double P     = (P_n[jC  ]+P_nm1[jC  ])/2 ;
		double P_jp1 = (P_n[jC+1]+P_nm1[jC+1])/2 ;
		double P_jm1 = (P_n[jC-1]+P_nm1[jC-1])/2 ;

		double Q     = (Q_n[jC  ]+Q_nm1[jC  ])/2 ;
		double Q_jp1 = (Q_n[jC+1]+Q_nm1[jC+1])/2 ;
		double Q_jm1 = (Q_n[jC-1]+Q_nm1[jC-1])/2 ;

		double t_Der_Q = D1_CrankNicolson_2ndOrder(Q_n[jC], Q_nm1[jC], dt) ;
		double t_Der_P = D1_CrankNicolson_2ndOrder(P_n[jC], P_nm1[jC], dt) ;

		double r_Der_Al = D1_center_2ndOrder(Al_jp1, Al_jm1, dr) ;
		double r_Der_Ze = D1_center_2ndOrder(Ze_jp1, Ze_jm1, dr) ;
		double r_Der_P  = D1_center_2ndOrder( P_jp1,  P_jm1, dr) ;
		double r_Der_Q  = D1_center_2ndOrder( Q_jp1,  Q_jm1, dr) ;

		double phi_Der_f = 1 ;
		double phiphi_Der_f = 0 ;

		double SE_LL_TR = (Al*(2*P*Q + (pow(P,2) + pow(Q,2))*Ze))/2. ;
		double SE_LL_ThTh = (pow(r_j,2)*(pow(P,2) - pow(Q,2)))/2. ;

		double P_Der_SE_LL_TR = (Al*(Q + P*Ze))/2. ;
		double P_Der_SE_LL_ThTh = (pow(r_j,2)*P)/2. ;
	/* Q field 
	*/
		double res_Q = t_Der_Q 
		-	D1_center_2ndOrder(
				Al_jp1*(P_jp1 + Ze_jp1*Q_jp1),
				Al_jm1*(P_jm1 + Ze_jm1*Q_jm1),
			dr)
		;
		double jac_Q = (1./dt) ;
	/* P field 
	*/
		double res_P = compute_res_P(
			c_gbs, r_j,
			phi_Der_f, phiphi_Der_f, 
			SE_LL_TR, SE_LL_ThTh, 
			Al,       Ze,       Q,       P,
			r_Der_Al, r_Der_Ze, r_Der_Q, r_Der_P,
			t_Der_P)
		;
		double jac_P = compute_jac_center_P(
			c_gbs, r_j, dr, dt,
			phi_Der_f, phiphi_Der_f, 
			SE_LL_TR,       SE_LL_ThTh, 
			P_Der_SE_LL_TR, P_Der_SE_LL_ThTh, 
			Al,       Ze,       Q,       P,
			r_Der_Al, r_Der_Ze, r_Der_Q, r_Der_P,
			t_Der_P)	
		;
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

		double r_j = stereographic_r(s_L, x_j) ;
		double dr = stereographic_dr(s_L, x_j, dx) ;

		double Al     = (Al_n[exc_jC+0]+Al_nm1[exc_jC+0])/2 ;
		double Al_jp1 = (Al_n[exc_jC+1]+Al_nm1[exc_jC+1])/2 ;
		double Al_jp2 = (Al_n[exc_jC+2]+Al_nm1[exc_jC+2])/2 ;

		double Ze     = (Ze_n[exc_jC+0]+Ze_nm1[exc_jC+0])/2 ;
		double Ze_jp1 = (Ze_n[exc_jC+1]+Ze_nm1[exc_jC+1])/2 ;
		double Ze_jp2 = (Ze_n[exc_jC+2]+Ze_nm1[exc_jC+2])/2 ;

		double P     = (P_n[exc_jC+0]+P_nm1[exc_jC+0])/2 ;
		double P_jp1 = (P_n[exc_jC+1]+P_nm1[exc_jC+1])/2 ;
		double P_jp2 = (P_n[exc_jC+2]+P_nm1[exc_jC+2])/2 ;

		double Q     = (Q_n[exc_jC+0]+Q_nm1[exc_jC+0])/2 ;
		double Q_jp1 = (Q_n[exc_jC+1]+Q_nm1[exc_jC+1])/2 ;
		double Q_jp2 = (Q_n[exc_jC+2]+Q_nm1[exc_jC+2])/2 ;

		double t_Der_P = D1_CrankNicolson_2ndOrder(P_n[exc_jC], P_nm1[exc_jC], dt) ;

		double r_Der_Al = D1_forward_2ndOrder(Al_jp2, Al_jp1, Al, dr) ;
		double r_Der_Ze = D1_forward_2ndOrder(Ze_jp2, Ze_jp1, Ze, dr) ;
		double r_Der_P  = D1_forward_2ndOrder( P_jp2,  P_jp1,  P, dr) ;
		double r_Der_Q  = D1_forward_2ndOrder( Q_jp2,  Q_jp1,  Q, dr) ;

		double phi_Der_f = 1 ;
		double phiphi_Der_f = 0 ;

		double SE_LL_TR = (Al*(2*P*Q + (pow(P,2) + pow(Q,2))*Ze))/2. ;
		double SE_LL_ThTh = (pow(r_j,2)*(pow(P,2) - pow(Q,2)))/2. ;

		double P_Der_SE_LL_TR = (Al*(Q + P*Ze))/2. ;
		double P_Der_SE_LL_ThTh = (pow(r_j,2)*P)/2. ;
	/* Q field 
	*/
		double res_Q = D1_CrankNicolson_2ndOrder(
				Q_n[exc_jC], 
				Q_nm1[exc_jC], 
				dt)
		-	D1_forward_2ndOrder(
				Al_jp2*(P_jp2 + Ze_jp2*Q_jp2),
				Al_jp1*(P_jp1 + Ze_jp1*Q_jp1),
				Al    *(P     + Ze*    Q    ),
				dr)
		;
		double jac_Q =	1/dt 
		- 	(1./2.)*(-3./(2*dr))*Al*Ze 
		;
	/* P field 
	*/
		double res_P = compute_res_P(
			c_gbs, r_j,
			phi_Der_f, phiphi_Der_f, 
			SE_LL_TR, SE_LL_ThTh, 
			Al,       Ze,       Q,       P,
			r_Der_Al, r_Der_Ze, r_Der_Q, r_Der_P,
			t_Der_P)
		;
		double jac_P = compute_jac_forward_P(
			c_gbs, r_j, dr, dt,
			phi_Der_f, phiphi_Der_f, 
			SE_LL_TR,       SE_LL_ThTh, 
			P_Der_SE_LL_TR, P_Der_SE_LL_ThTh, 
			Al,       Ze,       Q,       P,
			r_Der_Al, r_Der_Ze, r_Der_Q, r_Der_P,
			t_Der_P)
		;
	/****/
		Q_n[exc_jC] -= res_Q / jac_Q ;
		P_n[exc_jC] -= res_P / jac_P ;

		res_infty_norm = weighted_infty_norm(1-x_j/s_L, res_Q, res_infty_norm) ;
		res_infty_norm = weighted_infty_norm(1-x_j/s_L, res_P, res_infty_norm) ;
	}
/*--------------------------------------------------------------------------*/
/* Dirichlet outer boundary conditions if outermost level;
 * otherwise do not evolve outer boundary anyways (it is interpolated) */
/*--------------------------------------------------------------------------*/
	return res_infty_norm ;
} 
/*===========================================================================*/
void advance_tStep_massless_scalar_EdGB(
	double s_L,double c_gbs,
	int Nx, 
	double dt, double dx, 
	double err_tolerance, 
	bool excision_on,
	int exc_jC,
	double bbox[2], 
	bool perim_interior[2],
	double *Al_n, double *Al_nm1, double *Ze_n, double *Ze_nm1,
	double *P_n,  double  *P_nm1, double  *Q_n, double  *Q_nm1)
{ 
	if (exc_jC==Nx-1) { /* then no domain */
		return ;
	}
	double res = 0 ;

	do {
		res = compute_iteration_Crank_Nicolson_PQ(
			s_L,	c_gbs,
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

	Kreiss_Oliger_filter(Nx, P_n) ;
	Kreiss_Oliger_filter(Nx, Q_n) ;

	if (fabs(bbox[0])<machine_epsilon) {
		Kreiss_Oliger_filter_origin(P_n, "even") ;
		Kreiss_Oliger_filter_origin(Q_n, "odd") ;
	}
	return ;
}	
