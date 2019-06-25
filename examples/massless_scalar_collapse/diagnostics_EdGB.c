#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "stencils.h"
#include "diagnostics_EdGB.h"
#include "basic_matrix_computations.h"


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
	int Nx,       double dx,
	double s_L,   int exc_jC,
	double x_lower,
	double* Al_n, double* Ze_n,
	double* P_n,  double* Q_n,
	double* SE_LL_TR,
	double* SE_LL_ThTh)
{
	for (int jC=exc_jC; jC<Nx; jC++) {
		double x_j = x_lower + (dx*jC) ;
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
	double s_L, double c_gbs,
	int exc_jC, 
	double x_lower,
	double* Al_n, 
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double*  P_n, double*  P_nm1, double*  P_nm2,
	double*  Q_n,
	double* SE_LL_TR,
	double* eom_TR)
{
	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		double x_j = x_lower + (jC * dx) ;
                double r_j = stereographic_r(s_L, x_j) ; 

		double Al = Al_n[jC] ;
		double Ze = Ze_n[jC] ;

		double P = P_n[jC] ;
		double Q = Q_n[jC] ;

		double r_Der_Ze = D1_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1], Ze_n[jC-1], dx) ;
		double r_Der_P  = D1_stereographic_center_2ndOrder(s_L, x_j,  P_n[jC+1],  P_n[jC-1], dx) ;

		double t_Der_Ze = D1_backward_2ndOrder(Ze_n[jC], Ze_nm1[jC], Ze_nm2[jC], dt) ;
		double t_Der_P  = D1_backward_2ndOrder( P_n[jC],  P_nm1[jC],  P_nm2[jC], dt) ;

		double phi_Der_f = 1. ;
		double phiphi_Der_f = 0. ;

		eom_TR[jC] = 
		- SE_LL_TR[jC] - (8*phiphi_Der_f*c_gbs*Al*P*Q*pow(Ze,2))/pow(r_j,2) - (8*phi_Der_f*t_Der_P*c_gbs*pow(Ze,3))/pow(r_j,2) - (Al*pow(Ze,3))/pow(r_j,2) - (8*phiphi_Der_f*c_gbs*Al*pow(P,2)*pow(Ze,3))/pow(r_j,2) + t_Der_Ze*((2*Ze)/r_j - (16*phi_Der_f*c_gbs*Q*Ze)/pow(r_j,2) - (16*phi_Der_f*c_gbs*P*pow(Ze,2))/pow(r_j,2)) + r_Der_Ze*((-2*Al*pow(Ze,2))/r_j + (8*phi_Der_f*c_gbs*Al*Q*pow(Ze,2))/pow(r_j,2) + (16*phi_Der_f*c_gbs*Al*P*pow(Ze,3))/pow(r_j,2)) + r_Der_P*((-8*phi_Der_f*c_gbs*Al*pow(Ze,2))/pow(r_j,2) + (8*phi_Der_f*c_gbs*Al*pow(Ze,4))/pow(r_j,2))
		;
	}
	return ;
}/*==========================================================================*/
static void compute_eom_ThTh(
	int Nx, 
	double dt, double dx,
	double s_L, double c_gbs,
	int exc_jC, 
	double x_lower,
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double*  P_n, double*  P_nm1, double*  P_nm2,
	double*  Q_n,
	double* SE_LL_ThTh,
	double* eom_ThTh)
{
	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		double x_j = x_lower + (jC * dx) ;
                double r_j = stereographic_r(s_L, x_j) ; 

		double Al = Al_n[jC] ;
		double Ze = Ze_n[jC] ;

		double P = P_n[jC] ;
		double Q = Q_n[jC] ;

		double r_Der_Al = D1_stereographic_center_2ndOrder(s_L, x_j, Al_n[jC+1], Al_n[jC-1], dx) ; 
		double r_Der_Ze = D1_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1], Ze_n[jC-1], dx) ;
		double r_Der_P  = D1_stereographic_center_2ndOrder(s_L, x_j,  P_n[jC+1],  P_n[jC-1], dx) ;
		double r_Der_Q  = D1_stereographic_center_2ndOrder(s_L, x_j,  Q_n[jC+1],  Q_n[jC-1], dx) ;

		double rr_Der_Al = D2_stereographic_center_2ndOrder(s_L, x_j, Al_n[jC+1], Al_n[jC], Al_n[jC-1], dx) ;
		double rr_Der_Ze = D2_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1], Ze_n[jC], Ze_n[jC-1], dx) ;

		double t_Der_Al = D1_backward_2ndOrder(Al_n[jC], Al_nm1[jC], Al_nm2[jC], dt) ;
		double t_Der_Ze = D1_backward_2ndOrder(Ze_n[jC], Ze_nm1[jC], Ze_nm2[jC], dt) ;
		double t_Der_P  = D1_backward_2ndOrder( P_n[jC],  P_nm1[jC],  P_nm2[jC], dt) ;

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
		double phi_Der_f = 1. ;
		double phiphi_Der_f = 0. ;
		eom_ThTh[jC] = 
		- SE_LL_ThTh[jC] + (r_j*tr_Der_Ze*(r_j - 8*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))/Al + (r_j*tr_Der_Al*Ze*(r_j - 8*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))/pow(Al,2) + r_j*rr_Der_Ze*Ze*(-r_j + 8*phi_Der_f*c_gbs*Q + 8*phi_Der_f*c_gbs*P*Ze) + (r_j*r_Der_Al*t_Der_Al*Ze*(-r_j + 8*phi_Der_f*c_gbs*Q + 8*phi_Der_f*c_gbs*P*Ze))/pow(Al,3) + r_j*pow(r_Der_Ze,2)*(-r_j + 8*phi_Der_f*c_gbs*Q + 16*phi_Der_f*c_gbs*P*Ze) - (8*phi_Der_f*r_j*pow(r_Der_Al,2)*c_gbs*Ze*(P + Q*Ze))/pow(Al,2) + (r_j*rr_Der_Al*(-r_j + 8*phi_Der_f*c_gbs*Q + 8*phi_Der_f*c_gbs*P*Ze)*(-1 + pow(Ze,2)))/Al + r_Der_Ze*(8*phi_Der_f*r_j*r_Der_Q*c_gbs*Ze - 2*r_j*(1 + 4*phiphi_Der_f*c_gbs*pow(P,2) - 4*phiphi_Der_f*c_gbs*pow(Q,2))*Ze + 8*phi_Der_f*r_j*r_Der_P*c_gbs*pow(Ze,2)) + t_Der_P*((-8*phi_Der_f*r_j*r_Der_Ze*c_gbs*Ze)/Al - (8*phi_Der_f*r_j*r_Der_Al*c_gbs*pow(Ze,2))/pow(Al,2)) + t_Der_Ze*((-8*phi_Der_f*r_j*r_Der_Q*c_gbs)/Al - (8*phi_Der_f*r_j*r_Der_Ze*c_gbs*P)/Al + (r_j - 8*phiphi_Der_f*r_j*c_gbs*pow(Q,2))/Al + (r_j*r_Der_Al*(r_j - 8*phi_Der_f*c_gbs*Q - 16*phi_Der_f*c_gbs*P*Ze))/pow(Al,2)) + r_Der_Al*((-8*phi_Der_f*r_j*r_Der_Q*c_gbs)/Al + (8*phi_Der_f*r_j*r_Der_P*c_gbs*Ze*(-2 + pow(Ze,2)))/Al + (r_j*r_Der_Ze*(-8*phi_Der_f*c_gbs*P + (-3*r_j + 16*phi_Der_f*c_gbs*Q)*Ze + 32*phi_Der_f*c_gbs*P*pow(Ze,2)))/Al - (r_j*(-1 + 8*phiphi_Der_f*c_gbs*pow(Q,2) + 16*phiphi_Der_f*c_gbs*P*Q*Ze + (1 + 8*phiphi_Der_f*c_gbs*pow(P,2))*pow(Ze,2)))/Al)
		;
		eom_ThTh[jC] /= pow(0.1+r_j,2) ; 
	}
	return ;
}
/*==========================================================================*/
static void compute_eom_scalar(
	int Nx, 
	double dt, double dx,
	double s_L, double c_gbs,
	int exc_jC, 
	double x_lower,
	double* Al_n, 
	double* Ze_n,
	double* P_n, double* P_nm1, double* P_nm2,
	double* Q_n,
	double* SE_LL_TR_n, 
	double* SE_LL_ThTh_n, 
	double* eom_scalar)
{
	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		double x_j = x_lower + (jC * dx) ;
                double r_j = stereographic_r(s_L, x_j) ; 

		double Al = Al_n[jC] ;
		double Ze = Ze_n[jC] ;

		double P = P_n[jC] ;
		double Q = Q_n[jC] ;

		double r_Der_Al = D1_stereographic_center_2ndOrder(s_L, x_j, Al_n[jC+1], Al_n[jC-1], dx) ; 
		double r_Der_Ze = D1_stereographic_center_2ndOrder(s_L, x_j, Ze_n[jC+1], Ze_n[jC-1], dx) ;
		double r_Der_P  = D1_stereographic_center_2ndOrder(s_L, x_j,  P_n[jC+1],  P_n[jC-1], dx) ;
		double r_Der_Q  = D1_stereographic_center_2ndOrder(s_L, x_j,  Q_n[jC+1],  Q_n[jC-1], dx) ;

		double t_Der_P  = D1_backward_2ndOrder( P_n[jC],  P_nm1[jC],  P_nm2[jC], dt) ;

		;
		double SE_LL_TR = SE_LL_TR_n[jC] ;
		double SE_LL_ThTh = SE_LL_ThTh_n[jC] ;

		double phi_Der_f = 1. ;
		double phiphi_Der_f = 0. ;
		eom_scalar[jC] = 
			(-2*Al*Q)/r_j - (2*Al*P*Ze)/r_j + SE_LL_TR*((-4*phi_Der_f*c_gbs*Ze)/pow(r_j,2) + (32*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*pow(Q,2)*Ze)/pow(r_j,2)) + pow(phi_Der_f,2)*((-128*pow(c_gbs,2)*Al*pow(Q,3))/pow(r_j,3) - (384*pow(c_gbs,2)*Al*P*pow(Q,2)*Ze)/pow(r_j,3) - (384*pow(c_gbs,2)*Al*pow(P,2)*Q*pow(Ze,2))/pow(r_j,3) - (128*pow(c_gbs,2)*Al*pow(P,3)*pow(Ze,3))/pow(r_j,3)) + phi_Der_f*((32*c_gbs*Al*pow(Q,2))/pow(r_j,2) + (64*c_gbs*Al*P*Q*Ze)/pow(r_j,2) + (32*c_gbs*Al*pow(P,2)*pow(Ze,2))/pow(r_j,2) - (4*c_gbs*Al*pow(Ze,4))/pow(r_j,4)) + phiphi_Der_f*phi_Der_f*((-32*pow(c_gbs,2)*Al*P*Q*pow(Ze,3))/pow(r_j,4) - (32*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,4))/pow(r_j,4) + (32*pow(c_gbs,2)*Al*pow(Q,2)*pow(Ze,4))/pow(r_j,4)) + pow(phiphi_Der_f,2)*phi_Der_f*((256*pow(c_gbs,3)*Al*P*pow(Q,3)*pow(Ze,3))/pow(r_j,4) + (256*pow(c_gbs,3)*Al*pow(P,2)*pow(Q,2)*pow(Ze,4))/pow(r_j,4)) + SE_LL_ThTh*((8*phi_Der_f*c_gbs*Al*pow(Ze,2))/pow(r_j,4) + pow(phi_Der_f,2)*((-64*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,5) - (64*pow(c_gbs,2)*Al*P*pow(Ze,3))/pow(r_j,5))) + pow(r_Der_Ze,2)*((64*pow(phi_Der_f,2)*pow(c_gbs,2)*Al*Q*pow(Ze,2))/pow(r_j,3) + pow(phi_Der_f,3)*((-512*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,2))/pow(r_j,4) - (256*pow(c_gbs,3)*Al*P*Q*pow(Ze,3))/pow(r_j,4))) + r_Der_Q*(-Al + (32*pow(phi_Der_f,2)*pow(c_gbs,2)*SE_LL_TR*Ze)/pow(r_j,2) + phi_Der_f*((16*c_gbs*Al*Q)/r_j + (16*c_gbs*Al*P*Ze)/r_j) + pow(phi_Der_f,2)*((-64*pow(c_gbs,2)*Al*pow(Q,2))/pow(r_j,2) - (128*pow(c_gbs,2)*Al*P*Q*Ze)/pow(r_j,2) - (64*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,2))/pow(r_j,2) + (32*pow(c_gbs,2)*Al*pow(Ze,4))/pow(r_j,4)) + phiphi_Der_f*pow(phi_Der_f,2)*((256*pow(c_gbs,3)*Al*P*Q*pow(Ze,3))/pow(r_j,4) + (256*pow(c_gbs,3)*Al*pow(P,2)*pow(Ze,4))/pow(r_j,4))) + pow(r_Der_Al,2)*(pow(phi_Der_f,2)*((64*pow(c_gbs,2)*P*pow(Ze,3))/(pow(r_j,3)*Al) + (64*pow(c_gbs,2)*Q*pow(Ze,4))/(pow(r_j,3)*Al)) + pow(phi_Der_f,3)*((-512*pow(c_gbs,3)*P*Q*pow(Ze,3))/(pow(r_j,4)*Al) - (512*pow(c_gbs,3)*pow(P,2)*pow(Ze,4))/(pow(r_j,4)*Al) - (512*pow(c_gbs,3)*pow(Q,2)*pow(Ze,4))/(pow(r_j,4)*Al) - (512*pow(c_gbs,3)*P*Q*pow(Ze,5))/(pow(r_j,4)*Al))) + r_Der_P*(-(Al*Ze) + phi_Der_f*((16*c_gbs*Al*Q*Ze)/r_j + (16*c_gbs*Al*P*pow(Ze,2))/r_j) + pow(phi_Der_f,2)*((-64*pow(c_gbs,2)*Al*pow(Q,2)*Ze)/pow(r_j,2) - (128*pow(c_gbs,2)*Al*P*Q*pow(Ze,2))/pow(r_j,2) - (32*pow(c_gbs,2)*Al*pow(Ze,3))/pow(r_j,4) - (64*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,3))/pow(r_j,2) + (32*pow(c_gbs,2)*Al*pow(Ze,5))/pow(r_j,4)) + pow(phi_Der_f,3)*r_Der_Q*((256*pow(c_gbs,3)*Al*pow(Ze,3))/pow(r_j,4) - (256*pow(c_gbs,3)*Al*pow(Ze,5))/pow(r_j,4)) + phiphi_Der_f*pow(phi_Der_f,2)*((256*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,3))/pow(r_j,4) - (256*pow(c_gbs,3)*Al*pow(Q,2)*pow(Ze,5))/pow(r_j,4))) + t_Der_P*(1 + (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,4))/pow(r_j,4) + (256*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(Q,2)*pow(Ze,4))/pow(r_j,4) + phi_Der_f*((-16*c_gbs*Q)/r_j - (16*c_gbs*P*Ze)/r_j) + pow(phi_Der_f,2)*((64*pow(c_gbs,2)*pow(Q,2))/pow(r_j,2) + (128*pow(c_gbs,2)*P*Q*Ze)/pow(r_j,2) + (64*pow(c_gbs,2)*pow(P,2)*pow(Ze,2))/pow(r_j,2) - (32*pow(c_gbs,2)*pow(Ze,4))/pow(r_j,4)) + r_Der_Ze*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,3) + pow(phi_Der_f,3)*((-1024*pow(c_gbs,3)*Q*pow(Ze,3))/pow(r_j,4) - (768*pow(c_gbs,3)*P*pow(Ze,4))/pow(r_j,4))) + r_Der_Al*((128*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))/(pow(r_j,3)*Al) + pow(phi_Der_f,3)*((-1024*pow(c_gbs,3)*Q*pow(Ze,4))/(pow(r_j,4)*Al) - (768*pow(c_gbs,3)*P*pow(Ze,5))/(pow(r_j,4)*Al)))) + r_Der_Ze*(-(Al*P) + (256*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*Al*Q*pow(Ze,3))/pow(r_j,4) + phi_Der_f*((16*c_gbs*Al*P*Q)/r_j + (16*c_gbs*Al*pow(P,2)*Ze)/r_j + (16*c_gbs*Al*pow(Ze,3))/pow(r_j,3)) + pow(phi_Der_f,2)*((-64*pow(c_gbs,2)*Al*P*pow(Q,2))/pow(r_j,2) - (128*pow(c_gbs,2)*Al*pow(P,2)*Q*Ze)/pow(r_j,2) - (64*pow(c_gbs,2)*Al*pow(P,3)*pow(Ze,2))/pow(r_j,2) - (160*pow(c_gbs,2)*Al*Q*pow(Ze,3))/pow(r_j,4) - (96*pow(c_gbs,2)*Al*P*pow(Ze,4))/pow(r_j,4)) + SE_LL_TR*((8*phi_Der_f*c_gbs)/r_j + pow(phi_Der_f,2)*((-64*pow(c_gbs,2)*Q)/pow(r_j,2) - (32*pow(c_gbs,2)*P*Ze)/pow(r_j,2))) + phiphi_Der_f*(phi_Der_f*((64*pow(c_gbs,2)*Al*P*Q*pow(Ze,2))/pow(r_j,3) + (128*pow(c_gbs,2)*Al*pow(P,2)*pow(Ze,3))/pow(r_j,3)) + pow(phi_Der_f,2)*((-512*pow(c_gbs,3)*Al*P*pow(Q,2)*pow(Ze,2))/pow(r_j,4) - (1280*pow(c_gbs,3)*Al*pow(P,2)*Q*pow(Ze,3))/pow(r_j,4) + (256*pow(c_gbs,3)*Al*pow(Q,3)*pow(Ze,3))/pow(r_j,4) - (768*pow(c_gbs,3)*Al*pow(P,3)*pow(Ze,4))/pow(r_j,4))) + r_Der_P*(pow(phi_Der_f,2)*((64*pow(c_gbs,2)*Al*pow(Ze,2))/pow(r_j,3) - (128*pow(c_gbs,2)*Al*pow(Ze,4))/pow(r_j,3)) + pow(phi_Der_f,3)*((-512*pow(c_gbs,3)*Al*Q*pow(Ze,2))/pow(r_j,4) - (256*pow(c_gbs,3)*Al*P*pow(Ze,3))/pow(r_j,4) + (1024*pow(c_gbs,3)*Al*Q*pow(Ze,4))/pow(r_j,4) + (768*pow(c_gbs,3)*Al*P*pow(Ze,5))/pow(r_j,4)))) + r_Der_Al*(-Q - P*Ze + phi_Der_f*((16*c_gbs*pow(Q,2))/r_j + (32*c_gbs*P*Q*Ze)/r_j - (8*c_gbs*pow(Ze,2))/pow(r_j,3) + (16*c_gbs*pow(P,2)*pow(Ze,2))/r_j + (16*c_gbs*pow(Ze,4))/pow(r_j,3)) + pow(phi_Der_f,2)*((-64*pow(c_gbs,2)*pow(Q,3))/pow(r_j,2) - (192*pow(c_gbs,2)*P*pow(Q,2)*Ze)/pow(r_j,2) + (64*pow(c_gbs,2)*Q*pow(Ze,2))/pow(r_j,4) - (192*pow(c_gbs,2)*pow(P,2)*Q*pow(Ze,2))/pow(r_j,2) + (64*pow(c_gbs,2)*P*pow(Ze,3))/pow(r_j,4) - (64*pow(c_gbs,2)*pow(P,3)*pow(Ze,3))/pow(r_j,2) - (128*pow(c_gbs,2)*Q*pow(Ze,4))/pow(r_j,4) - (96*pow(c_gbs,2)*P*pow(Ze,5))/pow(r_j,4)) + SE_LL_TR*((8*phi_Der_f*c_gbs*Ze)/(r_j*Al) + pow(phi_Der_f,2)*((-64*pow(c_gbs,2)*Q*Ze)/(pow(r_j,2)*Al) - (32*pow(c_gbs,2)*P*pow(Ze,2))/(pow(r_j,2)*Al))) + r_Der_Q*((64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,2))/pow(r_j,3) + pow(phi_Der_f,3)*((-512*pow(c_gbs,3)*Q*pow(Ze,2))/pow(r_j,4) - (512*pow(c_gbs,3)*P*pow(Ze,3))/pow(r_j,4))) + r_Der_Ze*((16*phi_Der_f*c_gbs*Ze)/pow(r_j,2) + pow(phi_Der_f,2)*((-256*pow(c_gbs,2)*Q*Ze)/pow(r_j,3) - (192*pow(c_gbs,2)*P*pow(Ze,2))/pow(r_j,3) + (128*pow(c_gbs,2)*Q*pow(Ze,3))/pow(r_j,3)) + pow(phi_Der_f,3)*((1024*pow(c_gbs,3)*pow(Q,2)*Ze)/pow(r_j,4) + (1536*pow(c_gbs,3)*P*Q*pow(Ze,2))/pow(r_j,4) + (512*pow(c_gbs,3)*pow(P,2)*pow(Ze,3))/pow(r_j,4) - (1024*pow(c_gbs,3)*pow(Q,2)*pow(Ze,3))/pow(r_j,4) - (768*pow(c_gbs,3)*P*Q*pow(Ze,4))/pow(r_j,4))) + phiphi_Der_f*(phi_Der_f*((64*pow(c_gbs,2)*pow(Q,2)*pow(Ze,2))/pow(r_j,3) + (192*pow(c_gbs,2)*P*Q*pow(Ze,3))/pow(r_j,3) + (128*pow(c_gbs,2)*pow(P,2)*pow(Ze,4))/pow(r_j,3)) + pow(phi_Der_f,2)*((-512*pow(c_gbs,3)*pow(Q,3)*pow(Ze,2))/pow(r_j,4) - (2048*pow(c_gbs,3)*P*pow(Q,2)*pow(Ze,3))/pow(r_j,4) - (2304*pow(c_gbs,3)*pow(P,2)*Q*pow(Ze,4))/pow(r_j,4) - (768*pow(c_gbs,3)*pow(P,3)*pow(Ze,5))/pow(r_j,4))) + r_Der_P*(pow(phi_Der_f,2)*((192*pow(c_gbs,2)*pow(Ze,3))/pow(r_j,3) - (128*pow(c_gbs,2)*pow(Ze,5))/pow(r_j,3)) + pow(phi_Der_f,3)*((-1536*pow(c_gbs,3)*Q*pow(Ze,3))/pow(r_j,4) - (1280*pow(c_gbs,3)*P*pow(Ze,4))/pow(r_j,4) + (1024*pow(c_gbs,3)*Q*pow(Ze,5))/pow(r_j,4) + (768*pow(c_gbs,3)*P*pow(Ze,6))/pow(r_j,4))))
		;
	}
	return ;
}

/*===========================================================================*/
/* compute characteristics */
/*===========================================================================*/
static double compute_element_P_eom_tDer_P(
	double c_gbs,     double r_j, 
	double phi_Der_f, double phiphi_Der_f, 
	double rho,       
	double Ze,       double Q,       double P,
	double r_Der_Q)	
{
	return 
		(768*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*pow(Ze,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))/(pow(r_j,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 12*phi_Der_f*c_gbs*P*Ze)) + (pow(r_j,5) + 64*pow(phi_Der_f,2)*rho*pow(r_j,3)*pow(c_gbs,2)*pow(Ze,2) + 256*pow(phi_Der_f,2)*pow(r_j,3)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2) - 768*pow(phi_Der_f,3)*pow(r_j,2)*pow(c_gbs,3)*pow(P,3)*pow(Ze,3) - 96*pow(phi_Der_f,2)*r_j*pow(c_gbs,2)*pow(Ze,4) - 512*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Q,3)*(pow(r_j,2) + 12*phiphi_Der_f*c_gbs*pow(Ze,4)) - 8*phi_Der_f*c_gbs*Q*(3*pow(r_j,4) - 56*phi_Der_f*pow(r_j,3)*c_gbs*P*Ze + 64*pow(phi_Der_f,2)*rho*pow(r_j,2)*pow(c_gbs,2)*pow(Ze,2) + 256*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2) - 96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4)) - 4*P*(7*phi_Der_f*pow(r_j,4)*c_gbs*Ze + 96*pow(phi_Der_f,3)*rho*pow(r_j,2)*pow(c_gbs,3)*pow(Ze,3) - 192*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Ze,5)) + 64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*(-4*phi_Der_f*c_gbs*P*Ze*(7*pow(r_j,2) + 24*phiphi_Der_f*c_gbs*pow(Ze,4)) + 3*(pow(r_j,3) + 4*phiphi_Der_f*r_j*c_gbs*pow(Ze,4))))/(pow(r_j,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 12*phi_Der_f*c_gbs*P*Ze))
	; 
}
/*===========================================================================*/
static double compute_element_P_eom_tDer_Q(void)
{
	return 0 ; 
}
/*===========================================================================*/
static double compute_element_Q_eom_tDer_P(void)
{
	return 0 ; 
}
/*===========================================================================*/
static double compute_element_Q_eom_tDer_Q(void)
{
	return 1. ; 
}
/*===========================================================================*/
static double compute_element_P_eom_rDer_P(
	double c_gbs,     double r_j, 
	double phi_Der_f, double phiphi_Der_f, 
	double rho,       double Jr,       
	double Al,        double Ze,       double Q,       double P,
	double r_Der_P,   double r_Der_Q)	
{
	return
		(-1536*pow(phi_Der_f,3)*r_Der_P*pow(c_gbs,3)*Al*pow(Ze,4))/pow(r_j,4) - (768*pow(phi_Der_f,3)*r_Der_Q*pow(c_gbs,3)*Al*pow(Ze,5)*(r_j - 4*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))/(pow(r_j,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 12*phi_Der_f*c_gbs*P*Ze)) - (Al*Ze*(6144*pow(phi_Der_f,4)*pow(c_gbs,4)*pow(Q,4)*(pow(r_j,2) + 4*phiphi_Der_f*c_gbs*pow(Ze,4)) + 32*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*(15*pow(r_j,4) - 256*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*Jr*Ze + 24*phiphi_Der_f*pow(r_j,2)*c_gbs*pow(Ze,4) - 96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4) - 16*phi_Der_f*r_j*c_gbs*P*Ze*(17*pow(r_j,2) + 48*phiphi_Der_f*c_gbs*pow(Ze,2) + 24*phiphi_Der_f*c_gbs*pow(Ze,4)) + 64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2)*(19*pow(r_j,2) + 120*phiphi_Der_f*c_gbs*pow(Ze,2) + 24*phiphi_Der_f*c_gbs*pow(Ze,4))) - 256*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Q,3)*(11*pow(r_j,3) + 36*phiphi_Der_f*r_j*c_gbs*pow(Ze,4) - 4*phi_Der_f*c_gbs*P*Ze*(25*pow(r_j,2) + 96*phiphi_Der_f*c_gbs*pow(Ze,2) + 72*phiphi_Der_f*c_gbs*pow(Ze,4))) + (r_j - 8*phi_Der_f*c_gbs*P*Ze)*(pow(r_j,5) - 128*pow(phi_Der_f,2)*pow(r_j,3)*pow(c_gbs,2)*Jr*Ze + 64*pow(phi_Der_f,2)*rho*pow(r_j,3)*pow(c_gbs,2)*pow(Ze,2) + 256*pow(phi_Der_f,2)*pow(r_j,3)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2) - 768*pow(phi_Der_f,3)*pow(r_j,2)*pow(c_gbs,3)*pow(P,3)*pow(Ze,3) - 96*pow(phi_Der_f,2)*r_j*pow(c_gbs,2)*pow(Ze,4) + 4*phi_Der_f*c_gbs*P*Ze*(-7*pow(r_j,4) + 384*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*Jr*Ze - 96*pow(phi_Der_f,2)*rho*pow(r_j,2)*pow(c_gbs,2)*pow(Ze,2) + 192*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4))) - 4*phi_Der_f*c_gbs*Q*(9*pow(r_j,5) - 512*pow(phi_Der_f,2)*pow(r_j,3)*pow(c_gbs,2)*Jr*Ze + 128*pow(phi_Der_f,2)*rho*pow(r_j,3)*pow(c_gbs,2)*pow(Ze,2) - 288*pow(phi_Der_f,2)*r_j*pow(c_gbs,2)*pow(Ze,4) + 128*pow(phi_Der_f,2)*r_j*pow(c_gbs,2)*pow(P,2)*pow(Ze,2)*(17*pow(r_j,2) + 60*phiphi_Der_f*c_gbs*pow(Ze,2)) - 256*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(P,3)*pow(Ze,3)*(25*pow(r_j,2) + 144*phiphi_Der_f*c_gbs*pow(Ze,2)) + 4*phi_Der_f*c_gbs*P*Ze*(-61*pow(r_j,4) + 1280*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*Jr*Ze - 32*pow(r_j,2)*c_gbs*(3*phiphi_Der_f + 5*pow(phi_Der_f,2)*rho*c_gbs)*pow(Ze,2) + 576*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4)))))/(pow(r_j,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 12*phi_Der_f*c_gbs*P*Ze)*(r_j - 8*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))
	; 
}
/*===========================================================================*/
static double compute_element_P_eom_rDer_Q(
	double c_gbs,     double r_j, 
	double phi_Der_f, double phiphi_Der_f, 
	double rho,       double Jr,       double SE_LL_TR,       
	double Al,        double Ze,       double Q,       double P,
	double r_Der_P,   double t_Der_P)	
{
	return
		(768*pow(phi_Der_f,3)*t_Der_P*pow(c_gbs,3)*pow(Ze,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))/(pow(r_j,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 12*phi_Der_f*c_gbs*P*Ze)) - (768*pow(phi_Der_f,3)*r_Der_P*pow(c_gbs,3)*Al*pow(Ze,5)*(r_j - 4*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))/(pow(r_j,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 12*phi_Der_f*c_gbs*P*Ze)) + (64*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*SE_LL_TR*Ze*pow(r_j - 8*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze,2) + Al*(1024*pow(phi_Der_f,4)*pow(r_j,2)*pow(c_gbs,4)*pow(Q,4)*(-4 + pow(Ze,2)) + 256*pow(phi_Der_f,3)*pow(c_gbs,3)*pow(Q,3)*(8*pow(r_j,3) - 64*phi_Der_f*pow(r_j,2)*c_gbs*P*Ze - pow(r_j,3)*pow(Ze,2) + 8*phi_Der_f*pow(r_j,2)*c_gbs*P*pow(Ze,3) + 96*phiphi_Der_f*phi_Der_f*pow(c_gbs,2)*P*pow(Ze,5)) + 16*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Q,2)*(-24*pow(r_j,4) + 128*phi_Der_f*pow(r_j,2)*c_gbs*(2*phi_Der_f*c_gbs*Jr + 3*r_j*P)*Ze + (pow(r_j,4) - 1536*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*pow(P,2))*pow(Ze,2) - 16*phi_Der_f*pow(r_j,2)*c_gbs*(4*phi_Der_f*c_gbs*Jr + r_j*P)*pow(Ze,3) + 32*pow(phi_Der_f,2)*pow(c_gbs,2)*(12 - rho*pow(r_j,2) + 2*(pow(r_j,2) + 48*phiphi_Der_f*c_gbs)*pow(P,2))*pow(Ze,4) - 192*phiphi_Der_f*phi_Der_f*r_j*pow(c_gbs,2)*P*pow(Ze,5) + 1536*phiphi_Der_f*pow(phi_Der_f,2)*pow(c_gbs,3)*pow(P,2)*pow(Ze,6)) - pow(r_j - 8*phi_Der_f*c_gbs*P*Ze,2)*(pow(r_j,4) - 64*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*Jr*Ze - 16*phi_Der_f*pow(r_j,3)*c_gbs*P*Ze - 96*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4) + 64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2)*(pow(r_j,2) - 12*phiphi_Der_f*c_gbs*pow(Ze,2))) + 32*phi_Der_f*c_gbs*Q*(r_j - 8*phi_Der_f*c_gbs*P*Ze)*(pow(r_j,4) - 16*phi_Der_f*pow(r_j,3)*c_gbs*P*Ze - 48*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(Ze,4) + 4*pow(phi_Der_f,2)*pow(r_j,2)*pow(c_gbs,2)*Jr*Ze*(-8 + pow(Ze,2)) + 64*pow(phi_Der_f,2)*pow(c_gbs,2)*pow(P,2)*pow(Ze,2)*(pow(r_j,2) - 6*phiphi_Der_f*c_gbs*pow(Ze,2)))))/(pow(r_j,4)*(r_j - 8*phi_Der_f*c_gbs*Q - 12*phi_Der_f*c_gbs*P*Ze)*(r_j - 8*phi_Der_f*c_gbs*Q - 8*phi_Der_f*c_gbs*P*Ze))
	; 
}
/*===========================================================================*/
static double compute_element_Q_eom_rDer_P(double Al)
{
	return 
		- Al 
	; 
}
/*===========================================================================*/
static double compute_element_Q_eom_rDer_Q(double Al, double Ze)
{
	return 
		- Al*Ze
	; 
}
/*===========================================================================*/
static void computeArray_center_EdGB_characteristics(
	int Nx, 
	double dt, double dx,
	double s_L, double c_gbs,
	int exc_jC, 
	double x_lower,
	double* Al_n, 
	double* Ze_n, 
	double*  P_n, double* P_nm1, double* P_nm2, 
	double*  Q_n,
	double* SE_LL_TR,
	int* elliptic_pts,
	double* ingoing, double* outgoing)
{
	double phi_Der_f = 1. ;
	double phiphi_Der_f = 0. ; 

	for (int jC=exc_jC+1; jC<Nx-1; jC++) {
		double x_j = x_lower + (dx * jC) ;

		double r_j = stereographic_r( s_L, x_j) ;
		double dr  = stereographic_dr(s_L, x_j, dx) ;

		double Al = Al_n[jC] ;
		double Ze = Ze_n[jC] ;
		double P  = P_n[jC] ;
		double Q  = Q_n[jC] ;

		double r_Der_P = D1_center_2ndOrder(P_n[jC+1], P_n[jC-1], dr) ;
		double r_Der_Q = D1_center_2ndOrder(Q_n[jC+1], Q_n[jC-1], dr) ;

		double t_Der_P = D1_backward_2ndOrder(P_n[jC], P_nm1[jC], P_nm2[jC], dt) ;

		double rho = 0 ;
		double Jr = 0 ;

		double matrix_A[2][2] = {0} ;
		matrix_A[0][0] = compute_element_P_eom_tDer_P(
			c_gbs,     r_j, 
			phi_Der_f, phiphi_Der_f, 
			rho,       
			Ze,       Q,       P,
			r_Der_Q)	
		;
		matrix_A[0][1] = compute_element_P_eom_tDer_Q()	;
		matrix_A[1][0] = compute_element_Q_eom_tDer_P()	;
		matrix_A[1][1] = compute_element_Q_eom_tDer_Q()	;

		double matrix_B[2][2] = {0} ;
		matrix_B[0][0] = compute_element_P_eom_rDer_P(
			c_gbs,     r_j, 
			phi_Der_f, phiphi_Der_f, 
			rho,       Jr,       
			Al,        Ze,       Q,       P,
			r_Der_P,   r_Der_Q)	
		;
		matrix_B[0][1] = compute_element_P_eom_rDer_Q(
			c_gbs,     r_j, 
			phi_Der_f, phiphi_Der_f, 
			rho,       Jr,       SE_LL_TR[jC],
			Al,        Ze,       Q,       P,
			r_Der_P,   t_Der_P)	
		;
		matrix_B[1][0] = compute_element_Q_eom_rDer_P(Al) ;
		matrix_B[1][1] = compute_element_Q_eom_rDer_Q(Al, Ze) ;

		double matrix_inv_A[2][2] = {0} ;
		compute_2dMatrix_inverse(matrix_A, matrix_inv_A) ;
		
		double matrix_inv_AB[2][2] = {0} ;
		compute_2dMatrix_multiply(matrix_inv_A, matrix_B, matrix_inv_AB) ;

		double eigenvalues[2] = {0} ;
		bool is_elliptic = 
			compute_2dMatrix_eigenvalues(matrix_inv_AB, eigenvalues) 
		;
		if (is_elliptic == true) {
			*elliptic_pts += 1 ;
		}
		outgoing[jC] =
			eigenvalues[0] 
		;
		ingoing[jC] =
			eigenvalues[1] 
		;
	}
	return ;
}
/*==========================================================================*/
void compute_diagnostics_massless_scalar_EdGB(
	int Nx, int exc_jC, 
	double s_L,double c_gbs,
	double dt, double dx,
	double x_lower,
	double* Al_n, double* Al_nm1, double* Al_nm2,
	double* Ze_n, double* Ze_nm1, double* Ze_nm2,
	double* P_n, double* P_nm1, double* P_nm2,
	double* Q_n, 
	double* eom_TR,
	double* eom_ThTh,
	double* eom_scalar,
	double* ingoing_characteristic,
	double* outgoing_characteristic)
{
	double* SE_LL_TR   = calloc(Nx, sizeof(double)) ;
	double* SE_LL_ThTh = calloc(Nx, sizeof(double)) ;

	compute_SE_LL_components_massless_scalar(
		Nx, dx,
		s_L, exc_jC,
		x_lower,
		Al_n, Ze_n,
		P_n,  Q_n,
		SE_LL_TR,
		SE_LL_ThTh)
	;
	compute_eom_TR(
		Nx, 
		dt, dx,
		s_L, c_gbs,
		exc_jC, 
		x_lower,
		Al_n, 
		Ze_n, Ze_nm1, Ze_nm2,
		 P_n,  P_nm1,  P_nm2,
		 Q_n,
		SE_LL_TR,
		eom_TR)
	;
	compute_eom_ThTh(
		Nx, 
		dt, dx,
		s_L, c_gbs,
		exc_jC, 
		x_lower,
		Al_n, Al_nm1, Al_nm2,
		Ze_n, Ze_nm1, Ze_nm2,
		 P_n,  P_nm1,  P_nm2,
		 Q_n,
		SE_LL_ThTh,
		eom_ThTh)
	;
	compute_eom_scalar(
		Nx, 
		dt, dx,
		s_L, c_gbs,
		exc_jC, 
		x_lower,
		Al_n, 
		Ze_n,
		P_n, P_nm1, P_nm2,
		Q_n, 
		SE_LL_TR,
		SE_LL_ThTh,
		eom_scalar)
	;
	int elliptic_pts = 0 ;
	computeArray_center_EdGB_characteristics(
		Nx, 
		dt, dx,
		s_L, c_gbs,
		exc_jC, 
		x_lower,
		Al_n,
		Ze_n, 
		P_n,   P_nm1,  P_nm2,
		Q_n,
		SE_LL_TR,
		&elliptic_pts,
		ingoing_characteristic, 
		outgoing_characteristic)
	;
	set_array_val(0, exc_jC, 0, P_n) ;
	set_array_val(0, exc_jC, 0, Q_n) ;
	set_array_val(0, exc_jC, 0, eom_TR) ;
	set_array_val(0, exc_jC, 0, eom_ThTh) ;
	set_array_val(0, exc_jC, 0, eom_scalar) ;
	set_array_val(0, exc_jC, 0, ingoing_characteristic) ;
	set_array_val(0, exc_jC, 0, outgoing_characteristic) ;

	free(SE_LL_TR) ;
	free(SE_LL_ThTh) ;

	return ;
}
