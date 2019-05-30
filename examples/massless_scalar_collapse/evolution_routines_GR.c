#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "evolution_routines_GR.h"
#include "stereographic_routines.h"

#define ERR_TOLERANCE ((double)1e-10)
#define MACHINE_EPSILON ((double)1e-14)

/*==========================================================================*/
/* 	Notations:
	s_L: stereographic length
	exc_jC: excision grid point (index jC)
*/
/*==========================================================================*/

/*==========================================================================*/
inline double weighted_infty_norm(double weight, double val_1, double val_2)
{
        return (fabs(val_2)>fabs(val_1)) ? fabs(weight*val_2) : fabs(weight*val_1) ;
}
/*==========================================================================*/
static inline double max_fabs(double var_1, double var_2) 
{
	return (fabs(var_1)>fabs(var_2)) ? fabs(var_1) : fabs(var_2) ;
}
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
/* see gr-qc/0302072 */
/*==========================================================================*/
void Kreiss_Oliger_Filter(
	int Nx,
	double* field)
{
	double epsilon_ko = 0.4 ;
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
/* for outer excision boundary */
	epsilon_ko = 1.0 ;
	field[1] += (epsilon_ko/16.) * (
			field[4] 
		+ 	(-4.*field[3]) 
		+ 	(6.*field[2]) 
		+ 	(-4.*field[1]) 
		+ 	field[0] 
		) ;
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
/* ODE solvers for lapse and shift */
/*==========================================================================*/
static double compute_iteration_GR_Al(
	double s_L,
	int Nx,
	double dt, 	double dx,
	int start_jC,
	double bbox[2],
	double* Al, 	double* Ze, 
	double*  P, 	double*  Q)
{
        double
                x_j1h, r_j1h, dr
        ;
        double
                Al_j1h, Ze_j1h, Q_j1h, P_j1h
        ;
        double
                r_Der_Al_j1h
        ;
        double
                Jr_j1h
        ;
	int size = 0 ;
	if (fabs(bbox[1]-s_L)<MACHINE_EPSILON) size = Nx-1 ;
	else size = Nx ;
        double res_Al = 0 ;
        double jac_Al = 1 ;
        double res_infty_norm = 0 ; /* returning this */
 /* scalar field functions */
        for (int jC=start_jC; jC<Nx-1; jC++) {
                x_j1h = ((jC+1) + jC) * dx / 2 ;

                r_j1h = stereographic_r(s_L, x_j1h) ;

                dr = stereographic_dr(s_L, x_j1h, dx) ;
/*--------------------------------------------------------------------------*/
/* there is Jr/Ze term that is zero in empty
 * space but undefined if below machine precision.
 * The slope is zero though in empty spacce so
 * to that precision we do this. */
/*--------------------------------------------------------------------------*/
                Al_j1h = (Al[jC+1] + Al[jC]) / 2 ;
                Ze_j1h = (Ze[jC+1] + Ze[jC]) / 2 ;

                Q_j1h = (Q[jC+1] + Q[jC]) / 2 ;
                P_j1h = (P[jC+1] + P[jC]) / 2 ;

                r_Der_Al_j1h = (Al[jC+1] - Al[jC]) / dr ;

                Jr_j1h = - Q_j1h * P_j1h ;

                if ((fabs(Jr_j1h) < 10*MACHINE_EPSILON)
                &&  (fabs(Ze_j1h) < 10*MACHINE_EPSILON)
                ) {
                        Al[jC+1] = Al[jC] ;
                } else {
                        res_Al =
                        +       r_Der_Al_j1h*Ze_j1h
                        -       (r_j1h*Al_j1h*Jr_j1h)/2.
                        ;
                        jac_Al =
                                Ze_j1h/dr
                        -       r_j1h*Jr_j1h/4.
                        ;
                        Al[jC+1] -= res_Al / jac_Al ;
                }
/*--------------------------------------------------------------------------*/
                if ((isnan(res_Al) != 0)
                ||  (isnan(jac_Al) != 0)
                ) {
                        printf("jC:%d\tres_Al:%.e\tjac_Al:%.e\n", jC, res_Al, jac_Al) ;
                        exit(EXIT_FAILURE) ;
                }
                res_infty_norm = weighted_infty_norm(1-x_j1h/s_L, res_Al, res_infty_norm) ;
        }
	if (size==Nx-1) Al[Nx-1] = Al[Nx-2] ;

        return res_infty_norm ;
}
/*==========================================================================*/
static double compute_iteration_GR_Ze(
	double s_L,
	int Nx,
	double dt, 	double dx,
	int start_jC,
	double bbox[2],
	double* Al, 	double* Ze, 
	double*  P, 	double*  Q)
{
        double
                x_j1h, x_jp1, x_j,
                r_j1h, r_jp1, r_j,
                dr
        ;
        double
                Al_sqrd_jp1, Al_sqrd_j,
                Ze_sqrd_jp1, Ze_sqrd_j
        ;
        double
                Al_j1h, Q_j1h, P_j1h
        ;
        double
                rho_j1h
        ;
	int size = 0 ;
	if (fabs(bbox[1]-s_L)<MACHINE_EPSILON) size = Nx-1 ;
	else size = Nx ;
        double res_Ze_sqrd = 0 ; 
        double jac_Ze_sqrd = 1 ; 
 /* scalar field functions */   

        double res_infty_norm = 0 ; /* returning this */

        for (int jC=start_jC; jC<size-1; jC++) {
                x_j1h = ((jC+1) + jC) * dx / 2 ; 

                x_jp1 = (jC+1) * dx ;
                x_j   = (jC+0) * dx ;

                r_j1h = stereographic_r(s_L, x_j1h) ;
                r_jp1 = stereographic_r(s_L, x_jp1) ;
                r_j   = stereographic_r(s_L, x_j  ) ; 

                dr = stereographic_dr(s_L, x_j1h, dx) ;

                Al_j1h = (Al[jC+1] + Al[jC]) / 2. ;

                Al_sqrd_jp1 = pow(Al[jC+1], 2) ;
                Al_sqrd_j   = pow(Al[jC+0], 2) ;

                Ze_sqrd_jp1 = pow(Ze[jC+1], 2) ;
                Ze_sqrd_j   = pow(Ze[jC+0], 2) ;

                Q_j1h = (Q[jC+1] + Q[jC]) / 2. ;
                P_j1h = (P[jC+1] + P[jC]) / 2. ;

                rho_j1h = (1./2) * (pow(Q_j1h,2) + pow(P_j1h,2)) ;
/*---------------------------------------------------------------------------*/
                res_Ze_sqrd = 
                (
                        (       (r_jp1)*Al_sqrd_jp1*Ze_sqrd_jp1
                        -       (r_j  )*Al_sqrd_j  *Ze_sqrd_j
                        )/dr
                )
                -       pow(r_j1h,2)*pow(Al_j1h,2)*rho_j1h
                ;

                jac_Ze_sqrd = 
                (
                        (
                                (r_jp1)*Al_sqrd_jp1
                        )/dr
                )
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
                res_infty_norm = weighted_infty_norm(1-x_j1h/s_L, res_Ze_sqrd, res_infty_norm) ;
        }
	if (size==Nx-1) Ze[Nx-1] = 0 ;

        return res_infty_norm ;
}
/*==========================================================================*/
static double compute_iteration_GR_excision_boundary_condition_Ze(
	double s_L,
	int Nx,
	double dt, 	double dx,
	int exc_jC,
	double bbox[2],
	double* Al_n,  double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double* P_n,   double* P_nm1,  double* Q_n,  double* Q_nm1)
{
        double x_j = dx * exc_jC ;
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
void solve_Al_Ze(
	double s_L,
	int Nx,
	double dt, 	double dx,
	bool excision_on,
	int exc_jC,
	int child_perim_coords[2],
	double bbox[2],
	double* Al_n, 	double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double*  P_n, 	double*  P_nm1, double*  Q_n, double*  Q_nm1)
{
	if (exc_jC==Nx-1) {
		return ;
	}
	double res = 0 ;
	int start_jC = 0 ;
	do {
		res = 0 ;
		if ((excision_on==true)
		&&  (exc_jC>0)
		&&  	(   (child_perim_coords[1]==Nx-1) /* i.e. if this is the shadow grid */
		     	||  (exc_jC>child_perim_coords[1])
			)
		) {	
			start_jC = exc_jC ;
			res += compute_iteration_GR_excision_boundary_condition_Ze(
				s_L,
				Nx,
				dt, 	dx,
				start_jC,
				bbox,
				Al_n,  Al_nm1, Ze_n, Ze_nm1,
				 P_n,   P_nm1,  Q_n,  Q_nm1)
			;
		} else if (
		    (child_perim_coords[1]>0)
		&&  (child_perim_coords[1]>exc_jC)
		&&  (child_perim_coords[1]!=Nx-1) /* i.e. if this is not the shadow grid */
		) {
			start_jC = child_perim_coords[1] ;
		} else {
			start_jC = 0 ;
		}
		res += compute_iteration_GR_Ze(
			s_L,
			Nx,
			dt, 	dx,
			start_jC,
			bbox,
			Al_n, 	Ze_n, 
			P_n, 	Q_n)
		;
		res += compute_iteration_GR_Al(
			s_L,
			Nx,
			dt, 	dx,
			start_jC,
			bbox,
			Al_n, 	Ze_n, 
			P_n, 	Q_n)
		;
	} while (res>ERR_TOLERANCE) ;
	return ;
}
/*==========================================================================*/
static double compute_iteration_GR_Crank_Nicolson_PQ(
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

	double  
		x_jm1, x_j, x_jp1, x_jp2, 
		r_jm1, r_j, r_jp1, r_jp2, dr 
	;
	double res_Q = 0 ;
	double res_P = 0 ;
	double jac_Q = 1 ;
	double jac_P = 1 ;

	int size = 0 ;
	if (fabs(bbox[1]-s_L)<MACHINE_EPSILON) size = Nx-1 ;
	else size = Nx ;
	double res_infty_norm = 0 ; /* returning this */
/*--------------------------------------------------------------------------*/
/* interior: we go to Nx-2 as we do not want to actually include the point
   at infinity in our computational domain */
/*--------------------------------------------------------------------------*/
	for (int jC=exc_jC+1;jC<size-1;jC++) {
		x_j   = lower_x + (dx * (jC)  ) ;
		x_jp1 = lower_x + (dx * (jC+1)) ;
		x_jm1 = lower_x + (dx * (jC-1)) ;

		r_j   = stereographic_r(s_L, x_j  ) ;
		r_jp1 = stereographic_r(s_L, x_jp1) ;
		r_jm1 = stereographic_r(s_L, x_jm1) ;

		dr = stereographic_dr(s_L, x_j, dx) ;
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
	/* one iteration */
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
		x_j = 0 ;
		dr  = pow(1. - (x_j/s_L), -2) * dx;

		res_P = D1_forward_2ndOrder(
			P_n[2], P_n[1], P_n[0], 
		dr) ;
		jac_P = - (3./2.) / dr ;

		P_n[0] -= res_P / jac_P ;
	}
	if ((exc_jC > 0) 
	&&  (excision_on==true)
	) {
		x_j   = lower_x + (dx * (exc_jC))   ;
		x_jp1 = lower_x + (dx * (exc_jC+1)) ;
		x_jp2 = lower_x + (dx * (exc_jC+2)) ;

		r_j   = stereographic_r(s_L, x_j  ) ;
		r_jp1 = stereographic_r(s_L, x_jp1) ;
		r_jp2 = stereographic_r(s_L, x_jp2) ;

		dr = stereographic_dr(s_L, x_j, dx) ;
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
void advance_tStep_massless_scalar(
	double s_L,
	int Nx, 
	double dt, double dx, 
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
		res = compute_iteration_GR_Crank_Nicolson_PQ(
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
	} while (res>ERR_TOLERANCE) ;

	Kreiss_Oliger_Filter(Nx, P_n) ;
	Kreiss_Oliger_Filter(Nx, Q_n) ;

	return ;
}	
/*==========================================================================*/
/* leftgoing Gaussian pulse */
/*==========================================================================*/
void initial_data_Gaussian(
	double s_L,
	int Nx, 	
	double dt, double dx,
	double bbox[2],
	char* direction,
	double amp, double width, double center,
	double* Al, double* Ze, 
	double*  P, double*  Q)
{
	double left_point = bbox[0] ;

	double x = 0 ;
	double r = 0 ;

	set_array_val(0, Nx, 1.0, Al) ;
	set_array_val(0, Nx, 0.0, Ze) ;

	int end_jC = (fabs(bbox[1]-s_L)<ERR_TOLERANCE) ? Nx-1 : Nx ;

	printf("bbox[1]\t%f\ts_L\t%f\twhole grid\t%d\tNx\t%d\tend_jC\t%d\n", bbox[1], s_L, (fabs(bbox[1]-s_L)<ERR_TOLERANCE), Nx, end_jC) ;

	for (int jC=0; jC<end_jC; jC++) {
		x = (jC * dx) + left_point ;
		r = stereographic_r(s_L, x) ;
		Q[jC] = amp * exp(-pow((r-center)/width,2)) * (
			(-(r-center)/pow(width,2)) * pow(r,2) 
		+	2*r
		) ;
		if (strcmp(direction,"ingoing")==0) {
			P[jC] = + Q[jC] ;
		} else if (strcmp(direction,"outgoing")==0) {
			P[jC] = - Q[jC] ;
		} else {
			fprintf(stderr,"ERROR(initial_data_Gaussian): direction!={{in,out}going}\n") ;
			exit(EXIT_FAILURE) ;
		}
	}
	if (end_jC == Nx-1) {
		P[Nx-1] = 0 ;
		Q[Nx-1] = 0 ;
	}
	return ;
}
