#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include "evolution_routines_GR.h"
#include "stencils.h"

#define ERR_TOLERANCE ((double)1e-10)
#define MACHINE_EPSILON ((double)1e-14)

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
/* Gaussian pulse */
/*==========================================================================*/
void initial_data_Gaussian(
	double s_L,
	int Nx, 	
	double dx,
	double bbox[2],
	char *direction,
	double amp, double width, double center,
	double *Al, double *Ze, 
	double  *P, double  *Q,
	double *phi_n, double *phi_nm1)
{
	double x_lower = bbox[0] ;

	double x = 0 ;
	double r = 0 ;

	set_array_val(0, Nx, 1.0, Al) ;
	set_array_val(0, Nx, 0.0, Ze) ;

	int end_jC = (fabs(bbox[1]-s_L)<ERR_TOLERANCE) ? Nx-1 : Nx ;

	printf("bbox[1]\t%f\ts_L\t%f\twhole grid\t%d\tNx\t%d\tend_jC\t%d\n", bbox[1], s_L, (fabs(bbox[1]-s_L)<ERR_TOLERANCE), Nx, end_jC) ;

	for (int jC=0; jC<end_jC; jC++) {
		x = (jC * dx) + x_lower ;
		r = stereographic_r(s_L, x) ;

		double gaussian = pow(width,-4) * amp * exp(-pow((r-center)/width,2)) ;

		phi_n[jC] = pow(r,4) * gaussian ; 
		phi_nm1[jC] = phi_n[jC] ; 

		Q[jC] = gaussian * (
		-	2 * pow(r,4) * ((r-center)/pow(width,2))  
		+	4 * pow(r,3)
		) ;
		if (strcmp(direction,"ingoing")==0) {
			P[jC] = + Q[jC] ;
		} else if (strcmp(direction,"outgoing")==0) {
			P[jC] = - Q[jC] ;
		} else if (strcmp(direction,"stationary")==0) {
			P[jC]=0;
		} else {
			fprintf(stderr,"ERROR(initial_data_Gaussian): direction!={{in,out}going, stationary}\n") ;
			exit(EXIT_FAILURE) ;
		}
	}
	if (end_jC == Nx-1) {
		phi_n[Nx-1] = 0 ;
		phi_nm1[Nx-1] = phi_n[Nx-1] ;
		Q[Nx-1] = 0 ;
		P[Nx-1] = 0 ;
	}
	return ;
}
/*==========================================================================*/
/* initial black hole (in Painleve-Gullstrand coordinates) */
/*==========================================================================*/
void initial_data_black_hole(
	double s_L,
	int Nx, 	
	double dx,
	double bbox[2],
	double mass,
	int exc_jC,
	double* Al, double* Ze)
{
	double x_lower = bbox[0] ;

	double x = 0 ;
	double r = 0 ;

	set_array_val(0, Nx, 1.0, Al) ;
	set_array_val(0, Nx, 0.0, Ze) ;

	int end_jC = (fabs(bbox[1]-s_L)<ERR_TOLERANCE) ? Nx-1 : Nx ;

	printf("exc_jC= %d\n", exc_jC) ;

	for (int jC=exc_jC; jC<end_jC; jC++) {
		x = (jC * dx) + x_lower ;
		r = stereographic_r(s_L, x) ;

		Ze[jC] = sqrt(2*mass/r) ;
	}
	if (end_jC == Nx-1) {
		Ze[Nx-1] = 0 ;
	}
	return ;
}
