#ifndef _EVOLUTION_ROUTINES_H_
#define _EVOLUTION_ROUTINES_H_

#include <stdbool.h>

void initial_data_Gaussian(
	double s_L,
	int Nx, 	
	double dt, double dx,
	int exc_jC,
	double bbox[2],
	bool perim_interior[2],
	double* Al, double* Ze, 
	double*  P, double*  Q)
;
void initial_data_Gaussian_GR(
	double s_L,
	int Nx, 	
	double dt, double dx,
	int exc_jC,
	double bbox[2],
	bool perim_interior[2],
	double* Al, double* Ze, 
	double*  P, double*  Q)
;
/* we set Dirichlet boundary conditions for now */
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
;
void solve_Al_Ze(
	double s_L,
	int Nx,
	double dt, 	double dx,
	bool excision_on,
	int exc_jC,
	double bbox[2],
	bool perim_interior[2],
	double* Al_n, 	double* Al_nm1, double* Ze_n, double* Ze_nm1,
	double*  P_n, 	double*  P_nm1, double*  Q_n, double*  Q_nm1)
;
#endif
