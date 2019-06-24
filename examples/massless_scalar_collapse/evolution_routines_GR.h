#ifndef _EVOLUTION_ROUTINES_H_
#define _EVOLUTION_ROUTINES_H_

#include <stdbool.h>

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
	double *Al_n, double *Al_nm1, double *Ze_n, double *Ze_nm1,
	double  *P_n, double  *P_nm1, double  *Q_n, double  *Q_nm1)
;
/*===========================================================================*/
void solve_Al_Ze_GR(
	double s_L,
	int Nx,
	double dt, 	double dx,
	double err_tolerance,
	bool excision_on,
	int start_jC,
	double bbox[2],
	double *Al_n, 	double *Al_nm1, double *Ze_n, double *Ze_nm1,
	double  *P_n, 	double  *P_nm1, double  *Q_n, double  *Q_nm1)
;
#endif
