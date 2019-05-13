#ifndef _EVOLUTION_ROUTINES_H_
#define _EVOLUTION_ROUTINES_H_

#include <stdbool.h>

void copy_to_2nd_array(int Nx, double* array_1, double* array_2) 
;
void initial_Data(
	int Nx, 	double dx,
	double left_point,
	double* P_n, 	double* Q_n)
;
/* we set Dirichlet boundary conditions for now */
void advance_tStep_wave(
	int Nx,
	double dt, 	double dx,
	bool perim_interior[2],
	double* P_n, 	double* P_nm1,
	double* Q_n,	double* Q_nm1)
;
void save_to_txt_file(
	int Nx,
	FILE* output,
	double* field)
;
#endif
