#ifndef _FILE_IO_H_
#define _FILE_IO_H_

#include <stdbool.h>

/*==========================================================================*/
/*	these max lengths are essentially arbitrary at the moment 
*/
#define MAX_NAME_LEN 1024
#define MAX_LINE_LEN 4096

/*==========================================================================*/
void get_run_data(
	char** theory,
	char** output_dir,
	int* Nx, int* Nt, int* t_step_save,
	double* cfl_num, 
	double* bbox_0, double* bbox_1,
	double* stereographic_L,
	double* dt, double* dx,
	double* err_tolerance) 
; 
/*==========================================================================*/
void get_initial_data(
	char** initial_data,
	char** direction,
	double* amp, double* width, double* center,
	double* initial_black_hole_mass) 
;
/*==========================================================================*/
#endif
