#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <bbhutil.h>

#include "amr_evolve.h"
#include "amr_grid_hierarchy.h"

#include "evolution_routines.h"

#define MAX_FILE_NAME 1024
#define OUTPUT_DIR "/home/jripley/one-dim-amr/examples/massless_scalar_collapse/output/"

/*===========================================================================*/
/* global variables for evolution-convenient for function calls */
/*===========================================================================*/
double* Al_n ;
double* Ze_n ;
double* P_n ;
double* Q_n ;

double* Al_nm1 ;
double* Ze_nm1 ;
double* P_nm1 ;
double* Q_nm1 ;

double cfl_num ;
double bbox[2] ;
double dt, dx ;
double time ;
double stereographic_L ; /* stereographic length: for compactification */

int Nx, Nt, t_step_save ;
int excised_jC ;

int num_fields ;
int 
	Al_n_index, Al_nm1_index, Al_nm2_index,
	Ze_n_index, Ze_nm1_index, Ze_nm2_index,

	P_n_index, P_nm1_index, P_nm2_index,
	Q_n_index, Q_nm1_index, Q_nm2_index
;
int perim_coords[2] ;

char output_name_P_sdf[MAX_FILE_NAME+1] ;
char output_name_Q_sdf[MAX_FILE_NAME+1] ;

bool perim_interior[2] ;
bool excision_on = false ;
bool made_files  = false ;
/*===========================================================================*/
/* call after variables have been defined */
/*===========================================================================*/
void set_initial_run_data(void)
{
	Nx = pow(2,8)+1 ;
	Nt = pow(2,13)+1 ;
	t_step_save = 1 ;

	perim_interior[0] = false ;
	perim_interior[1] = false ;

	cfl_num = 0.05 ;

	bbox[0] =   0 ;
	bbox[1] =  50 ;
	stereographic_L = bbox[1] ;

	dx = (bbox[1] - bbox[0]) / (Nx-1) ;
	dt = cfl_num * dx ;

	return ;
}
void set_field_indices(void) 
{

	P_n_index   = 0 ;
	P_nm1_index = 1 ;
	Q_n_index   = 2 ;
	Q_nm1_index = 3 ;

	Al_n_index   = 4 ;
	Al_nm1_index = 5 ;
	Ze_n_index   = 6 ;
	Ze_nm1_index = 7 ;

	num_grid_funcs = 8 ;

	return ;
}
/*===========================================================================*/
/* call to set global variables for field evolution */
/*===========================================================================*/
void set_globals(struct amr_grid* grid)
{	

	P_n   = grid->grid_funcs[P_n_index  ] ;
	P_nm1 = grid->grid_funcs[P_nm1_index] ;
	Q_n   = grid->grid_funcs[Q_n_index  ] ;
	Q_nm1 = grid->grid_funcs[Q_nm1_index] ;

	Al_n   = grid->grid_funcs[Al_n_index  ] ;
	Al_nm1 = grid->grid_funcs[Al_nm1_index] ;
	Ze_n   = grid->grid_funcs[Ze_n_index  ] ;
	Ze_nm1 = grid->grid_funcs[Ze_nm1_index] ;

	Nx = grid->Nx ;

	dt = grid->dt ;
	dx = grid->dx ;

	time = grid->time ;

	bbox[0] = grid->bbox[0] ;
	bbox[1] = grid->bbox[1] ;

	perim_interior[0] = grid->perim_interior[0] ;
	perim_interior[1] = grid->perim_interior[1] ;
	
	perim_coords[0] = grid->perim_coords[0] ;
	perim_coords[1] = grid->perim_coords[1] ;

	excised_jC = grid->excised_jC ;

	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void initial_data(struct amr_grid* grid)
{
	set_globals(grid) ;

	initial_data_Gaussian(
		Nx, 	dx,
		bbox,
		Al_n, Al_nm1, Ze_n, Ze_nm1,
		 P_n,  P_nm1,  Q_n,  Q_nm1)
	;	
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void wave_evolve(struct amr_grid* grid)
{
	set_globals(grid) ;

	advance_tStep_massless_scalar(
		stereographic_L,
		Nx, dt, dx, bbox, perim_interior,
		Al_n, Al_nm1, Ze_n, Ze_nm1,
		 P_n,  P_nm1,  Q_n,  Q_nm1
	) ;	
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void save_to_file(struct amr_grid* grid)
{
	set_globals(grid) ;

	if (made_files == false) {
		snprintf(output_name_P_sdf, MAX_FILE_NAME, "%sP.sdf", OUTPUT_DIR) ;
		snprintf(output_name_Q_sdf, MAX_FILE_NAME, "%sQ.sdf", OUTPUT_DIR) ;
	}
	made_files = true ;

	gft_out_bbox(output_name_P_sdf, time, &Nx, 1, bbox, P_n) ;
	gft_out_bbox(output_name_Q_sdf, time, &Nx, 1, bbox, Q_n) ;
	
	return ;
}
/*===========================================================================*/
/* call amr evolution */
/*===========================================================================*/

int main(int argc, char* argv[])
{
	set_initial_run_data() ;
	set_field_indices() ; 

	struct amr_grid_hierarchy* gh 
	= amr_init_grid_hierarchy(
        	num_grid_funcs,
        	Nt, Nx, t_step_save,
		cfl_num,
		bbox,
		excision_on)
	;
	amr_main(
		gh, 
		initial_data,
		wave_evolve,
		save_to_file)
	;
	amr_destroy_grid_hierarchy(gh) ;

	return EXIT_SUCCESS ;
}
