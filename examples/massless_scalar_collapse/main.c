#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <bbhutil.h>
#include <assert.h>

#include "amr_evolve.h"
#include "amr_grid_hierarchy.h"

#include "evolution_routines.h"

#define MAX_FILE_NAME 1024
#define OUTPUT_DIR "/home/jripley/one-dim-amr/examples/massless_scalar_collapse/output/"

/*===========================================================================*/
/* global variables for evolution-convenient for function calls */
/*===========================================================================*/
double* Al_n ;
double* Al_nm1 ;
double* Al_nm2 ;
double* Ze_n ;
double* Ze_nm1 ;
double* Ze_nm2 ;
double* P_n ;
double* P_nm1 ;
double* P_nm2 ;
double* Q_n ;
double* Q_nm1 ;
double* Q_nm2 ;

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

char output_name_Al_sdf[MAX_FILE_NAME+1] ;
char output_name_Ze_sdf[MAX_FILE_NAME+1] ;
char output_name_P_sdf[MAX_FILE_NAME+1] ;
char output_name_Q_sdf[MAX_FILE_NAME+1] ;

bool perim_interior[2] ;
bool excision_on = false ;
bool made_files  = false ;
/*===========================================================================*/
/* call after variables have been defined */
/*===========================================================================*/
void set_run_data(void)
{
	Nx = pow(2,8)+1 ;
	Nt = pow(2,14)+1 ;
	t_step_save = 2 ;

	perim_interior[0] = false ;
	perim_interior[1] = false ;

	cfl_num = 0.2 ;

	bbox[0] =   0 ;
	bbox[1] =  50 ;
	stereographic_L = bbox[1] ;

	dx = (bbox[1] - bbox[0]) / (Nx-1) ;
	dt = cfl_num * dx ;

	return ;
}
/*==========================================================================*/
/*number of time steps for evolution fields: 3 */
/*==========================================================================*/
amr_field* set_fields(void) 
{
	amr_field* fields 
	= amr_init_fields("Al", "ode", 3) ;
	
	amr_add_field(fields, "Ze", "ode", 3) ;

	amr_add_field(fields, "P", "hyperbolic", 3) ;
	amr_add_field(fields, "Q", "hyperbolic", 3) ;

	return fields ;
}
void find_field_indices(amr_field* fields) 
{
	Al_n_index   = amr_find_field_index(fields, "Al") ;
	Al_nm1_index = Al_n_index + 1 ;
	Al_nm2_index = Al_n_index + 2 ;

	Ze_n_index   = amr_find_field_index(fields, "Ze") ;
	Ze_nm1_index = Ze_n_index + 1 ;
	Ze_nm2_index = Ze_n_index + 2 ;

	P_n_index   = amr_find_field_index(fields, "P") ;
	P_nm1_index = P_n_index + 1 ;
	P_nm2_index = P_n_index + 2 ;

	Q_n_index   = amr_find_field_index(fields, "Q") ;
	Q_nm1_index = Q_n_index + 1 ;
	Q_nm2_index = Q_n_index + 2 ;

	assert(Al_n_index >= 0) ;
	assert(Ze_n_index >= 0) ;
	assert(P_n_index >= 0) ;
	assert(Q_n_index >= 0) ;

	return ;
}
/*===========================================================================*/
/* call to set global variables for field evolution */
/*===========================================================================*/
void set_globals(amr_grid* grid)
{	

	P_n   = grid->grid_funcs[P_n_index  ] ;
	P_nm1 = grid->grid_funcs[P_nm1_index] ;
	P_nm2 = grid->grid_funcs[P_nm2_index] ;
	Q_n   = grid->grid_funcs[Q_n_index  ] ;
	Q_nm1 = grid->grid_funcs[Q_nm1_index] ;
	Q_nm2 = grid->grid_funcs[Q_nm2_index] ;

	Al_n   = grid->grid_funcs[Al_n_index  ] ;
	Al_nm1 = grid->grid_funcs[Al_nm1_index] ;
	Al_nm2 = grid->grid_funcs[Al_nm2_index] ;
	Ze_n   = grid->grid_funcs[Ze_n_index  ] ;
	Ze_nm1 = grid->grid_funcs[Ze_nm1_index] ;
	Ze_nm2 = grid->grid_funcs[Ze_nm2_index] ;

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
void initial_data(amr_grid* grid)
{
	set_globals(grid) ;

	initial_data_Gaussian(
		stereographic_L,
		Nx, 	
		dt, dx,
		excised_jC,
		bbox,
		perim_interior,
		Al_n, Ze_n, 
		P_n,   Q_n)
	;
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void wave_evolve(amr_grid* grid)
{
	set_globals(grid) ;

	advance_tStep_massless_scalar(
		stereographic_L,
		Nx, dt, dx, 
		excised_jC,
		bbox, perim_interior,
		Al_n, Al_nm1, Ze_n, Ze_nm1,
		 P_n,  P_nm1,  Q_n,  Q_nm1
	) ;	
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void save_to_file(amr_grid* grid)
{
	set_globals(grid) ;

	if (made_files == false) {
		snprintf(output_name_Al_sdf, MAX_FILE_NAME, "%sAl.sdf", OUTPUT_DIR) ;
		snprintf(output_name_Ze_sdf, MAX_FILE_NAME, "%sZe.sdf", OUTPUT_DIR) ;
		snprintf(output_name_P_sdf, MAX_FILE_NAME, "%sP.sdf", OUTPUT_DIR) ;
		snprintf(output_name_Q_sdf, MAX_FILE_NAME, "%sQ.sdf", OUTPUT_DIR) ;
	}
	made_files = true ;
	gft_out_bbox(output_name_Al_sdf, time, &Nx, 1, bbox, Al_n) ;
	gft_out_bbox(output_name_Ze_sdf, time, &Nx, 1, bbox, Ze_n) ;
	gft_out_bbox(output_name_P_sdf, time, &Nx, 1, bbox, P_n) ;
	gft_out_bbox(output_name_Q_sdf, time, &Nx, 1, bbox, Q_n) ;
	
	return ;
}
/*===========================================================================*/
/* call amr evolution */
/*===========================================================================*/

int main(int argc, char* argv[])
{
	set_run_data() ;

	amr_field* fields = set_fields() ; 
	find_field_indices(fields) ; 

	amr_grid_hierarchy* gh 
	= amr_init_grid_hierarchy(
		fields,
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
