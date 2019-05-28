#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <bbhutil.h>
#include <assert.h>
/* for checking if output directory exists */
#include <dirent.h>
#include <errno.h>

#include "amr_evolve.h"
#include "amr_grid_hierarchy.h"
#include "diagnostics_general.h"
#include "diagnostics_GR.h"
#include "evolution_routines_GR.h"
#include "file_io.h"


/*===========================================================================*/
/* global variables: set here (will eventually set by reading from
 * initial data file */
/*===========================================================================*/
/*===========================================================================*/
/* global variables for evolution-convenient for function calls */
/*===========================================================================*/
/*---------------------------------------------------------------------------*/
/* ODEs: extr: for extrapolation (saved when levels merge) */
/*---------------------------------------------------------------------------*/
double *Al_n, *Al_nm1, *Al_nm2, *Al_extr_m1, *Al_extr_m2 ;
double *Ze_n, *Ze_nm1, *Ze_nm2, *Ze_extr_m1, *Ze_extr_m2 ;
/*---------------------------------------------------------------------------*/
/* Hyperbolic variables */
/*---------------------------------------------------------------------------*/
double  *P_n,  *P_nm1,  *P_nm2 ;
double  *Q_n,  *Q_nm1,  *Q_nm2 ;
/*---------------------------------------------------------------------------*/
/* diagnostics */
/*---------------------------------------------------------------------------*/
double *mass_aspect ;
double *ingoing_null_characteristic, *outgoing_null_characteristic ;
double *Ricci_scalar, *Gauss_Bonnet_scalar ;
double *eom_TR, *eom_ThTh, *eom_scalar ;
/*---------------------------------------------------------------------------*/
/* run data */
/*---------------------------------------------------------------------------*/
char *theory, *output_dir ;

double cfl_num ;
double bbox[2] ;
double dt, dx ;
double time ;
double stereographic_L ; /* stereographic length: for compactification */

int Nx, Nt, t_step_save ;
int excised_jC ;
/*---------------------------------------------------------------------------*/
/* initial data */
/*---------------------------------------------------------------------------*/
char* id_type ;

double amp, width, center ;
/*---------------------------------------------------------------------------*/
int num_fields ;

int 	Al_n_index, Al_nm1_index, Al_nm2_index, Al_extr_m1_index, Al_extr_m2_index,
	Ze_n_index, Ze_nm1_index, Ze_nm2_index, Ze_extr_m1_index, Ze_extr_m2_index,

	P_n_index, P_nm1_index, P_nm2_index,
	Q_n_index, Q_nm1_index, Q_nm2_index,

	mass_aspect_index,
	ingoing_null_characteristic_index, outgoing_null_characteristic_index,
	Ricci_scalar_index, Gauss_Bonnet_scalar_index,
	eom_TR_index, eom_ThTh_index, eom_scalar_index
;
int perim_coords[2] ;
/*---------------------------------------------------------------------------*/
char output_name_Al[MAX_NAME_LEN+1] ;
char output_name_Ze[MAX_NAME_LEN+1] ;
char output_name_P[MAX_NAME_LEN+1] ;
char output_name_Q[MAX_NAME_LEN+1] ;

char output_name_ingoing_null_characteristic[MAX_NAME_LEN+1] ;
char output_name_outgoing_null_characteristic[MAX_NAME_LEN+1] ;

char output_name_Ricci_scalar[MAX_NAME_LEN+1] ;
char output_name_Gauss_Bonner_scalar[MAX_NAME_LEN+1] ;

char output_name_mass_aspect[MAX_NAME_LEN+1] ;

char output_name_eom_TR[MAX_NAME_LEN+1] ;
char output_name_eom_ThTh[MAX_NAME_LEN+1] ;
char output_name_eom_scalar[MAX_NAME_LEN+1] ;
/*---------------------------------------------------------------------------*/
bool perim_interior[2] ;
bool excision_on = true ;
bool made_output_files  = false ;
/*==========================================================================*/
/*number of time steps for evolution fields: 3.
 * Two extrapolation levels for linear extrapolation of ode variables */
/*==========================================================================*/
amr_field* set_fields(void) 
{
	int time_levels = 3 ;
	int ode_extrap_levels = 2 ;

	amr_field* fields 
	= amr_init_fields("Al", "ode", time_levels, ode_extrap_levels) ;
	
	amr_add_field(fields, "Ze", "ode", time_levels, ode_extrap_levels) ;

	amr_add_field(fields, "P", "hyperbolic", time_levels, 0) ;
	amr_add_field(fields, "Q", "hyperbolic", time_levels, 0) ;

	amr_add_field(fields, "mass_aspect",  "diagnostic", time_levels, 0) ;

	amr_add_field(fields, "ingoing_null_characteristic",  "diagnostic", 1, 0) ;
	amr_add_field(fields, "outgoing_null_characteristic", "diagnostic", 1, 0) ;

	amr_add_field(fields, "Ricci_scalar",        "diagnostic", 1, 0) ;
	amr_add_field(fields, "Gauss_Bonnet_scalar", "diagnostic", 1, 0) ;

	amr_add_field(fields, "eom_TR",     "diagnostic", 1, 0) ;
	amr_add_field(fields, "eom_ThTh",   "diagnostic", 1, 0) ;
	amr_add_field(fields, "eom_scalar", "diagnostic", 1, 0) ;

	return fields ;
}
void find_field_indices(amr_field* fields) 
{
	Al_n_index   = amr_return_field_index(fields, "Al") ;
	Al_nm1_index = Al_n_index + 1 ;
	Al_nm2_index = Al_n_index + 2 ;
	Al_extr_m1_index = Al_n_index + 3 ;
	Al_extr_m2_index = Al_n_index + 4 ;

	Ze_n_index   = amr_return_field_index(fields, "Ze") ;
	Ze_nm1_index = Ze_n_index + 1 ;
	Ze_nm2_index = Ze_n_index + 2 ;
	Ze_extr_m1_index = Ze_n_index + 3 ;
	Ze_extr_m2_index = Ze_n_index + 4 ;

	P_n_index   = amr_return_field_index(fields, "P") ;
	P_nm1_index = P_n_index + 1 ;
	P_nm2_index = P_n_index + 2 ;	

	Q_n_index   = amr_return_field_index(fields, "Q") ;
	Q_nm1_index = Q_n_index + 1 ;
	Q_nm2_index = Q_n_index + 2 ;

	mass_aspect_index  = amr_return_field_index(fields, "mass_aspect")  ;

	ingoing_null_characteristic_index  = amr_return_field_index(fields, "ingoing_null_characteristic")  ;
	outgoing_null_characteristic_index = amr_return_field_index(fields, "outgoing_null_characteristic") ;

	Ricci_scalar_index = amr_return_field_index(fields, "Ricci_scalar")  ;
	Gauss_Bonnet_scalar_index = amr_return_field_index(fields, "Gauss_Bonnet_scalar") ;

	eom_TR_index = amr_return_field_index(fields, "eom_TR") ;
	eom_ThTh_index = amr_return_field_index(fields, "eom_ThTh") ;
	eom_scalar_index = amr_return_field_index(fields, "eom_scalar") ;

	assert(Ze_n_index == Al_extr_m2_index+1) ;
	assert(P_n_index == Ze_extr_m2_index+1) ;
	assert(Q_n_index == P_nm2_index+1) ;
	assert(mass_aspect_index == Q_nm2_index+1) ;

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

	mass_aspect = grid->grid_funcs[mass_aspect_index] ;

	ingoing_null_characteristic  = grid->grid_funcs[ingoing_null_characteristic_index] ;
	outgoing_null_characteristic = grid->grid_funcs[outgoing_null_characteristic_index] ;

	Ricci_scalar = grid->grid_funcs[Ricci_scalar_index] ;
	Gauss_Bonnet_scalar = grid->grid_funcs[Gauss_Bonnet_scalar_index] ;

	eom_TR = grid->grid_funcs[eom_TR_index] ;
	eom_ThTh = grid->grid_funcs[eom_ThTh_index] ;
	eom_scalar = grid->grid_funcs[eom_scalar_index] ;

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

	excision_on = grid->excision_on ;
	excised_jC  = grid->excised_jC ;

	return ;
}
/*===========================================================================*/
/* free initial data */
/*===========================================================================*/
void set_free_initial_data(amr_grid* grid)
{
	set_globals(grid) ;

	initial_data_Gaussian(
		stereographic_L,
		Nx, 	
		dt, dx,
		bbox,
		amp, width, center,
		Al_n, Ze_n, 
		P_n, Q_n)
	;
	return ;
}
/*===========================================================================*/
/* rescale Al on all levels so value at spatial infinity is unity.
 * level is the finest grid the spans the entire computational domain */
/*===========================================================================*/
void rescale_Al(amr_grid* grid)
{
	double rescale_param = 1 ;
	int level = grid->level ;
	if (level==1) {
		rescale_param = grid->grid_funcs[Al_n_index][Nx-1] ;
		for (amr_grid* iter=(grid->parent); iter!=NULL; iter=iter->child) {
			for (int iC=0; iC<(iter->Nx); iC++) {
				iter->grid_funcs[Al_n_index][iC] /= rescale_param ;
			}
		}
	}
}
/*===========================================================================*/
/*===========================================================================*/
void solve_ode(amr_grid* grid)
{
	int child_perim_coords[2] = {-1,-1} ;
	if (grid->child!=NULL) {
		child_perim_coords[0] = grid->child->perim_coords[0] ;
		child_perim_coords[1] = grid->child->perim_coords[1] ;
	}
	if (strcmp(theory, "massless_scalar_GR") == 0) {
		set_globals(grid) ;
		solve_Al_Ze(
			stereographic_L,
			Nx,
			dt, 	dx,
			excision_on,
			excised_jC,
			child_perim_coords,
			bbox,
			Al_n, 	Al_nm1, Ze_n, Ze_nm1,
			 P_n, 	 P_nm1,  Q_n,  Q_nm1)
		;
		rescale_Al(grid) ;
	}
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void evolve_hyperbolic_pde(amr_grid* grid)
{
	set_globals(grid) ;
	advance_tStep_massless_scalar(
		stereographic_L,
		Nx, dt, dx, 
		excision_on,
		excised_jC,
		bbox, perim_interior,
		Al_n, Al_nm1, Ze_n, Ze_nm1,
		 P_n,  P_nm1,  Q_n,  Q_nm1
	) ;	
	return ;
}
/*===========================================================================*/
/* independent residuals, curvature scalars, etc */
/*===========================================================================*/
void compute_diagnostics(amr_grid* grid)
{
	set_globals(grid) ;
/*--------------------------------------------------------------------------*/
	compute_diagnostics_massless_scalar_GR(
		Nx, excised_jC, 
		stereographic_L,
		dt, dx,
		Al_n, Al_nm1, Al_nm2,
		Ze_n, Ze_nm1, Ze_nm2,
		P_n, P_nm1, P_nm2,
		Q_n, Q_nm1, Q_nm2,
		eom_TR,
		eom_ThTh,
		eom_scalar)
	;
/*--------------------------------------------------------------------------*/
/* ''_general: includes setting excised_jC */ 
/*--------------------------------------------------------------------------*/
	compute_diagnostics_general(
		Nx, 
		stereographic_L,
		dt, dx,
		Al_n, Al_nm1, Al_nm2,
		Ze_n, Ze_nm1, Ze_nm2,
		&excised_jC,
		mass_aspect,
		ingoing_null_characteristic,
		outgoing_null_characteristic,
		Ricci_scalar,
		Gauss_Bonnet_scalar)
	;
	if (grid->parent==NULL) { 
		grid->excised_jC = excised_jC ;
		printf("t\t%f\tMA\t%f\texc_jC\t%d\n", (grid->tC)*dt, mass_aspect[Nx-1], excised_jC) ;
	}
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void save_to_file(amr_grid* grid)
{
	set_globals(grid) ;

	if (made_output_files == false) {
	/* 	confirm that output directory exists
	*/
		DIR* dir = opendir(output_dir) ;
		if (dir!=NULL) {
			closedir(dir) ;
		} else if (ENOENT==errno) {
			printf("ERROR(main.c): unable to open %s\n", output_dir) ;
			printf("directory does not exist?\n") ;
			exit(EXIT_FAILURE) ;
		} else {
			printf("ERROR(main.c): unable to open %s\n", output_dir) ;
			exit(EXIT_FAILURE) ;
		}
	/* 	print output file names
	*/
		snprintf(output_name_Al, MAX_NAME_LEN, "%sAl.sdf", output_dir) ;
		snprintf(output_name_Ze, MAX_NAME_LEN, "%sZe.sdf", output_dir) ;
		snprintf(output_name_P,  MAX_NAME_LEN, "%sP.sdf",  output_dir) ;
		snprintf(output_name_Q,  MAX_NAME_LEN, "%sQ.sdf",  output_dir) ;

		snprintf(output_name_mass_aspect, MAX_NAME_LEN, "%smass_aspect.sdf", output_dir) ;

		snprintf(output_name_ingoing_null_characteristic,  MAX_NAME_LEN, "%singoing_null_characteristic.sdf", output_dir) ;
		snprintf(output_name_outgoing_null_characteristic, MAX_NAME_LEN, "%soutgoing_null_characteristic.sdf",  output_dir) ;
		snprintf(output_name_Ricci_scalar,        MAX_NAME_LEN, "%sRicci_scalar.sdf",        output_dir) ;
		snprintf(output_name_Gauss_Bonner_scalar, MAX_NAME_LEN, "%sGauss_Bonner_scalar.sdf", output_dir) ;

		snprintf(output_name_eom_TR, MAX_NAME_LEN, "%seom_TR.sdf", output_dir) ;
		snprintf(output_name_eom_ThTh, MAX_NAME_LEN, "%seom_ThTh.sdf", output_dir) ;
		snprintf(output_name_eom_scalar, MAX_NAME_LEN, "%seom_scalar.sdf", output_dir) ;
	}
	made_output_files = true ;
/* 	save to file-see man pages for bbhutil utilities on Choptuik's website 
*/
	gft_out_bbox(output_name_Al, time, &Nx, 1, bbox, Al_n) ;
	gft_out_bbox(output_name_Ze, time, &Nx, 1, bbox, Ze_n) ;
	gft_out_bbox(output_name_P,  time, &Nx, 1, bbox,  P_n) ;
	gft_out_bbox(output_name_Q,  time, &Nx, 1, bbox,  Q_n) ;

	gft_out_bbox(output_name_mass_aspect,  time, &Nx, 1, bbox, mass_aspect) ;

	gft_out_bbox(output_name_ingoing_null_characteristic,  time, &Nx, 1, bbox, ingoing_null_characteristic) ;
	gft_out_bbox(output_name_outgoing_null_characteristic, time, &Nx, 1, bbox, outgoing_null_characteristic) ;
	gft_out_bbox(output_name_Ricci_scalar,        time, &Nx, 1, bbox,  Ricci_scalar) ;
	gft_out_bbox(output_name_Gauss_Bonner_scalar, time, &Nx, 1, bbox,  Gauss_Bonnet_scalar) ;

	gft_out_bbox(output_name_eom_TR, time, &Nx, 1, bbox, eom_TR) ;
	gft_out_bbox(output_name_eom_ThTh, time, &Nx, 1, bbox, eom_ThTh) ;
	gft_out_bbox(output_name_eom_scalar, time, &Nx, 1, bbox, eom_scalar) ;
	
	return ;
}
/*===========================================================================*/
/* call amr evolution */
/*===========================================================================*/

int main(int argc, char* argv[])
{
	get_run_data(
		&theory,
		&output_dir,
		&Nx, &Nt, &t_step_save,
		&cfl_num, 
		&bbox[0], &bbox[1],
		&stereographic_L,
		&dt, &dx) 
	; 
	get_initial_data(&id_type, &amp, &width, &center) ; 

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
		set_free_initial_data,
		evolve_hyperbolic_pde,
		solve_ode,
		compute_diagnostics,
		save_to_file)
	;
	amr_destroy_grid_hierarchy(gh) ;
/*--------------------------------------------------------------------------*/
/* theory and id_type are malloc'ed in get_run_data and get_initial_data,
 * respectively. */ 
/*--------------------------------------------------------------------------*/
	free(theory) ;
	free(id_type) ;
/*--------------------------------------------------------------------------*/

	return EXIT_SUCCESS ;
}
