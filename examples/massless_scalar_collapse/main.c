#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <bbhutil.h>
#include <assert.h>
#include <time.h>
/* for checking if output directory exists */
#include <dirent.h>
#include <errno.h>

#include "amr_evolve.h"
#include "amr_grid_hierarchy.h"
#include "diagnostics_general.h"
#include "diagnostics_GR.h"
#include "diagnostics_EdGB.h"
#include "free_initial_data.h"
#include "evolution_routines_GR.h"
#include "evolution_routines_EdGB.h"
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
double *phi_n, *phi_nm1, *phi_nm2 ;
double *mass_aspect ;
double *ingoing_null_characteristic, *outgoing_null_characteristic ;
double *ingoing_scalar_characteristic, *outgoing_scalar_characteristic ;
double *Ricci_scalar, *Gauss_Bonnet_scalar ;
double *eom_TR, *eom_ThTh, *eom_scalar ;
/*---------------------------------------------------------------------------*/
/* run data */
/*---------------------------------------------------------------------------*/
char *theory, *output_dir ;
char *solver_Ze ;

double cfl_num ;
double bbox[2] ;
double dt, dx ;
double grid_time ;
double stereographic_L ; /* stereographic length: for compactification */
double coupling_gbs ; /* for EdGB theory. gbs: Gauss-Bonnet scalar */

double err_tolerance ;

int Nx, Nt, t_step_save ;
int excised_jC ;
/*---------------------------------------------------------------------------*/
/* initial data */
/*---------------------------------------------------------------------------*/
char *initial_data, *direction ;

double amp, width, center ;
double initial_black_hole_mass ;
/*---------------------------------------------------------------------------*/
int num_fields ;

int 	Al_n_index, Al_nm1_index, Al_nm2_index, Al_extr_m1_index, Al_extr_m2_index,
	Ze_n_index, Ze_nm1_index, Ze_nm2_index, Ze_extr_m1_index, Ze_extr_m2_index,

	P_n_index, P_nm1_index, P_nm2_index,
	Q_n_index, Q_nm1_index, Q_nm2_index,

	phi_n_index, phi_nm1_index, phi_nm2_index,

	mass_aspect_index,
	ingoing_null_characteristic_index, outgoing_null_characteristic_index,
	ingoing_scalar_characteristic_index, outgoing_scalar_characteristic_index,
	Ricci_scalar_index, Gauss_Bonnet_scalar_index,
	eom_TR_index, eom_ThTh_index, eom_scalar_index
;
int perim_coords[2] ;
/*---------------------------------------------------------------------------*/
char output_name_Al[MAX_NAME_LEN+1] ;
char output_name_Ze[MAX_NAME_LEN+1] ;
char output_name_P[MAX_NAME_LEN+1] ;
char output_name_Q[MAX_NAME_LEN+1] ;

char output_name_phi[MAX_NAME_LEN+1] ;

char output_name_ingoing_null_characteristic[MAX_NAME_LEN+1] ;
char output_name_outgoing_null_characteristic[MAX_NAME_LEN+1] ;

char output_name_ingoing_scalar_characteristic[MAX_NAME_LEN+1] ;
char output_name_outgoing_scalar_characteristic[MAX_NAME_LEN+1] ;

char output_name_Ricci_scalar[MAX_NAME_LEN+1] ;
char output_name_Gauss_Bonnet_scalar[MAX_NAME_LEN+1] ;

char output_name_mass_aspect[MAX_NAME_LEN+1] ;

char output_name_eom_TR[MAX_NAME_LEN+1] ;
char output_name_eom_ThTh[MAX_NAME_LEN+1] ;
char output_name_eom_scalar[MAX_NAME_LEN+1] ;

char output_name_test_1[MAX_NAME_LEN+1] ;
char output_name_test_2[MAX_NAME_LEN+1] ;
/*---------------------------------------------------------------------------*/
bool perim_interior[2] ;
bool excision_on = true ;
bool made_output_files  = false ;
/*==========================================================================*/
/*number of time steps for evolution fields: 3.
 * Two extrapolation levels for linear extrapolation of ode variables */
/*==========================================================================*/
amr_field *set_fields(void) 
{
	int time_levels = 3 ;
	int ode_extrap_levels = 2 ;

	amr_field* fields 
	= amr_init_fields("Al", "ode", time_levels, ode_extrap_levels) ;
	
	amr_add_field(fields, "Ze", "ode", time_levels, ode_extrap_levels) ;

	amr_add_field(fields, "P", "hyperbolic", time_levels, 0) ;
	amr_add_field(fields, "Q", "hyperbolic", time_levels, 0) ;

	amr_add_field(fields, "phi", "diagnostic", time_levels, 0) ;

	amr_add_field(fields, "mass_aspect",  "diagnostic", 1, 0) ;

	amr_add_field(fields, "ingoing_null_characteristic",  "diagnostic", 1, 0) ;
	amr_add_field(fields, "outgoing_null_characteristic", "diagnostic", 1, 0) ;

	amr_add_field(fields, "ingoing_scalar_characteristic",  "diagnostic", 1, 0) ;
	amr_add_field(fields, "outgoing_scalar_characteristic", "diagnostic", 1, 0) ;

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

	phi_n_index   = amr_return_field_index(fields, "phi") ;
	phi_nm1_index = phi_n_index + 1 ;
	phi_nm2_index = phi_n_index + 2 ;

	mass_aspect_index  = amr_return_field_index(fields, "mass_aspect")  ;

	ingoing_null_characteristic_index  = amr_return_field_index(fields, "ingoing_null_characteristic")  ;
	outgoing_null_characteristic_index = amr_return_field_index(fields, "outgoing_null_characteristic") ;

	ingoing_scalar_characteristic_index  = amr_return_field_index(fields, "ingoing_scalar_characteristic")  ;
	outgoing_scalar_characteristic_index = amr_return_field_index(fields, "outgoing_scalar_characteristic") ;

	Ricci_scalar_index = amr_return_field_index(fields, "Ricci_scalar")  ;
	Gauss_Bonnet_scalar_index = amr_return_field_index(fields, "Gauss_Bonnet_scalar") ;

	eom_TR_index = amr_return_field_index(fields, "eom_TR") ;
	eom_ThTh_index = amr_return_field_index(fields, "eom_ThTh") ;
	eom_scalar_index = amr_return_field_index(fields, "eom_scalar") ;

	assert(Ze_n_index == Al_extr_m2_index+1) ;
	assert(P_n_index == Ze_extr_m2_index+1) ;
	assert(Q_n_index == P_nm2_index+1) ;
	assert(mass_aspect_index == phi_nm2_index+1) ;

	return ;
}
/*===========================================================================*/
/* call to set global variables for field evolution */
/*===========================================================================*/
void set_globals(amr_grid *grid)
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
	Al_extr_m1 = grid->grid_funcs[Al_extr_m1_index] ;
	Al_extr_m2 = grid->grid_funcs[Al_extr_m2_index] ;

	Ze_n   = grid->grid_funcs[Ze_n_index  ] ;
	Ze_nm1 = grid->grid_funcs[Ze_nm1_index] ;
	Ze_nm2 = grid->grid_funcs[Ze_nm2_index] ;
	Ze_extr_m1 = grid->grid_funcs[Ze_extr_m1_index] ;
	Ze_extr_m2 = grid->grid_funcs[Ze_extr_m2_index] ;

	phi_n   = grid->grid_funcs[phi_n_index  ] ;
	phi_nm1 = grid->grid_funcs[phi_nm1_index] ;
	phi_nm2 = grid->grid_funcs[phi_nm2_index] ;

	mass_aspect = grid->grid_funcs[mass_aspect_index] ;

	ingoing_null_characteristic  = grid->grid_funcs[ingoing_null_characteristic_index] ;
	outgoing_null_characteristic = grid->grid_funcs[outgoing_null_characteristic_index] ;

	ingoing_scalar_characteristic  = grid->grid_funcs[ingoing_scalar_characteristic_index] ;
	outgoing_scalar_characteristic = grid->grid_funcs[outgoing_scalar_characteristic_index] ;

	Ricci_scalar = grid->grid_funcs[Ricci_scalar_index] ;
	Gauss_Bonnet_scalar = grid->grid_funcs[Gauss_Bonnet_scalar_index] ;

	eom_TR = grid->grid_funcs[eom_TR_index] ;
	eom_ThTh = grid->grid_funcs[eom_ThTh_index] ;
	eom_scalar = grid->grid_funcs[eom_scalar_index] ;

	Nx = grid->Nx ;

	dt = grid->dt ;
	dx = grid->dx ;

	grid_time = grid->time ;

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

	if (strcmp(initial_data,"r4Exp")==0) {
		initial_data_Gaussian(
			stereographic_L,
			Nx, 	
			dx,
			bbox,
			direction,
			amp, width, center,
			Al_n, Ze_n, 
			P_n, Q_n,
			phi_n, phi_nm1)
		;
	} else 
	if (strcmp(initial_data,"initial_black_hole")==0) {
		if ((grid->level)==0) {
			(grid->excised_jC) = (int)round(initial_black_hole_mass/(2.*dx)) ;
		} else {
			set_excision_point(grid) ;
		}
		initial_data_black_hole(
			stereographic_L,
			Nx, 	
			dx,
			bbox,
			initial_black_hole_mass,
			grid->excised_jC,
			Al_n, Ze_n)
		;
	} else {
		fprintf(stderr, "ERROR(set_free_initial_data): initial_data = %s\n", initial_data) ;
		exit(EXIT_FAILURE) ;
	}
	return ;
}
/*===========================================================================*/
/* rescale Al to be unity at spatial infinity at the coarsest level, then
 * synchronize all higher levels to this value */
/*===========================================================================*/
void rescale_Al(amr_grid *grid)
{
	if ((grid->level)==0) {
		double rescale_param = grid->grid_funcs[Al_n_index][Nx-1] ;
		for (int iC=0; iC<Nx; iC++) {
			grid->grid_funcs[Al_n_index][iC] /= rescale_param ;
		}
	} else {
		double rescale_param = 
			grid->grid_funcs[Al_n_index][Nx-1]
		/	grid->parent->grid_funcs[Al_n_index][grid->perim_coords[1]] 
		;
		for (int iC=0; iC<Nx; iC++) {
			grid->grid_funcs[Al_n_index][iC] /= rescale_param ;
		}
	}
	return ;
}
/*===========================================================================*/
/*===========================================================================*/
void solve_ode(amr_grid *grid)
{
	set_globals(grid) ;
	if (strcmp(theory, "massless_scalar_GR") == 0) {
		solve_Al_Ze_GR(
			stereographic_L,
			Nx,
			dt, 	dx,
			err_tolerance,
			excision_on,
			excised_jC,
			bbox,
			Al_n, 	Al_nm1, Ze_n, Ze_nm1,
			 P_n, 	 P_nm1,  Q_n,  Q_nm1)
		;
	} else 
	if (strcmp(theory,"EdGB")==0) { 
		solve_Al_Ze_EdGB(
			stereographic_L, coupling_gbs,
			Nx,
			dt, 	dx,
			err_tolerance, 
			excision_on,
			excised_jC,
			bbox,
			Al_n, Al_nm1, Ze_n, Ze_nm1,
			 P_n,  P_nm1,  Q_n,  Q_nm1)
		;
	} else {} ;
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void evolve_hyperbolic_pde(amr_grid *grid)
{
	set_globals(grid) ;
	rescale_Al(grid) ;
	if (strcmp(theory,"massless_scalar_GR")==0) {
		if (strcmp(solver_Ze,"constrained")==0) {
			advance_tStep_PQ_massless_scalar_GR(
				stereographic_L,
				Nx, dt, dx, 
				err_tolerance,
				excision_on,
				excised_jC,
				bbox, perim_interior,
				Al_n, Al_nm1, Ze_n, Ze_nm1,
				 P_n,  P_nm1,  Q_n,  Q_nm1
			) ;	
		} else
		if (strcmp(solver_Ze,"free")==0) {
			/* add routine later */
		} else {
			fprintf(stderr,"ERROR(main.c): solver_Ze=%s\n",solver_Ze) ;
			exit(EXIT_FAILURE) ;
		}
	} else 
	if (strcmp(theory,"EdGB")==0) {
		if (strcmp(solver_Ze,"constrained")==0) {
			advance_tStep_PQ_massless_scalar_EdGB(
				stereographic_L, coupling_gbs,
				Nx, 
				dt, dx, 
				err_tolerance, 
				excision_on,
				excised_jC,
				bbox, 
				perim_interior,
				Al_n, Al_nm1, Ze_n, Ze_nm1,
				 P_n,  P_nm1,  Q_n,  Q_nm1)
			;
		} else
		if (strcmp(solver_Ze,"free")==0) {
			/* add routine later */
		} else {
			fprintf(stderr,"ERROR(main.c): solver_Ze=%s\n",solver_Ze) ;
			exit(EXIT_FAILURE) ;
		}
	} else { ; }
	
	advance_tStep_phi(
		Nx, 
		dt, 
		Al_n, Al_nm1, Ze_n, Ze_nm1,
		 P_n,  P_nm1,  Q_n,  Q_nm1,
		phi_n, phi_nm1)
	;
	return ;
}
/*===========================================================================*/
/* independent residuals, curvature scalars, etc */
/*===========================================================================*/
void compute_diagnostics(amr_grid *grid)
{
	set_globals(grid) ;

	double x_lower = bbox[0] ;
/*--------------------------------------------------------------------------*/
	if (strcmp(theory,"massless_scalar_GR")==0) {
		compute_diagnostics_massless_scalar_GR(
			Nx, excised_jC, 
			stereographic_L,
			dt, dx,
			x_lower,
			Al_n, Al_nm1, Al_nm2,
			Ze_n, Ze_nm1, Ze_nm2,
			P_n, P_nm1, P_nm2,
			Q_n, 
			eom_TR,
			eom_ThTh,
			eom_scalar)
		;
	} else 
	if (strcmp(theory,"EdGB")==0) {
		compute_diagnostics_massless_scalar_EdGB(
			Nx, excised_jC, 
			stereographic_L,
			coupling_gbs,
			dt, dx,
			x_lower,
			Al_n, Al_nm1, Al_nm2,
			Ze_n, Ze_nm1, Ze_nm2,
			P_n, P_nm1, P_nm2,
			Q_n, 
			eom_TR,
			eom_ThTh,
			eom_scalar,
			ingoing_scalar_characteristic,
			outgoing_scalar_characteristic)
		;
	} else {
		fprintf(stderr,"ERROR(main.c): theory=%s",solver_Ze) ;
		exit(EXIT_FAILURE) ;
	}
/*--------------------------------------------------------------------------*/
/* ''_general: includes setting excised_jC */ 
/*--------------------------------------------------------------------------*/
	compute_diagnostics_general(
		Nx, 
		stereographic_L,
		dt, dx,
		x_lower,
		Al_n, Al_nm1, Al_nm2,
		Ze_n, Ze_nm1, Ze_nm2,
		&excised_jC,
		mass_aspect,
		ingoing_null_characteristic,
		outgoing_null_characteristic,
		Ricci_scalar,
		Gauss_Bonnet_scalar)
	;
	if ((grid->level)==0) { 
		grid->excised_jC = excised_jC ;
		printf("t\t%f\tMA\t%f\texc_jC\t%d\n", (grid->time), mass_aspect[Nx-1], excised_jC) ;
	}
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void save_to_file(amr_grid *grid)
{
	set_globals(grid) ;

	if (!made_output_files) {
	/* 	confirm that output directory exists
	*/
		DIR* dir = opendir(output_dir) ;
		if (dir!=NULL) {
			closedir(dir) ;
		} else 
		if (ENOENT==errno) {
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

		snprintf(output_name_phi, MAX_NAME_LEN, "%sphi.sdf", output_dir) ;

		snprintf(output_name_mass_aspect, MAX_NAME_LEN, "%smass_aspect.sdf", output_dir) ;

		snprintf(output_name_ingoing_null_characteristic,  MAX_NAME_LEN, "%singoing_null_characteristic.sdf", output_dir) ;
		snprintf(output_name_outgoing_null_characteristic, MAX_NAME_LEN, "%soutgoing_null_characteristic.sdf",  output_dir) ;
		snprintf(output_name_ingoing_scalar_characteristic,  MAX_NAME_LEN, "%singoing_scalar_characteristic.sdf", output_dir) ;
		snprintf(output_name_outgoing_scalar_characteristic, MAX_NAME_LEN, "%soutgoing_scalar_characteristic.sdf",  output_dir) ;
		snprintf(output_name_Ricci_scalar,        MAX_NAME_LEN, "%sRicci_scalar.sdf",        output_dir) ;
		snprintf(output_name_Gauss_Bonnet_scalar, MAX_NAME_LEN, "%sGauss_Bonnet_scalar.sdf", output_dir) ;

		snprintf(output_name_eom_TR, MAX_NAME_LEN, "%seom_TR.sdf", output_dir) ;
		snprintf(output_name_eom_ThTh, MAX_NAME_LEN, "%seom_ThTh.sdf", output_dir) ;
		snprintf(output_name_eom_scalar, MAX_NAME_LEN, "%seom_scalar.sdf", output_dir) ;

		snprintf(output_name_test_1, MAX_NAME_LEN, "%stest_1.sdf", output_dir) ;
		snprintf(output_name_test_2, MAX_NAME_LEN, "%stest_2.sdf", output_dir) ;
	}
	made_output_files = true ;
/* 	save to file-see man pages for bbhutil utilities on Choptuik's website 
*/
	int rank = 1 ;

	gft_out_bbox(output_name_Al, grid_time, &Nx, rank, bbox, Al_n) ;
	gft_out_bbox(output_name_Ze, grid_time, &Nx, rank, bbox, Ze_n) ;
	gft_out_bbox(output_name_P,  grid_time, &Nx, rank, bbox,  P_n) ;
	gft_out_bbox(output_name_Q,  grid_time, &Nx, rank, bbox,  Q_n) ;

	double *test_1 = calloc(Nx,sizeof(double)) ;
	double *test_2 = calloc(Nx,sizeof(double)) ;

	for (int jC=0; jC<Nx; jC++) {
		test_1[jC] = P_nm1[jC] ;
		test_2[jC] = P_nm2[jC] ;
	}

	gft_out_bbox(output_name_test_1, grid_time, &Nx, rank, bbox, test_1) ;
	gft_out_bbox(output_name_test_2, grid_time, &Nx, rank, bbox, test_2) ;

	free(test_1) ;
	free(test_2) ;

	gft_out_bbox(output_name_phi,  grid_time, &Nx, rank, bbox, phi_n) ;

	gft_out_bbox(output_name_mass_aspect,  grid_time, &Nx, rank, bbox, mass_aspect) ;

	gft_out_bbox(output_name_ingoing_null_characteristic,  grid_time, &Nx, rank, bbox, ingoing_null_characteristic) ;
	gft_out_bbox(output_name_outgoing_null_characteristic, grid_time, &Nx, rank, bbox, outgoing_null_characteristic) ;
	gft_out_bbox(output_name_ingoing_scalar_characteristic,  grid_time, &Nx, rank, bbox, ingoing_scalar_characteristic) ;
	gft_out_bbox(output_name_outgoing_scalar_characteristic, grid_time, &Nx, rank, bbox, outgoing_scalar_characteristic) ;
	gft_out_bbox(output_name_Ricci_scalar,        grid_time, &Nx, rank, bbox,  Ricci_scalar) ;
	gft_out_bbox(output_name_Gauss_Bonnet_scalar, grid_time, &Nx, rank, bbox,  Gauss_Bonnet_scalar) ;

	gft_out_bbox(output_name_eom_TR, grid_time, &Nx, rank, bbox, eom_TR) ;
	gft_out_bbox(output_name_eom_ThTh, grid_time, &Nx, rank, bbox, eom_ThTh) ;
	gft_out_bbox(output_name_eom_scalar, grid_time, &Nx, rank, bbox, eom_scalar) ;

	fflush(NULL) ;
	
	return ;
}
/*===========================================================================*/
/* call amr evolution */
/*===========================================================================*/

int main(void) 
{
	get_run_data(
		&theory,
		&output_dir,
		&solver_Ze,
		&Nx, &Nt, &t_step_save,
		&cfl_num, 
		&bbox[0], &bbox[1],
		&stereographic_L,
		&coupling_gbs,
		&dt, &dx,
		&err_tolerance) 
	; 
	get_initial_data(
		&initial_data, &direction, &amp, &width, &center, 
		&initial_black_hole_mass) 
	; 
	amr_field *fields = set_fields() ; 
	find_field_indices(fields) ; 

	amr_grid_hierarchy *gh 
	= amr_init_grid_hierarchy(
		fields,
        	Nt, Nx, t_step_save,
		cfl_num,
		bbox,
		excision_on)
	;
	clock_t begin_time = clock() ;
	amr_main(
		gh, 
		set_free_initial_data,
		evolve_hyperbolic_pde,
		solve_ode,
		compute_diagnostics,
		save_to_file)
	;
	clock_t end_time = clock() ;
	double time_spent = (end_time-begin_time)/CLOCKS_PER_SEC ;
	printf("time evolving sim (sec): %f\n", time_spent) ;

	amr_destroy_grid_hierarchy(gh) ;
/*--------------------------------------------------------------------------*/
/* theory and initial_data are malloc'ed in get_run_data and get_initial_data,
 * respectively. */ 
/*--------------------------------------------------------------------------*/
	free(theory) ;
	free(initial_data) ;
/*--------------------------------------------------------------------------*/
	return EXIT_SUCCESS ;
}
