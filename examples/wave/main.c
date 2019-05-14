#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "amr_evolve.h"
#include "amr_grid_hierarchy.h"

#include "evolution_routines.h"

#define MAX_FILE_NAME 1024
#define OUTPUT_DIR "/home/jripley/one-dim-amr/examples/wave/output/"

/*===========================================================================*/
/* global variables for evolution-convenient for function calls */
/*===========================================================================*/
double* P_n ;
double* P_nm1 ;
double* Q_n ;
double* Q_nm1 ;
double dx, dt ;
double bbox[2] ;
int Nx ;
int excised_jC ;
int num_grid_funcs = 4 ;
int P_n_index, P_nm1_index, Q_n_index, Q_nm1_index ;
int perim_coords[2] ;
bool perim_interior[2] ;

char output_name_P[MAX_FILE_NAME+1] ;
char output_name_Q[MAX_FILE_NAME+1] ;

FILE* output_file_P ;
FILE* output_file_Q ;

bool made_files = false ;
/*===========================================================================*/
/* call after variables have been defined */
/*===========================================================================*/
void set_fields_index(void)
{
	P_n_index   = 0 ;
	P_nm1_index = 1 ;
	Q_n_index   = 2 ;
	Q_nm1_index = 3 ;

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

	Nx = grid->Nx ;

	bbox[0] = grid->bbox[0] ;
	bbox[1] = grid->bbox[1] ;

	perim_interior[0] = grid->perim_interior[0] ;
	perim_interior[1] = grid->perim_interior[1] ;
	
	perim_coords[0] = grid->perim_coords[0] ;
	perim_coords[1] = grid->perim_coords[1] ;

	dt = grid->dt ;
	dx = grid->dx ;

	excised_jC = grid->excised_jC ;

	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void initial_data(struct amr_grid* grid)
{
	set_globals(grid) ;

	initial_data_Gaussian(Nx, dx, bbox[0], P_n, P_nm1, Q_n, Q_nm1) ;
	
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void wave_evolve(struct amr_grid* grid)
{
	set_globals(grid) ;

	advance_tStep_wave(Nx, dt, dx, perim_interior, P_n, P_nm1, Q_n, Q_nm1) ;
	
	return ;
}
/*===========================================================================*/
/* computes one time step (to tolerance) of wave equation */
/*===========================================================================*/
void wave_to_file_output(struct amr_grid* grid)
{
	set_globals(grid) ;

	if (made_files == false) {
		FILE* output_file_P = fopen(output_name_P, "w") ;
		if (output_file_P == NULL ) {
			printf("ERROR(main.c): output_file_P == NULL\n") ;
			exit(EXIT_FAILURE) ;
		}
		FILE* output_file_Q = fopen(output_name_Q, "w") ;
		if (output_file_Q == NULL ) {
			printf("ERROR(main.c): output_file_Q == NULL\n") ;
			exit(EXIT_FAILURE) ;
		}	
		made_files = true ;
	} else {
		FILE* output_file_P = fopen(output_name_P, "a") ;
		if (output_file_P == NULL ) {
			printf("ERROR(main.c): output_file_P == NULL\n") ;
			exit(EXIT_FAILURE) ;
		}
		FILE* output_file_Q = fopen(output_name_Q, "a") ;
		if (output_file_Q == NULL ) {
			printf("ERROR(main.c): output_file_Q == NULL\n") ;
			exit(EXIT_FAILURE) ;
		}
	}
	save_to_txt_file(Nx, output_file_P, P_n) ;
	save_to_txt_file(Nx, output_file_Q, Q_n) ;

	fclose(output_file_P) ;
	fclose(output_file_Q) ;

	output_file_P = NULL ;
	output_file_Q = NULL ;
	
	return ;
}
/*===========================================================================*/
/* call amr evolution */
/*===========================================================================*/
int main(int argc, char* argv[])
{
	int Nx = pow(2,9)+1 ;
	int Nt = pow(2,12)+1 ;
	int tss = 4 ;

	bool perim_interior[2] = {false,false} ;

	double cfl_num = 0.25 ;

	double bbox[2] = {-50,50} ;

	double dx = (bbox[1] - bbox[0]) / (Nx-1) ;
	double dt = cfl_num * dx ;

	snprintf(output_name_P, MAX_FILE_NAME, "%sP.txt", OUTPUT_DIR) ;
	snprintf(output_name_Q, MAX_FILE_NAME, "%sQ.txt", OUTPUT_DIR) ;

	FILE* output_file_P = fopen(output_name_P, "w") ;
	if (output_file_P == NULL ) {
		printf("ERROR(main.c): output_file_P == NULL\n") ;
		exit(EXIT_FAILURE) ;
	}
	FILE* output_file_Q = fopen(output_name_Q, "w") ;
	if (output_file_Q == NULL ) {
		printf("ERROR(main.c): output_file_Q == NULL\n") ;
		exit(EXIT_FAILURE) ;
	}

	double* P_n   = calloc(Nx,sizeof(double)) ;
	double* P_nm1 = calloc(Nx,sizeof(double)) ;
	double* Q_n   = calloc(Nx,sizeof(double)) ;
	double* Q_nm1 = calloc(Nx,sizeof(double)) ;

	initial_data_Gaussian(Nx, dx, bbox[0], P_n, P_nm1, Q_n, Q_nm1) ;
	save_to_txt_file(Nx, output_file_P, P_n) ;
	save_to_txt_file(Nx, output_file_Q, Q_n) ;
	fflush(NULL) ;

	for (int tC=1; tC<Nt; tC++) {
		advance_tStep_wave(Nx, dt, dx, perim_interior, P_n, P_nm1, Q_n, Q_nm1) ;
		if (tC%tss == 0) {
			save_to_txt_file(Nx, output_file_P, P_n) ;
			save_to_txt_file(Nx, output_file_Q, Q_n) ;
			fflush(NULL) ;
		}
	}	

	free(P_n) ;
	free(P_nm1) ;
	free(Q_n) ;
	free(Q_nm1) ;

	fclose(output_file_P) ;
	fclose(output_file_Q) ;

	return EXIT_SUCCESS ;
}
