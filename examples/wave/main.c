#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grid_hierarchy.h"

#include "evolution_routines.h"

#define MAX_FILE_NAME 1024
#define OUTPUT_DIR "/home/jripley/one-dim-amr/examples/wave/output/"

int main(int argc, char* argv[])
{
	int Nx = pow(2,9)+1 ;
	int Nt = pow(2,12)+1 ;
	int tss = 4 ;

	double cfl_num = 0.25 ;

	double bbox[2] = {-50,50} ;

	double dx = (bbox[1] - bbox[0]) / (Nx-1) ;
	double dt = cfl_num * dx ;

	char output_name_P[MAX_FILE_NAME+1] ;
	char output_name_Q[MAX_FILE_NAME+1] ;

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

	initial_Data(Nx, dx, bbox[0], P_n, Q_n) ;
	copy_to_2nd_array(Nx, P_n, P_nm1) ;
	copy_to_2nd_array(Nx, Q_n, Q_nm1) ;
	save_to_txt_file(Nx, output_file_P, P_n) ;
	save_to_txt_file(Nx, output_file_Q, Q_n) ;
	fflush(NULL) ;

	for (int tC=1; tC<Nt; tC++) {
		advance_tStep_wave(Nx, dt, dx, P_n, P_nm1, Q_n, Q_nm1) ;

		copy_to_2nd_array(Nx, P_n, P_nm1) ;
		copy_to_2nd_array(Nx, Q_n, Q_nm1) ;

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
