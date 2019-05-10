#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "initial_data.h"
#include "PDE_solvers.h"

#define MAX_FILE_NAME 1024
#define OUTPUT_DIR "/home/jripley/one-dim-amr/examples/wave/output/"

int main(int argc, char* argv[])
{
	int Nx = pow(2,10)+1 ;
	int Nt = pow(2,10)+1 ;
	int tss = 4 ;

	char output_file_P[MAX_FILE_NAME+1] ;
	char output_file_Q[MAX_FILE_NAME+1] ;

	snprintf(output_file_P, MAX_FILE_NAME, "%sP.txt", OUTPUT_DIR) ;
	snprintf(output_file_Q, MAX_FILE_NAME, "%sQ.txt", OUTPUT_DIR) ;

	double* P_n   = calloc(Nx,sizeof(double)) ;
	double* P_nm1 = calloc(Nx,sizeof(double)) ;
	double* Q_n   = calloc(Nx,sizeof(double)) ;
	double* Q_nm1 = calloc(Nx,sizeof(double)) ;

	set_initial_data(P_n, Q_n) ;
	copy_array(Nx, P_n, P_nm1) ;
	copy_array(Nx, Q_n, Q_nm1) ;

	for (int tC=0; iC<Nt; tC++) {
		advance_tStep_wave(Nx, dt, dx, P_n, P_nm1, Q_n, Q_nm1) ;

		Kreiss_Oliger_Filer(P_n) ;
		Kreiss_Oliger_Filer(Q_n) ;

		copy_data(Nx, P_n, P_nm1) ;
		copy_data(Nx, Q_n, Q_nm1) ;

		if (tC%tss == 0) {
			save_to_file(output_file_P,P_n) ;
			save_to_file(output_file_Q,Q_n) ;
		}
	}	

	free(P_n) ;
	free(P_nm1) ;
	free(Q_n) ;
	free(Q_nm1) ;

	return EXIT_SUCCESS ;
}
