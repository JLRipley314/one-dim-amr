#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

#include "file_io.h"

#define HOME_DIR "/home/jripley/one-dim-amr/examples/massless_scalar_collapse/" 

/*==========================================================================*/
static void get_int_val(char* line_val, char* token, char* comparison, char* delimeter, int* val)
{
	if ((token != NULL) && (strcmp(token,comparison)==0)) {
		*val = strtod(line_val, NULL) ;
		printf("%s %d\n", comparison, *val) ;
	}
	return ;
}
/*==========================================================================*/
static void get_double_val(char* line_val, char* token, char* comparison, char* delimeter, double* val)
{
	if ((token != NULL) && (strcmp(token,comparison)==0)) {
		*val = strtod(line_val, NULL) ;
		printf("%s %f\n", comparison, *val) ;
	}
	return ;
}
/*==========================================================================*/
static void get_string_val(char* line_val, char* token, char* comparison, char* delimeter, char** val)
{
	if ((token != NULL) && (strcmp(token,comparison)==0)) {
		char* inter_val = strsep(&line_val, "\n") ;
		*val = malloc(strlen(inter_val)+1) ;
		strcpy(*val,inter_val) ;
		printf("%s %s\n", comparison, *val) ;
	}
	return ;
}
/*==========================================================================*/
static void check_if_file_exists(FILE* file, char* file_name) 
{
	if (file==NULL) {
		printf("ERROR(file_io.c): cannot find %s\n", file_name) ;
		exit(EXIT_FAILURE) ;
	}
	return ;
}
/*==========================================================================*/
void get_run_data(
	char **theory,
	char **output_dir,
	char **solver_Ze,
	int *Nx, int *Nt, int *t_step_save,
	double *cfl_num, 
	double *bbox_0, double *bbox_1,
	double *stereographic_L,
	double *coupling_gbs,
	double *dt, double *dx, 
	double *err_tolerance)
{
	char run_data_name[MAX_NAME_LEN+1] ;
	snprintf(run_data_name, MAX_NAME_LEN, "%srun_data.txt", HOME_DIR) ;
	FILE *run_data_file = fopen(run_data_name, "r") ;
	check_if_file_exists(run_data_file, run_data_name) ; 

        char *line  = NULL ;
        char *token = NULL ;
        size_t len = 0 ; 
        ssize_t read = 0 ; 

	char *delimeter = "=" ;
	
        while ((read = getline(&line, &len, run_data_file)) > 0) {
		char *line_val = line ;
                token = strsep(&line_val, delimeter) ;

		get_string_val(line_val, token, "theory", delimeter, theory) ;
		get_string_val(line_val, token, "output_dir", delimeter, output_dir) ;
		get_string_val(line_val, token, "solver_Ze", delimeter, solver_Ze) ;

		get_int_val(line_val, token, "Nx", delimeter, Nx) ;
		get_int_val(line_val, token, "Nt", delimeter, Nt) ;
		get_int_val(line_val, token, "t_step_save", delimeter, t_step_save) ;

		get_double_val(line_val, token, "cfl_num", delimeter, cfl_num) ;
		get_double_val(line_val, token, "dx", delimeter, dx) ;
		get_double_val(line_val, token, "dt", delimeter, dt) ;

		get_double_val(line_val, token, "stereographic_L", delimeter, stereographic_L) ;

		get_double_val(line_val, token, "coupling_gbs", delimeter, coupling_gbs) ;

		get_double_val(line_val, token, "err_tolerance", delimeter, err_tolerance) ;
	}
/*	
	this is for the coarsest grid so the bounding box (bbox)
	is from r=0 to the stereographic radius (r=infty) 
*/
	*bbox_0 = 0 ;
	*bbox_1 = *stereographic_L ;

	free(line) ;
	fclose(run_data_file) ;

	return ;
}
/*==========================================================================*/
void get_initial_data(
	char **initial_data,
	char **direction,
	double *amp, double *width, double *center,
	double *initial_black_hole_mass) 
{
	char initial_data_name[MAX_NAME_LEN+1] ;
	snprintf(initial_data_name, MAX_NAME_LEN, "%sinitial_data.txt", HOME_DIR) ;
	FILE *initial_data_file = fopen(initial_data_name, "r") ;
	check_if_file_exists(initial_data_file, initial_data_name) ; 

        char *line  = NULL ;
        char *token = NULL ;
        size_t len = 0 ; 
        ssize_t read = 0 ; 

	char *delimeter = "=" ;

        while ((read = getline(&line, &len, initial_data_file)) > 0) {
		char *line_val = line ;
                token = strsep(&line_val, delimeter) ;

		get_string_val(line_val, token, "initial_data", delimeter, initial_data) ;
		get_string_val(line_val, token, "direction", delimeter, direction) ;

		get_double_val(line_val, token, "amp", delimeter, amp) ;
		get_double_val(line_val, token, "width", delimeter, width) ;
		get_double_val(line_val, token, "center", delimeter, center) ;

		get_double_val(line_val, token, "initial_black_hole_mass", delimeter, initial_black_hole_mass) ;

	}
	free(line) ;
	fclose(initial_data_file) ;

	return ;
}
