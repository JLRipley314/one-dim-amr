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
		*val = strsep(&line_val, "\n") ;
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
run_data* init_run_data(void) 
{
	run_data* rd = malloc(sizeof(run_data)) ;

	char run_data_name[MAX_NAME_LEN+1] ;
	snprintf(run_data_name, MAX_NAME_LEN, "%srun_data.txt", HOME_DIR) ;
	FILE* run_data_file = fopen(run_data_name, "r") ;
	check_if_file_exists(run_data_file, run_data_name) ; 

        char* line  = NULL ;
        char* token = NULL ;
        size_t len = 0 ; 
        ssize_t read = 0 ; 

	char* delimeter = ":" ;

        while ((read = getline(&line, &len, run_data_file)) > 0) {
		char* line_val = line ;
                token = strsep(&line_val, delimeter) ;

		get_string_val(line_val, token, "theory", delimeter, &(rd->theory)) ;

		get_int_val(line_val, token, "Nx", delimeter, &(rd->Nx)) ;
		get_int_val(line_val, token, "Nt", delimeter, &(rd->Nt)) ;
		get_int_val(line_val, token, "t_step_save", delimeter, &(rd->t_step_save)) ;

		get_double_val(line_val, token, "bbox[0]", delimeter, &(rd->bbox[0])) ;
		get_double_val(line_val, token, "bbox[1]", delimeter, &(rd->bbox[1])) ;

		get_double_val(line_val, token, "cfl_num", delimeter, &(rd->cfl_num)) ;
		get_double_val(line_val, token, "dx", delimeter, &(rd->dx)) ;
		get_double_val(line_val, token, "dt", delimeter, &(rd->dt)) ;

		get_double_val(line_val, token, "stereographic_L", delimeter, &(rd->stereographic_L)) ;

	}
	free(line) ;
	fclose(run_data_file) ;

	return rd ;
}
/*==========================================================================*/
initial_data* init_initial_data(void) 
{
	initial_data* id = malloc(sizeof(initial_data)) ;

	char initial_data_name[MAX_NAME_LEN+1] ;
	snprintf(initial_data_name, MAX_NAME_LEN, "%sinitial_data.txt", HOME_DIR) ;
	FILE* initial_data_file = fopen(initial_data_name, "r") ;
	check_if_file_exists(initial_data_file, initial_data_name) ; 

        char* line  = NULL ;
        char* token = NULL ;
        size_t len = 0 ; 
        ssize_t read = 0 ; 

	char* delimeter = ":" ;

        while ((read = getline(&line, &len, initial_data_file)) > 0) {
		char* line_val = line ;
                token = strsep(&line_val, delimeter) ;

		get_string_val(line_val, token, "type", delimeter, &(id->type)) ;

		get_double_val(line_val, token, "amp", delimeter, &(id->amp)) ;
		get_double_val(line_val, token, "width", delimeter, &(id->width)) ;
		get_double_val(line_val, token, "center", delimeter, &(id->center)) ;

	}
	free(line) ;
	fclose(initial_data_file) ;

	return id ;
}
