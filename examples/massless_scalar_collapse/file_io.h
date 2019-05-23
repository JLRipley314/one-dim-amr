#ifndef _FILE_IO_H_
#define _FILE_IO_H_

#include <stdbool.h>

/*==========================================================================*/

#define MAX_NAME_LEN 1024
#define MAX_LINE_LEN 4096

typedef struct run_data
{
	char* theory ;

	int Nx ;
	int Nt ;

	int t_step_save ;

	double cfl_num ;

	double dx ;
	double dt ;

	double bbox[2] ;

	bool perim_interior[2] ;

	double stereographic_L ;

} run_data
;
/*==========================================================================*/
typedef struct initial_data
{
	char* type ;

	double amp ;
	double width ;
	double center ;	

} initial_data
;
/*==========================================================================*/
run_data* init_run_data(void) ; 
/*==========================================================================*/
initial_data* init_initial_data(void) ;
/*==========================================================================*/
void free_run_data(run_data* rd) ;
/*==========================================================================*/
void free_initial_data(initial_data* rd) ;
/*==========================================================================*/
#endif
