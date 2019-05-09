#ifndef _GRID_HIERARCHY_H_
#define _GRID_HIERARCHY_H_
/*============================================================================*/
/* 
	grid hierarchy structure definitions and basic functions for
	searching through the hierarchy.

	Much of this is a simplified version of what may be found in
	the gh.h and gh.c codes 

	In particular, we assume the fields are grid centered, and we
	do not add anything for multigrid
*/ 
/*============================================================================*/


#define AMR_MAX_DIM 1 /* one dim amr...at least for now*/ 
#define AMR_MAX_LEVELS 2 

#include <stdbool.h>

/*============================================================================*/
/* 	keeps track of where variable is located in storage on each grid,
	along with instructions for how to inject and interpolate */ 
/*============================================================================*/
struct var 
{
	struct var *next;
	char *name;    
	int in_amr_hier;     /* whether the variable exists in the AMR hierarchy or not */
	int status;          /* one of STATUS_... flages above */
	int num_time_level;    /* number of time-levels (in AMR hierarchy), from 1 .. num_time_level */
	int start_grid_func; /* starting grid-function number */

	int amr_inject;        /* injection operator to use within AMR hiearchy */
	int amr_interp;        /*  interpolation " */
	int amr_bndry_interp;  /*  interpolation " when setting AMR boundary conditions */

	int regrid_transfer; /* wether to transfer after a regrid */
	int phys_bdy_type[2*AMR_MAX_DIM]; /* to guide the interpolation */
	int c_to_v;          /* is variable involved in c_to_v - V2 */
	int v_to_c;          /* is variable involved in v_to_c - V2 */
			/* added the above 2 to the end of the data structure for easier cp versioning */   
}
;
/*============================================================================*/
/*	each grid holds storage for ``grid functions'': i.e. the fields
	at each refinement level */
/*============================================================================*/
struct grid
{
	struct grid* next ; /* to finer grid */
	struct grid* prev ; /* to coarser grid */

	int dim ; /* number of grid points */

	double bbox[2*AMR_MAX_DIM] ; /* physical coordinates bounding the grid */

	double time ;
	double* perimeter_coords[AMR_MAX_DIM] ; /* pointer to perimeter coordinate arrays */
	
	int num_grid_funcs ;
	double** grid_funcs ; /* pointer to array of pointer to grid function data */
}
;
/*============================================================================*/
/*	pointers to multiple grids at this refinement */
/*============================================================================*/
struct level
{
	double time, dt ;
	double dx[AMR_MAX_DIM] ; /* discretization for each dimension */

	struct grid* grids ;
}
;
/*============================================================================*/
/* contains all levels */
/*============================================================================*/
struct comp_grid_hierarchy 
{
	/* the following define the structure of the grid-hierarchy */

	int dim;                         /* spatial dimension */
	int rho_sp[AMR_MAX_LEVELS];       /* spatial refinement ratio --- same for all dimensions  */
	int rho_tm[AMR_MAX_LEVELS];       /* temporal refinement ratio */
	double dx[AMR_MAX_LEVELS][AMR_MAX_DIM]; /* level defined by discretization scale dx of first coordinate */
	double dt[AMR_MAX_LEVELS];
	double shape[AMR_MAX_DIM];        /* geometry of base level   */
	double  bbox[2*AMR_MAX_DIM];
	double* seq_grid_hier_bboxes[AMR_MAX_LEVELS]; /* a copy of the sequential, AMR grid-hierarchy */
	int num_seq_grid_hier_bboxes[AMR_MAX_LEVELS];
	double cfl_num;                  /* courant factor */
	int num_vars;                    /* total number of variables */
	struct var* vars;                /* pointer to a linked list of (num_vars) variable structures */
	int num_grid_funcs;              /* total number of grid functions */

/* 	big difference from Frans' code:  we do not have a ``context'' for different
	instances of grid hieracrhies as our application is simpler (1+1D, no parallel)
*/
	int min_level ;
	int max_level ;
	struct level* levels[AMR_MAX_LEVELS] ;
/* 	excision functions */
	int excision_on ;
/*	user defined function to mask excised regions */
	void (*app_fill_exc_mask)(double* mask, int dim, int* shape, double* bbox, double excised) ;
	double excised ; /* value denoting excised grid point */
	int amr_mask_grid_func_num ; /* amr mask grid function numbers */
}
;
/*============================================================================*/



#endif /*_GRID_HIERARCHY_H_*/
