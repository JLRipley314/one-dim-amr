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
	int num_time_level;  /* number of time-levels (in AMR hierarchy), from 1 .. num_time_level */
}
;
/*============================================================================*/
/*	each grid holds storage for ``grid functions'': i.e. the fields
	at each refinement level */
/*============================================================================*/
struct grid
{
	struct grid* child  ; /* to finer grid */
	struct grid* parent ; /* to coarser grid */

	int dim ; /* number of grid points */

	double bbox[2] ; /* physical coordinates bounding the grid */

	double time ;
	double perimeter_coords[2] ; /* coordinates with respect to parent grid */
	
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
	double dx ;          /* discretization for each dimension */

	struct grid* grids ;
}
;
/*============================================================================*/
/* contains all levels */
/*============================================================================*/
struct comp_grid_hierarchy 
{
	/* the following define the structure of the grid-hierarchy */

	int dim;                        /* spatial dimension */
	int refinement_ratio;           /* refine space and time steps the same amount and the same for each level */
	double dx[AMR_MAX_LEVELS];       /* level defined by discretization scale dx of first coordinate */
	double dt[AMR_MAX_LEVELS];
	double shape;                    /* geometry of base level   */
	double  bbox[2];
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
	bool excision_on ;
/*	user defined function to mask excised regions */
	void (*app_fill_exc_mask)(double* mask, int dim, int* shape, double* bbox, double excised) ;
	double excised ; /* value denoting excised grid point */
	int excised_jC ; /* value denoting excised grid point */
	int amr_mask_grid_func_num ; /* amr mask grid function numbers */
}
;
/*============================================================================*/

#endif /*_GRID_HIERARCHY_H_*/
