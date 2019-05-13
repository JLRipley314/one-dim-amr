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
#define REFINEMENT 4
#define REGRID 8

#include <stdbool.h>

/*============================================================================*/
/* 	keeps track of where variable is located in storage on each grid,
	along with instructions for how to inject and interpolate */ 
/*============================================================================*/
struct amr_var 
{
	struct amr_var *next;
	char *name;    
	int in_amr_hier;     /* whether the variable exists in the AMR hierarchy or not */
	int num_time_level;  /* number of time-levels (in AMR hierarchy), from 1 .. num_time_level */
}
;
/*============================================================================*/
/*	each grid holds storage for ``grid functions'': i.e. the fields
	at each refinement level */
/*============================================================================*/
struct amr_grid
{
	struct amr_grid* child  ; /* to finer grid */
	struct amr_grid* parent ; /* to coarser grid */
	struct amr_grid* sibling ; /* to coarser grid */

	int level ;
	int Nx ; /* number of grid points */
	
	double time ;
	int tC ;
	double dt, dx ;


	double bbox[2] ; /* physical coordinates bounding the grid */

	double perim_coord[2] ; /* coordinates with respect to parent grid */
	
	int num_grid_funcs  ;
	double** grid_funcs ; /* pointer to array of pointer to grid function data */

	bool perim_interior[2] ;
}
;
/*============================================================================*/
/* contains all levels */
/*============================================================================*/
struct amr_grid_hierarchy 
{
	/* the following define the structure of the grid-hierarchy */

	int dim;                        /* spatial dimension */
	double dx[AMR_MAX_LEVELS];       /* level defined by discretization scale dx of first coordinate */
	double dt[AMR_MAX_LEVELS];
	double shape;                    /* geometry of base level   */
	double  bbox[2];
	double* seq_grid_hier_bboxes[AMR_MAX_LEVELS]; /* a copy of the sequential, AMR grid-hierarchy */
	int num_seq_grid_hier_bboxes[AMR_MAX_LEVELS];
	double cfl_num;                  /* courant factor */
	int num_vars;                    /* total number of variables */
	struct amr_var* vars;            /* pointer to a linked list of (num_vars) variable structures */
	int num_grid_funcs;              /* total number of grid functions */

/* 	big difference from Frans' code:  we do not have a ``context'' for different
	instances of grid hieracrhies as our application is simpler (1+1D, no parallel)
*/
	int min_level ;
	int max_level ;
	struct amr_grid* grid ;
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
/* return 0 then no errors */
/*============================================================================*/
int amr_get_grid_funcs(struct amr_grid* grid, double** grid_funcs) ;

int amr_find_grid(int level, struct amr_grid_hierarchy* gh, struct amr_grid* grid) ;

int amr_add_grid(struct amr_grid* grid) ;

#endif /*_GRID_HIERARCHY_H_*/
