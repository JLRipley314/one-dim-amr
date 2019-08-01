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
#include <stdbool.h>

extern const int amr_max_levels; /* maximum number of grid levels */
extern const int refinement ; /* spatial and temporal refinement scale */
extern const int regrid ; /* regrid every regrid time steps */
extern const int buffer_coord ; /* regrid every regrid time steps */

extern const double trunc_err_tolerance ; /* when above this flag for finer grid */

extern const char HYPERBOLIC[] ;
extern const char ELLIPTIC[] ;
extern const char ODE[] ;
extern const char DIAGNOSTIC[] ;

/*============================================================================*/
/* 	keeps track of where variable is located in storage on each grid,
	along with instructions for how to inject and interpolate */ 
/*============================================================================*/
typedef struct amr_field 
{
	struct amr_field *prev;
	struct amr_field *next;
	char *name;    
	int index;           /* index of 0th time level */
	int time_levels; /* number of time-levels (in AMR hierarchy), from 1 .. num_time_level */
	int extrap_levels; /* for extrapolation (for ode/ELLIPTIC equations) */
	char* pde_type;      /* either HYPERBOLIC or ELLIPTIC */
	int flagged_coords[2] ; /* upper (0) and lower (1) boundary of flagged region for adding finer grid */

} amr_field 
;
/*============================================================================*/
/*	each grid holds storage for ``grid functions'': i.e. the fields
	at each refinement level */
/*============================================================================*/
typedef struct amr_grid
{
	struct amr_grid* child ; 	/* to finer grid */
	struct amr_grid* parent ;	/* to coarser grid */

	int level ;
	int Nx ; /* number of grid points */

	int excised_jC ; /* excision point in this grids coordinates (if 0 then no excision in this grid) */
	bool excision_on ;
	
	double time ;
	double dt, dx ;
	int tC ;

	double bbox[2] ; /* physical coordinates bounding the grid */

	int flagged_coords[2] ; /* for regridding-new coords for finer (child) grid */
	int perim_coords[2] ; /* coordinates with respect to parent grid */
	
	int num_grid_funcs;	/* total number of grid functions */
	int time_levels;	/* total number of grid functions */
	double** grid_funcs ;	/* pointer to array of pointer to grid function data */

	bool perim_interior[2] ; /* true: then grid boundaries not physical boundaries; else false */
} amr_grid
;
/*============================================================================*/
/* contains base grid, which then points to all subsequent grids */
/*============================================================================*/
typedef struct amr_grid_hierarchy 
{
	/* the following define the structure of the grid-hierarchy */

	double cfl_num;		/* courant factor */
	int time_levels;	/* total number of grid functions */
	int t_step_save ;
	int Nt ;
	struct amr_field* fields; /* pointer to a linked list of (num_vars) variable structures */

	int* field_indices ; /* labels grid function idex that corresponds to a field */ 

/* 	big difference from Frans' code:  we do not have a ``context'' for different
	instances of grid hieracrhies as our application is simpler (1+1D, no parallel)
*/
	struct amr_grid* grids ;

	bool excision_on ;
/*	user defined function to mask excised regions */
	void (*app_fill_exc_mask)(double* mask, int dim, int* shape, double* bbox, double excised) ;
} amr_grid_hierarchy
;
/*============================================================================*/
/* return 0 then no errors */
/*============================================================================*/
void amr_set_to_tail(amr_grid** grid) ; 
void amr_set_to_head(amr_grid** grid) ;

void amr_find_grid(int level, amr_grid_hierarchy* gh, amr_grid* grid) ;

void amr_add_finer_grid(int left_coord, int right_coord, amr_grid* parent) ;

void amr_insert_grid(amr_grid *grid_to_insert, amr_grid* grid) ;

void amr_destroy_grid(amr_grid* grid) ;

void flag_field_regridding_coords(amr_field *fields, amr_grid *parent, amr_grid *grid) ;

amr_grid_hierarchy* amr_init_grid_hierarchy(
	amr_field* fields,
	int Nt, int Nx, int t_step_save,
	double cfl_num,
	double bbox[2],
	bool excision_on)
;
void inject_overlaping_fields(amr_field *fields, amr_grid *grid, amr_grid *parent) ;

void regrid_all_finer_levels(amr_field *fields, amr_grid *base_grid) ;

void add_initial_grids(amr_grid_hierarchy* gh) ;

void amr_destroy_grid_hierarchy(amr_grid_hierarchy* gh) ;

void amr_add_field(amr_field* field, char* name, char* pde_type, int time_levels, int extrap_levels) ;

amr_field* amr_init_fields(char* name, char* pde_type, int time_levels, int extrap_levels) ;

int amr_return_field_index(amr_field* fields, char* name) ; 

 
#endif /*_GRID_HIERARCHY_H_*/
