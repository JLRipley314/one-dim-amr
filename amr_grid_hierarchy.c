/*===========================================================================*/
/*===========================================================================*/
#include "amr_grid_hierarchy.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*============================================================================*/
const int amr_max_levels = 6 ; 
const int refinement = 2 ; 
const int regrid = 64 ; 
const int buffer_coord = 64 ; 
const int min_grid_size = 32 ;

const double trunc_err_tolerance = 1e-5 ; 

const char HYPERBOLIC[] = "hyperbolic" ;
const char ELLIPTIC[] = "elliptic" ;
const char ODE[] = "ode" ;
const char DIAGNOSTIC[] = "diagnostic" ;
/*============================================================================*/
/* arrays called as array[row][column] */
/*============================================================================*/
double** allocate_double_2DArray(int rows, int columns, double initializeValue)  
{
        double** array   = malloc(rows * sizeof(double* )) ;
        assert(array    != NULL) ;    
        array[0]         = malloc(rows * columns * sizeof(double)) ; 
        assert(array[0] != NULL) ;    
        for (int iC=1; iC<rows; iC++) { 
                array[iC] = array[0] + (iC * columns) ;    
        }
	for (int iC=0; iC<rows; iC++) {
                for (int jC=0; jC<columns; jC++) {
                        array[iC][jC] = initializeValue ;    
                }
	}
        return array ;
}
/*============================================================================*/
void free_double_2DArray(double** array)  
{
        free(array[0]) ;    
        array[0] = NULL ;    
        free(array) ;    
        array    = NULL ;    

        return ;
}
/*============================================================================*/
/* head of the list of fields */
/*============================================================================*/
amr_field* amr_init_fields(char* name, char* pde_type, int time_levels, int extrap_levels)
{
	amr_field* field = malloc(sizeof(amr_field)) ;

	field->name = name ;
	field->pde_type = pde_type ;
	field->time_levels = time_levels ;
	field->extrap_levels = extrap_levels ;
	field->index = 0 ;
	field->prev = NULL ;
	field->next = NULL ;
	return field ;
}
/*============================================================================*/
/* add new field to end of fields list */
/*============================================================================*/
void amr_add_field(amr_field* fields, char* name, char* pde_type, int time_levels, int extrap_levels)
{
	amr_field* current = fields ;

	while (current->next != NULL) {
		current = current->next ;
	}
	current->next = malloc(sizeof(amr_field)) ;

	current->next->name = name ;
	current->next->pde_type = pde_type ;
	current->next->time_levels = time_levels ;
	current->next->extrap_levels = extrap_levels ;
	current->next->index = (current->index + current->time_levels + current->extrap_levels) ;

	current->next->prev = current ;
	current->next->next = NULL ;

	return ;
}
/*============================================================================*/
void amr_reset_field_pde_type(amr_field *fields, char *name, char *new_pde_type)
{
	amr_field *current = fields ;
	bool found_name = false ;
	while (current->next != NULL) {
		if (strcmp(current->name,name)==0) {
			found_name = true ;
			break ;
		}
	}
	if (found_name==false) {
		assert(found_name==true) ;
	}
	current->pde_type = new_pde_type ;

	current = NULL ;
}
/*============================================================================*/
void amr_set_to_head(amr_grid** grid) 
{
	assert(grid!=NULL) ;
	while (((*grid)->parent)!=NULL) (*grid)=(*grid)->parent ;
	return ;
}
/*============================================================================*/
void amr_set_to_tail(amr_grid** grid) 
{
	assert(grid!=NULL) ;
	while ((*grid)->child!=NULL) (*grid)=(*grid)->child ;
	return ;
}
/*============================================================================*/
int amr_return_field_index(amr_field* fields, char* name) 
{
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->name,name) == 0) {
			return field->index ;
		}
	}
	return -1 ;
}
/*============================================================================*/
static int amr_delete_fields(amr_field** fields)
{
	if (*fields==NULL) {
		return -1 ;
	}
	amr_field* iter = NULL ;
	do {
		iter = (*fields)->next ;
		free(*fields) ;
		(*fields) = iter ;
	} while (iter!=NULL) ;
	
	return 0 ;
}
/*============================================================================*/
/* finds grid at specified level then sets grid pointer to that grid */
/*============================================================================*/
void amr_find_grid(int level, amr_grid_hierarchy* gh, amr_grid* grid) 
{
	assert(level<amr_max_levels-1) ;
				
	grid = gh->grids ;
	while ((grid->level)!=level) {
		grid = grid->child ;
		assert(grid!=NULL) ;
	}
	return ;
}
/*============================================================================*/
/* finds grid at finest level then sets grid pointer to that grid */
/*============================================================================*/
int amr_find_finest_grid(amr_grid_hierarchy* gh, amr_grid* grid) 
{
	grid = gh->grids ;
	while (grid->child != NULL) {
		grid = grid->child ;
	}
	printf("finest level %d\n", grid->level) ;
	return 0 ;
}
/*============================================================================*/
/* initializes new grid. all the grid functions ('grid_funcs')
   are initialized to zero */
/*============================================================================*/
amr_grid *amr_make_finer_grid(int left_coord, int right_coord, amr_grid* grid)
{
	assert((grid->level)+1<amr_max_levels) ;
	amr_grid *new_grid = malloc(sizeof(amr_grid)) ;
	assert(new_grid!=NULL) ;

	new_grid->level = (grid->level)+1 ;

	new_grid->dt = (grid->dt)/refinement ;
	new_grid->dx = (grid->dx)/refinement ;
	new_grid->time = grid->time ;
	new_grid->tC = 0 ;

	new_grid->perim_coords[0] = left_coord ;
	new_grid->perim_coords[1] = right_coord ;

	new_grid->Nx = refinement*(right_coord-left_coord) + 1 ;

	if ((grid->excised_jC)>left_coord) {
		new_grid->excised_jC = refinement*((grid->excised_jC)-left_coord) ;
	} else {
		new_grid->excised_jC = 0 ;
	}
	new_grid->bbox[0] = grid->bbox[0] + (left_coord*(grid->dx)) ;
	new_grid->bbox[1] = grid->bbox[0] + (right_coord*(grid->dx)) ;
	
	new_grid->num_grid_funcs  = grid->num_grid_funcs ;
	new_grid->grid_funcs = allocate_double_2DArray(new_grid->num_grid_funcs, new_grid->Nx, 0.) ;

	new_grid->excision_on = grid->excision_on ;

	if ((grid->perim_interior[0] == false)
	&&  (left_coord == 0)
	) {
		new_grid->perim_interior[0] = false ;
	} else {
		new_grid->perim_interior[0] = true ;
	}
	if ((grid->perim_interior[1] == false)
	&&  (right_coord == grid->Nx-1)
	) {
		new_grid->perim_interior[1] = false ;
	} else {
		new_grid->perim_interior[1] = true ;
	}
	printf("amr_make_finer_grid: made grid level %d\n", new_grid->level) ;
	printf("refinement %d\n", refinement) ;
	printf("bbox[0]\t%f\n", new_grid->bbox[0]) ;
	printf("bbox[1]\t%f\n", new_grid->bbox[1]) ;
	printf("perim_coords[0]\t%d\n", new_grid->perim_coords[0]) ;
	printf("perim_coords[1]\t%d\n", new_grid->perim_coords[1]) ;
	printf("perim_interior[0]\t%s\n", new_grid->perim_interior[0] ? "true" : "false") ;
	printf("perim_interior[1]\t%s\n", new_grid->perim_interior[1] ? "true" : "false") ;
	printf("Nx\t%d\n", new_grid->Nx) ;
	printf("dx\t%f\n", new_grid->dx) ;
	printf("dt\t%f\n", new_grid->dt) ;
	printf("num_grid_funcs\t%d\n", new_grid->num_grid_funcs) ;

	return new_grid ;
}
/*==========================================================================*/
/* linearly interpolate coarser grid to finer grid */
/*==========================================================================*/
static void interpolate_grid_func(
	int lower_coord, int Nx, double *parent_gf, double *child_gf)
{
	for (int jC=0; jC<(Nx-1); jC++) {
		double p0 = parent_gf[jC+lower_coord] ; 
		double p1 = (
			parent_gf[jC+lower_coord+1]
		-	parent_gf[jC+lower_coord]
		)/refinement ; 
		for (int kC=0; kC<refinement; kC++) {
			child_gf[2*jC+kC] = p0 + kC*p1 ;
		}
	}
	child_gf[2*(Nx-1)] = parent_gf[Nx-1] ;
}
/*==========================================================================*/
static void interpolate_all_grid_funcs(amr_grid *parent, amr_grid *child)
{
	int lower_coord = child->perim_coords[0] ;
	int upper_coord = child->perim_coords[1] ;
	int Nx = upper_coord - lower_coord + 1 ;
	for (int iC=0; iC<(parent->num_grid_funcs); iC++) {
		interpolate_grid_func(
			lower_coord, Nx, parent->grid_funcs[iC], child->grid_funcs[iC])
		;
	}
	return ;
}
/*==========================================================================*/
/* 	if the parent grid moved, then the finer grid should have different
	bounding coordinates to not move itself in physical space */
/*==========================================================================*/
static void reset_child_grid_perim_coords(
	int old_lower_coord, amr_grid *base_grid)
{
	for (amr_grid *grid=(base_grid->child); grid!=NULL; grid=(grid->child)) {
		int current_lower_coord = grid->parent->perim_coords[0] ;
		int grid_lower_coord = grid->perim_coords[0] ;

		grid->perim_coords[0] += refinement*(old_lower_coord-current_lower_coord) ;
		grid->perim_coords[1] += refinement*(old_lower_coord-current_lower_coord) ;

		old_lower_coord = grid_lower_coord ;
	}
	return ;
}
/*==========================================================================*/
/* 	already injected old child grid, so no need to copy again */
/*==========================================================================*/
static void add_flagged_child_grid(amr_grid *grid)
{
	int new_lower_coord = grid->flagged_coords[0] ;
	int new_upper_coord = grid->flagged_coords[1] ;

	amr_grid *old_child = grid->child ;
	int old_lower_coord = 0 ;
	if (old_child!=NULL) {
		old_lower_coord = old_child->perim_coords[0] ;
	}
	if ((new_upper_coord-new_lower_coord)<min_grid_size) {
		if ((old_child!=NULL)
		&& ((old_child->child)==NULL)
		) {
			amr_destroy_grid(old_child) ;
			old_child = NULL ;
		}
		return ;
	}
	amr_grid *new_child = amr_make_finer_grid(
		new_lower_coord, new_upper_coord, grid
	) ;
	interpolate_all_grid_funcs(grid, new_child) ;

	if (old_child!=NULL) {
		amr_destroy_grid(old_child) ; 
		old_child = NULL ;	
	} 
	amr_insert_grid(new_child, grid) ;
	reset_child_grid_perim_coords(old_lower_coord, grid->child) ;

	return ;
}
/*============================================================================*/
void amr_add_finer_grid(int left_coord, int right_coord, amr_grid* parent)
{
	amr_grid *new_grid = amr_make_finer_grid(left_coord, right_coord, parent) ;

	parent->child = new_grid ;
	new_grid->parent = parent ;
	new_grid->child = NULL ;

	return ;
}
/*============================================================================*/
/* set the base/level 0 (shadow) grid and the level one grid */
/*============================================================================*/
amr_grid_hierarchy* amr_init_grid_hierarchy(
	amr_field* fields,
	int Nt, int Nx, int t_step_save,
	double cfl_num,
	double bbox[2],
	bool excision_on)
{
	amr_grid_hierarchy* gh = malloc(sizeof(amr_grid_hierarchy)) ;

	gh->cfl_num = cfl_num ;
	gh->Nt  = Nt ;
	gh->t_step_save = t_step_save ;
	gh->fields = fields ;
	gh->excision_on = excision_on ;

	int num_grid_funcs = 0 ;
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		num_grid_funcs += field->time_levels + field->extrap_levels ;
	} 
/*---------------------------	
 * base (shadow) grid 
 *-------------------------*/
	amr_grid* base_grid = malloc(sizeof(amr_grid)) ;
	assert(base_grid != NULL) ;	
	base_grid->level = 0 ;

	base_grid->grid_funcs = allocate_double_2DArray(num_grid_funcs, Nx, 0.) ; 	
	base_grid->num_grid_funcs  = num_grid_funcs ;
	
	base_grid->Nx = Nx ;
	base_grid->excised_jC = 0 ;
	base_grid->tC = 0 ;

	base_grid->bbox[0] = bbox[0] ;
	base_grid->bbox[1] = bbox[1] ;

	base_grid->dx = (bbox[1]-bbox[0])/(Nx-1.) ;
	base_grid->dt = cfl_num*base_grid->dx ;
	base_grid->time = 0 ;

	base_grid->perim_interior[0] = false ;
	base_grid->perim_interior[1] = false ;
	base_grid->perim_coords[0] = 0 ;
	base_grid->perim_coords[1] = Nx-1 ;

	base_grid->child = NULL ;
	base_grid->parent = NULL ;

	base_grid->excision_on = excision_on ;
	
	gh->grids = base_grid ;

	printf("bbox[0]\t%f\n", gh->grids->bbox[0]) ;
	printf("bbox[1]\t%f\n", gh->grids->bbox[1]) ;
	printf("perim_coords[0]\t%d\n", gh->grids->perim_coords[0]) ;
	printf("perim_coords[1]\t%d\n", gh->grids->perim_coords[1]) ;
	printf("Nx\t%d\n", gh->grids->Nx) ;
	printf("dx\t%f\n", gh->grids->dx) ;
	printf("dt\t%f\n", gh->grids->dt) ;
	printf("num_grid_funcs\t%d\n", gh->grids->num_grid_funcs) ;
	fflush(NULL) ;
/*	level one grid */
	amr_add_finer_grid(0, Nx-1, base_grid) ;
	
	return gh ;
}
/*==========================================================================*/
/* restriction along shared grid points */
/*==========================================================================*/
static void inject_grid_func(
	int Nx, int perim_coord_left, double *gf_parent, double *gf)
{
	for (int iC=0; iC<Nx; iC++) {
		if (iC%refinement==0) {
			gf_parent[perim_coord_left+(iC/refinement)] = gf[iC] ;
		}
	}
	return ;
}
/*==========================================================================*/
/* injecting to parent grid */
/*==========================================================================*/
void inject_overlaping_fields(
		amr_field *fields, amr_grid *grid, amr_grid *parent)
{
	int index = 0 ;
	for (amr_field *field=fields; field!=NULL; field=field->next) {
		index = field->index ;
		inject_grid_func(
			grid->Nx, 
			grid->perim_coords[0], 
			parent->grid_funcs[index],
			grid->grid_funcs[index]
		) ;
	}
	return ;
}
/*==========================================================================*/
/* flagging all hyperbolic and ode fields */
/*==========================================================================*/
void flag_field_regridding_coords(amr_field *fields, amr_grid *parent, amr_grid *grid)
{
	double regrid_err_lim = trunc_err_tolerance*pow(grid->dt,2) ;
	for (amr_field *field=fields; field!=NULL; field=(field->next)) {
		if ((strcmp(field->pde_type,HYPERBOLIC)==0)
//		||  (strcmp(field->pde_type,ODE)==0) 
		) {
			/* keep on going */
		} else {
			field->flagged_coords[0] = (grid->Nx)-1 ;
			field->flagged_coords[1] = 0 ;
			continue ;
		}
		int field_index = field->index ;
		int lower_flagged_coords = (-1) ;
		int upper_flagged_coords = (-1) ;
		int lower_jC = grid->perim_coords[0] ;
		int upper_jC = grid->perim_coords[1] ;
		int start_jC = (grid->perim_interior[0]==true) 
		?	lower_jC + buffer_coord
		:	lower_jC 
		;
		int end_jC = (grid->perim_interior[1]==true) 
		?	upper_jC - buffer_coord
		:	upper_jC 
		;
		for (int jC=start_jC; jC<end_jC; jC++) {
			double parent_val = parent->grid_funcs[field_index][jC] ;

			int grid_index = refinement*(jC-lower_jC) ; 
			double grid_val = grid->grid_funcs[field_index][grid_index] ;

			double trunc_err = fabs(parent_val-grid_val) ;

			if (trunc_err > regrid_err_lim) {
				if (lower_flagged_coords==(-1)) {
					lower_flagged_coords = grid_index ;
					upper_flagged_coords = grid_index ;
				} else {
					upper_flagged_coords = grid_index ;
				} 
			}
		}
		field->flagged_coords[0] = lower_flagged_coords ;
		field->flagged_coords[1] = upper_flagged_coords ;
//		printf("field %s\t", field->name) ;
//		printf("lower %d\tupper %d\n", field->flagged_coords[0], field->flagged_coords[1]) ;
	}
	return ;
}
/*==========================================================================*/
static inline double max_double(double val_1, double val_2) 
{
	return (val_1>val_2) ? val_1 : val_2 ;
}
static inline double min_double(double val_1, double val_2) 
{
	return (val_1<val_2) ? val_1 : val_2 ;
}
/*==========================================================================*/
static void determine_grid_coords(
	amr_field *fields, amr_grid *grid)
{
	int Nx = grid->Nx ;
	int lower_child_grid_coord = Nx-1 ; 
	int upper_child_grid_coord = 0 ; 
/* find min and max cords */
	for (amr_field *field=fields; field!=NULL; field=(field->next)) {
		int lower_coord = field->flagged_coords[0] ;
		int upper_coord = field->flagged_coords[1] ;
		if ((	(strcmp(field->pde_type,HYPERBOLIC)==0) 
//		     || (strcmp(field->pde_type,ODE)==0)
		) &&  	(lower_coord!=(-1) && upper_coord!=(-1)) 
		) {
			lower_child_grid_coord = min_double(lower_coord,lower_child_grid_coord) ;
			upper_child_grid_coord = max_double(upper_coord,upper_child_grid_coord) ;
		}
	}
/* if too close to physical boundary place adjacent to boundary */
	if ((grid->perim_interior[0]==false)
	&&  (lower_child_grid_coord<buffer_coord)
	) {
		lower_child_grid_coord = 0 ;
	}
	if ((grid->perim_interior[1]==false)
	&&  (upper_child_grid_coord>(Nx-1-buffer_coord))
	) {
		upper_child_grid_coord = Nx-1 ;
	}

	grid->flagged_coords[0] = lower_child_grid_coord ;
	grid->flagged_coords[1] = upper_child_grid_coord ;

	return ;
}
/*==========================================================================*/
void regrid_all_finer_levels(amr_field *fields, amr_grid *base_grid)
{
	amr_grid *grid = base_grid ;
	amr_set_to_tail(&grid) ;
	while ((grid->level)>=(base_grid->level)) {
		if ((grid->child)!=NULL) {
			inject_overlaping_fields(fields, grid->child, grid) ;
		}
		flag_field_regridding_coords(fields, grid->parent, grid) ;
		determine_grid_coords(fields, grid) ;
		printf("level %d\n", grid->level) ;
		printf("lower %d\tupper %d\n", grid->flagged_coords[0], grid->flagged_coords[1]) ;
		fflush(NULL) ;
		if ((grid->level)<amr_max_levels-1) {
			add_flagged_child_grid(grid) ;
		}
		grid = grid->parent ;
	} 
	grid = NULL ;
	return ;
}
/*============================================================================*/
/* starting from level 1 grid, add self similar grid hierarchy from origin */
/*============================================================================*/
void add_self_similar_initial_grids(
	amr_grid_hierarchy* gh, int grid_size_ratio, int num_grids) 
{
	amr_grid* grid = gh->grids->child ; /* start at grid level=1 */
	int Nx = grid->Nx ; 

	for (int iC=0; iC<num_grids; iC++) {
		amr_add_finer_grid(0, (int)(Nx/grid_size_ratio), grid) ;
		grid = grid->child ;
		Nx = grid->Nx ;
	}
	grid = grid->child ;
	grid = NULL ;
	return ; 
}
/*============================================================================*/
void amr_insert_grid(amr_grid *grid_to_insert, amr_grid* grid)
{
	grid_to_insert->child  = grid->child ;
	grid_to_insert->parent = grid ;
	if ((grid->child)!=NULL) {
		(grid->child)->parent = grid_to_insert ;
	}
	grid->child = grid_to_insert ;

	return ;
} 
/*============================================================================*/
void amr_destroy_grid(amr_grid* grid) 
{	
	free_double_2DArray(grid->grid_funcs) ;

	if ((grid->parent)!=NULL) {
		(grid->parent)->child = grid->child  ;
	} 
	if ((grid->child)!=NULL) {
		(grid->child)->parent = grid->parent ;
	} 
	free(grid) ;
	grid = NULL ;

	return ;
}
/*============================================================================*/
void amr_destroy_grid_hierarchy(amr_grid_hierarchy* gh) 
{
	amr_grid* grid = gh->grids ;
	amr_grid* child = NULL ;

	do {
		child = grid->child ;
		amr_destroy_grid(grid) ;
		grid = child ;
	} while (grid != NULL) ;
	
	amr_delete_fields(&(gh->fields)) ;
	free(gh) ;
	gh = NULL ;
	
	return ;
}
