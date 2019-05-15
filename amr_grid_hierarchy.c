/*===========================================================================*/
/*===========================================================================*/
#include "amr_grid_hierarchy.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

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
int amr_get_grid_funcs(struct amr_grid* grid, double** grid_funcs) 
{
	for (int iC=0; iC<grid->num_grid_funcs; iC++) {
		grid_funcs[iC] = grid->grid_funcs[iC] ;
	}
	return 0 ;
}
/*============================================================================*/
/* finds grid at specified level then sets grid pointer to that grid */
/*============================================================================*/
int amr_find_grid(int level, struct amr_grid_hierarchy* gh, struct amr_grid* grid) 
{
	if (level>AMR_MAX_LEVELS-1) {
		printf("ERROR(find_grid): level>AMR_MAX_LEVELS-1\n") ;
		return -1 ;
	} 
	grid = gh->grid ;
	while (grid->level != level) {
		grid = grid->child ;
		if (grid == NULL) {
			printf("ERROR(find_grid): grid==NULL\n") ;
			return -1 ;
		}
	}
	return 0 ;
}
/*============================================================================*/
/* finds grid at finest level then sets grid pointer to that grid */
/*============================================================================*/
int amr_find_finest_grid(struct amr_grid_hierarchy* gh, struct amr_grid* grid) 
{
	grid = gh->grid ;
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
int amr_add_finer_grid(int left_coord, int right_coord, struct amr_grid* parent)
{
	if (parent->level+1 >= AMR_MAX_LEVELS) {
		printf("ERROR(amr_add_grid): level>=AMR_MAX_LEVELS\n") ;
		return -1 ;
	}
	struct amr_grid* new_grid = malloc(sizeof(struct amr_grid)) ;
	if (new_grid == NULL) {
		printf("ERROR(amr_add_grid): new_grid=NULL\n") ;
		return -1 ;
	}

	parent->child = new_grid ;
	new_grid->parent = parent ;
	new_grid->child = NULL ;
	new_grid->sibling = NULL ;
	new_grid->neighbor = NULL ;

	new_grid->level = parent->level+1 ;

	new_grid->dt = parent->dt/REFINEMENT ;
	new_grid->dx = parent->dx/REFINEMENT ;

	new_grid->perim_coords[0] = left_coord ;
	new_grid->perim_coords[1] = right_coord ;

	new_grid->Nx = REFINEMENT*(right_coord-left_coord) + 1 ;

	new_grid->bbox[0] = parent->bbox[0] + (left_coord *parent->dx) ;
	new_grid->bbox[1] = parent->bbox[0] + (right_coord*parent->dx) ;
	
	new_grid->num_grid_funcs = parent->num_grid_funcs ;
	new_grid->grid_funcs = allocate_double_2DArray(new_grid->num_grid_funcs, new_grid->Nx, 0.) ;

	if ((parent->perim_interior[0] == false)
	&&  (left_coord == 0)
	) {
		new_grid->perim_interior[0] = false ;
	}
	if ((parent->perim_interior[1] == false)
	&&  (right_coord == parent->Nx)
	) {
		new_grid->perim_interior[1] = false ;
	}
	printf("amr_add_finer_grid: made grid level %d\n", new_grid->level) ;
	printf("bbox[0]\t%f\n", new_grid->bbox[0]) ;
	printf("bbox[1]\t%f\n", new_grid->bbox[1]) ;
	printf("perim_coords[0]\t%d\n", new_grid->perim_coords[0]) ;
	printf("perim_coords[1]\t%d\n", new_grid->perim_coords[1]) ;
	printf("Nx\t%d\n", new_grid->Nx) ;
	printf("dx\t%f\n", new_grid->dx) ;
	printf("dt\t%f\n", new_grid->dt) ;
	printf("num_grid_funcs\t%d\n", new_grid->num_grid_funcs) ;

	return 0 ;
}
/*============================================================================*/
/* set the base/level 0 (shadow) grid and the level one grid */
/*============================================================================*/
struct amr_grid_hierarchy* amr_init_grid_hierarchy(
	int num_grid_funcs,
	int Nt, int Nx, int t_step_save,
	double cfl_num,
	double bbox[2],
	bool excision_on)
{
	struct amr_grid_hierarchy* gh = malloc(sizeof(struct amr_grid_hierarchy)) ;

	gh->cfl_num = cfl_num ;
	gh->Nt  = Nt ;
	gh->t_step_save = t_step_save ;

/*	base (shadow) grid */
	struct amr_grid* base_grid = malloc(sizeof(struct amr_grid)) ;
	assert(base_grid != NULL) ;	

	base_grid->level = 0 ;
	
	base_grid->grid_funcs = allocate_double_2DArray(num_grid_funcs, Nx, 1.) ; 
	
	base_grid->num_grid_funcs = num_grid_funcs ;
	base_grid->Nx = Nx ;

	base_grid->bbox[0] = bbox[0] ;
	base_grid->bbox[1] = bbox[1] ;

	base_grid->dx = (bbox[1]-bbox[0])/(Nx-1.) ;
	base_grid->dt = cfl_num*base_grid->dx ;

	base_grid->perim_interior[0] = false ;
	base_grid->perim_interior[1] = false ;
	base_grid->perim_coords[0] = 0 ;
	base_grid->perim_coords[1] = Nx-1 ;

	base_grid->child = NULL ;
	base_grid->parent = NULL ;
	base_grid->sibling = NULL ;
	base_grid->neighbor = NULL ;	
	
	gh->grid = base_grid ;

	printf("bbox[0]\t%f\n", gh->grid->bbox[0]) ;
	printf("bbox[1]\t%f\n", gh->grid->bbox[1]) ;
	printf("perim_coords[0]\t%d\n", gh->grid->perim_coords[0]) ;
	printf("perim_coords[1]\t%d\n", gh->grid->perim_coords[1]) ;
	printf("Nx\t%d\n", gh->grid->Nx) ;
	printf("dx\t%f\n", gh->grid->dx) ;
	printf("dt\t%f\n", gh->grid->dt) ;
	printf("num_grid_funcs\t%d\n", gh->grid->num_grid_funcs) ;
	fflush(NULL) ;
/*	level one grid */
	amr_add_finer_grid(0, Nx-1, base_grid) ;
	
	return gh ;
}
/*============================================================================*/
int amr_destroy_grid(struct amr_grid* grid) 
{	
	free_double_2DArray(grid->grid_funcs) ;

	if (grid->parent != NULL) {
		grid->parent->child = grid->child  ;
	} 
	if (grid->child != NULL) {
		grid->child->parent = grid->parent ;
	} 
	free(grid) ;
	grid = NULL ;

	return 0 ;
}
/*============================================================================*/
int amr_destroy_grid_hierarchy(struct amr_grid_hierarchy* gh) 
{
	struct amr_grid* grid = gh->grid ;
	struct amr_grid* child = NULL ;

	do {
		child = grid->child ;
		amr_destroy_grid(grid) ;
		grid = child ;
	} while (grid != NULL) ;

	free(gh) ;
	gh = NULL ;
	
	return 0 ;
}
