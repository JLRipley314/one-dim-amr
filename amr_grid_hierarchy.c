/*===========================================================================*/
/*===========================================================================*/
#include "amr_grid_hierarchy.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*============================================================================*/
static double** allocate_double_2DArray(int rows, int columns, double initializeValue)  
{
        double **pa_array   = malloc(rows * sizeof(double *)) ;
        assert(pa_array    != NULL) ;    
        pa_array[0]         = malloc(rows * columns * sizeof(double)) ; 
        assert(pa_array[0] != NULL) ;    
        for (int iC=1; iC<rows; iC++) { 
                pa_array[iC] = pa_array[0] + (iC * columns) ;    
                for (int jC=0; jC<columns; jC++) {
                        pa_array[iC][jC] = initializeValue ;    
                }
        }
        return pa_array ;
}
/*============================================================================*/
static void free_double_2DArray(double** array)  
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
int amr_add_grid(int perim_coords[2], struct amr_grid* grid)
{
	if (grid->level+1 >= AMR_MAX_LEVELS) {
		printf("ERROR(amr_add_grid): level>=AMR_MAX_LEVELS\n") ;
		return -1 ;
	}
	struct amr_grid* new_grid = malloc(sizeof(struct amr_grid)) ;
	if (new_grid == NULL) {
		printf("ERROR(amr_add_grid): new_grid=NULL\n") ;
		return -1 ;
	}

	new_grid->parent = grid ;
	new_grid->child = NULL ;
	new_grid->sibling = NULL ;

	new_grid->level = grid->level+1 ;

	new_grid->dt = grid->dt/REFINEMENT ;
	new_grid->dx = grid->dx/REFINEMENT ;

	new_grid->perim_coords[0] = perim_coords[0] ;
	new_grid->perim_coords[1] = perim_coords[1] ;

	new_grid->Nx = REFINEMENT*(perim_coords[1]-perim_coords[0]) + 1 ;

	new_grid->bbox[0] = grid->bbox[0] + (perim_coords[0]*grid->dx) ;
	new_grid->bbox[1] = grid->bbox[0] + (perim_coords[1]*grid->dx) ;
	
	new_grid->num_grid_funcs = grid->num_grid_funcs ;
	new_grid->grid_funcs = allocate_double_2DArray(new_grid->num_grid_funcs, new_grid->Nx, 0.) ; 


	if ((grid->perim_interior[0] == false)
	&&  (perim_coords[0] == 0)
	) {
		new_grid->perim_interior[0] = false ;
	}
	if ((grid->perim_interior[1] == false)
	&&  (perim_coords[1] == grid->Nx)
	) {
		new_grid->perim_interior[1] = false ;
	}
	return 0 ;
}
/*============================================================================*/
int amr_destroy_grid(struct amr_grid* grid) 
{
	
	free_double_2DArray(grid->grid_funcs) ;

	grid->parent->child = grid->child  ;
	grid->child->parent = grid->parent ;

	free(grid) ;

	return 0 ;
}
