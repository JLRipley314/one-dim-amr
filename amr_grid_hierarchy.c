/*===========================================================================*/
/*===========================================================================*/
#include "amr_grid_hierarchy.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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
int amr_add_grid(struct amr_grid* grid)
{
	if (grid->level+1 >= AMR_MAX_LEVELS) {
		printf("ERROR(amr_add_grid): level>=AMR_MAX_LEVELS\n") ;
		return -1 ;
	}
	struct amr_grid* new_grid = malloc(sizeof(struct amr_grid)) ;

	new_grid->parent = grid ;
	new_grid->child = NULL ;
	new_grid->sibling = NULL ;

	new_grid->level = grid->level+1 ;

	new_grid->dt = grid->dt/REFINEMENT ;
	new_grid->dx = grid->dx/REFINEMENT ;

	new_grid->bbox[0] = grid->bbox[0] + (new_grid->perim_coord[0]*grid->dx) ;
	new_grid->bbox[1] = grid->bbox[0] + (new_grid->perim_coord[1]*grid->dx) ;

	if ((grid->perim_interior[0] == false)
	&&  (new_grid->perim_coord[0] == 0)
	) {
		new_grid->perim_interior[0] = false ;
	}
	if ((grid->perim_interior[1] == false)
	&&  (new_grid->perim_coord[1] == grid->Nx)
	) {
		new_grid->perim_interior[1] = false ;
	}
	return 0 ;
}
