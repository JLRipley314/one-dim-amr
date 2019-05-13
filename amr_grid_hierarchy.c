/*===========================================================================*/
/*===========================================================================*/
#include "amr_grid_hierarchy.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*============================================================================*/
void amr_get_grid_funcs(struct grid* grid, double** grid_funcs) 
{
	for (int iC=0; iC<grid->num_grid_funcs; iC++) {
		grid_funcs[iC] = grid->grid_funcs[iC] ;
	}
}

