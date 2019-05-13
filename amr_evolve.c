#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "amr_grid_hierarchy.h"
#include "amr_evolve.h"

/*==========================================================================*/
/* for now its linear prolongation */
/*==========================================================================*/
/*static*/ void prolong(int Nx_grid, int perim_coord_left, double* parent, double* grid)
{
	double coef_0 = 0 ;
	double coef_1 = 0 ;

	for (int iC=0; iC<Nx_grid-REFINEMENT; iC++) {
		if (iC%REFINEMENT==0) {
			coef_0 = parent[perim_coord_left+(iC/REFINEMENT)] ;
			coef_1 = (
				parent[perim_coord_left+(iC/REFINEMENT)+1] 
			-	parent[perim_coord_left+(iC/REFINEMENT)+0] 
			)/REFINEMENT
			;
			for (int jC=0; jC<REFINEMENT; jC++) {
				grid[iC+jC] = coef_0 + (coef_1*jC) ;
			}
		}
	}
}
/*==========================================================================*/
/* restriction */
/*==========================================================================*/
/*static*/ void inject(int Nx_grid, int perim_coord_left, double* parent, double* grid)
{
	for (int iC=0; iC<Nx_grid; iC++) {
		if (iC%REFINEMENT==0) {
			parent[perim_coord_left+(iC/REFINEMENT)] = grid[iC] ;
		}
	}
}
/*==========================================================================*/
void amr_evolve( struct grid* grid,
	int num_t_steps,
	void (*evolve_pde)(void))
{
	for (int tC=0; tC<num_t_steps; tC++) {
		grid->tC += 1 ;
		if (grid->tC%REGRID == 0) {
			/* regrid all finer levels */
		}
		if (grid->perim_interior[0] == true) { 
			/*interpolate inner*/
		}
		if (grid->perim_interior[1] == true) { 
			/*interpolate outer*/
		}
		evolve_pde() ;
		if (grid->child != NULL) {
			amr_evolve(
				grid->child,
				REFINEMENT,
				evolve_pde)
			;
		}
	}
	if (grid->parent != NULL) {
		/* compute truncation error */
		/* inject */
	}
	return ;
}
