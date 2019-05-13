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
/*static*/ void prolong(
	int Nx, int perim_coord_left, double* gf, double* gf_child)
{
	double coef_0 = 0 ;
	double coef_1 = 0 ;

	for (int iC=0; iC<Nx-REFINEMENT; iC++) {
		if (iC%REFINEMENT==0) {
			coef_0 = gf[perim_coord_left+(iC/REFINEMENT)] ;
			coef_1 = (
				gf[perim_coord_left+(iC/REFINEMENT)+1] 
			-	gf[perim_coord_left+(iC/REFINEMENT)+0] 
			)/REFINEMENT
			;
			for (int jC=0; jC<REFINEMENT; jC++) {
				gf_child[iC+jC] = coef_0 + (coef_1*jC) ;
			}
		}
	}
	return ;
}
/*==========================================================================*/
/* restriction along shared grid points */
/*==========================================================================*/
/*static*/ void inject_grid_func(
	int Nx, int perim_coord_left, double* gf_parent, double* gf)
{
	for (int iC=0; iC<Nx; iC++) {
		if (iC%REFINEMENT==0) {
			gf_parent[perim_coord_left+(iC/REFINEMENT)] = gf[iC] ;
		}
	}
	return ;
}
/*==========================================================================*/
/*static*/ void inject_all_grid_funcs()
{
	return ;
}
/*==========================================================================*/
/* evolves all grids in hierarchy */
/*==========================================================================*/
static void amr_evolve_grid(
	struct amr_grid* grid,
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
			amr_evolve_grid(
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
/*==========================================================================*/
/* evolves all grids in hierarchy */
/*==========================================================================*/
void amr_main(
	struct amr_grid_hierarchy* gh, 
	int num_t_steps,
	int save_time,
	void (*evolve_pde)(void),
	void (*save_to_file)(void))
{
	for (int tC=0; tC<num_t_steps; tC++) {
		amr_evolve_grid(
			gh->grid,
			1,
			evolve_pde) 
		;
		if (tC%save_time) {
			save_to_file() ;
		}
	}
}
