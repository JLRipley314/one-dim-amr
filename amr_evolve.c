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
static void inject_grid_func(
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
static void inject_all_grid_funcs(struct amr_grid* parent, struct amr_grid* grid)
{
	for (int iC=0; iC<(grid->num_grid_funcs); iC++) {
		inject_grid_func(
			grid->Nx, 
			grid->perim_coords[0], 
			parent->grid_funcs[iC],
			grid->grid_funcs[iC]
		) ;
	}
	return ;
}
/*==========================================================================*/
/* {ordered field_n, field_nm1, ...} */
/*==========================================================================*/
static void set_interior_hyperbolic_boundary(
	struct amr_grid* parent,
	struct amr_grid* grid)
{
	double coef_0, coef_1 ; 

	int perim ;

	int num_time_levels = 2 ;
	int num_grid_funcs = grid->num_grid_funcs ;

	for (int iC=0; iC<num_grid_funcs; iC += num_time_levels) {

		if (grid->perim_interior[0] == true) {
			perim  =  grid->perim_coords[0] ;
			coef_0 =  parent->grid_funcs[iC+1][perim] ;
			coef_1 = (
				parent->grid_funcs[iC  ][perim]
			-	parent->grid_funcs[iC+1][perim]
			)/REFINEMENT ;

			grid->grid_funcs[iC  ][0] = coef_0 + coef_1*(((grid->tC)%REFINEMENT)+1) ;
			grid->grid_funcs[iC+1][0] = coef_0 + coef_1*(((grid->tC)%REFINEMENT)  ) ;
		}
		if (grid->perim_interior[1] == true) {
			perim  =  grid->perim_coords[1] ;
			coef_0 =  parent->grid_funcs[iC+1][perim] ;
			coef_1 = (
				parent->grid_funcs[iC  ][perim]
			-	parent->grid_funcs[iC+1][perim]
			)/REFINEMENT ;

			grid->grid_funcs[iC  ][grid->Nx-1] = coef_0 + coef_1*(((grid->tC)%REFINEMENT)+1) ;
			grid->grid_funcs[iC+1][grid->Nx-1] = coef_0 + coef_1*(((grid->tC)%REFINEMENT)  ) ;
		}
	}

	return ;
}
/*==========================================================================*/
/* evolves all grids in hierarchy */
/*==========================================================================*/
static void amr_evolve_grid(
	struct amr_grid* grid,
	int num_t_steps,
	void (*evolve_pde)(struct amr_grid*))
{
	for (int tC=0; tC<num_t_steps; tC++) {
		grid->tC   += 1 ;
		grid->time += grid->dt ;
		if (grid->tC%REGRID == 0) {
			/* regrid all finer levels */
		}
		if (grid->parent != NULL) {
			set_interior_hyperbolic_boundary(grid->parent, grid) ;
		}
		evolve_pde(grid) ;
		if (grid->child != NULL) {
			amr_evolve_grid(
				grid->child,
				REFINEMENT,
				evolve_pde)
			;
		}
	}
	if (grid->parent != NULL) {
		/* TO DO: compute truncation error */
		inject_all_grid_funcs(grid->parent, grid) ;
	}
	return ;
}
/*==========================================================================*/
/* sets initial data on all pre assigned grid levels */
/*==========================================================================*/
static void set_initial_data(
	struct amr_grid_hierarchy* gh,
	void (*initial_data)(struct amr_grid*))
{
	struct amr_grid* grid = gh->grid ;
	while (grid != NULL) {
		initial_data(grid) ;
		grid = grid->child ;
	}
	return ;
}
/*==========================================================================*/
/* evolves all grids in hierarchy */
/*==========================================================================*/
void amr_main(
	struct amr_grid_hierarchy* gh, 
	void (*initial_data)(struct amr_grid*),
	void (*evolve_pde)(struct amr_grid*),
	void (*save_to_file)(struct amr_grid*))
{
	set_initial_data(gh,initial_data) ;
	save_to_file(gh->grid->child) ;

	for (int tC=1; tC<(gh->Nt); tC++) {
		amr_evolve_grid(
			gh->grid,
			1,
			evolve_pde) 
		;
		if (tC%(gh->t_step_save)==0) {
			for (struct amr_grid* grid = gh->grid; grid != NULL; grid=grid->child) {
				save_to_file(grid) ;
			}
		}
	}
	return ;	
}
