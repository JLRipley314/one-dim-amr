#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "amr_grid_hierarchy.h"
#include "amr_evolve.h"

/*==========================================================================*/
/* routines for copying time steps to next level */
/*==========================================================================*/
static void copy_to_2nd_array(int Nx, double* field_1, double* field_2) 
{
	for (int iC=0; iC<Nx; iC++) {
		field_2[iC] = field_1[iC] ;
	}
	return ;
}
/*==========================================================================*/
/* shifting fields from farthest back in time to most recent in time */
/*==========================================================================*/
static void shift_field(int index, int time_levels, amr_grid* grid)
{
	int Nx = grid->Nx ;
	for (int iC=index+time_levels-1; iC>index; iC--) {
		copy_to_2nd_array(Nx, grid->grid_funcs[iC-1], grid->grid_funcs[iC]) ;
	}
}
/*==========================================================================*/
static void shift_fields_one_time_level(
	amr_field* fields,
	amr_grid* grid)
{
	int index, time_levels ;
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		index = field->index ;
		time_levels = field->time_levels ;
		shift_field(index,time_levels,grid) ;
	}
	return ;
}
/*==========================================================================*/
static void shift_grids_one_time_level(amr_grid_hierarchy* gh)
{
	for (amr_grid* grid=gh->grids; grid!=NULL; grid=grid->child) {
		shift_fields_one_time_level(gh->fields, grid) ;
	}
}
/*==========================================================================*/
/* for now its linear prolongation */
/*==========================================================================*/
/*static*/ void prolong_grid_func(
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
/* compute truncation error for given grid function*/
/*==========================================================================*/
/*static*/ void compute_truncation_error_grid_func(
	int Nx, int perim_coord_left, double* gf_parent, double* gf, int* trunc_err_flags)
{
	trunc_err_flags[0] = -1 ;
	trunc_err_flags[1] = -1 ;
	double trunc_err = 0 ;

	for (int iC=0; iC<Nx; iC++) {
		if (iC%REFINEMENT==0) {
			trunc_err = fabs(gf_parent[perim_coord_left+(iC/REFINEMENT)] - gf[iC]) ;
			if ((trunc_err > TRUNC_ERR_TOLERANCE)
			&&  (trunc_err_flags[0] == -1)
			) {
				trunc_err_flags[0] = iC ;
			}
			if ((trunc_err > TRUNC_ERR_TOLERANCE)
			&&  (trunc_err_flags[0] > -1)
			) {
				trunc_err_flags[1] = iC ;
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
static void inject_overlaping_fields(
	amr_field* fields, amr_grid* parent, amr_grid* grid)
{
	int index = 0 ;
	for (amr_field* iter=fields; iter!=NULL; iter=iter->next) {
		index = iter->index ;
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
/* {ordered field_n, field_nm1, ...} */
/*==========================================================================*/
/*static*/ void set_interior_hyperbolic_boundary_linear_interp(
	amr_field* field,
	amr_grid* parent,
	amr_grid* grid)
{
	double coef_0, coef_1 ; 

	int perim ;

	int tC = (grid->tC)%REFINEMENT ;

	int index = field->index ;

	if (grid->perim_interior[0] == true) {
		perim  =  grid->perim_coords[0] ;
		coef_0 =  parent->grid_funcs[index][perim] ;
		coef_1 = (
			parent->grid_funcs[index  ][perim]
		-	parent->grid_funcs[index+1][perim]
		)/REFINEMENT ;

		grid->grid_funcs[index  ][0] = coef_0 + coef_1*(tC+1) ;
		grid->grid_funcs[index+1][0] = coef_0 + coef_1*(tC  ) ;
	}
	if (grid->perim_interior[1] == true) {
		perim  =  grid->perim_coords[1] ;
		coef_0 =  parent->grid_funcs[index+1][perim] ;
		coef_1 = (
			parent->grid_funcs[index  ][perim]
		-	parent->grid_funcs[index+1][perim]
		)/REFINEMENT ;

		grid->grid_funcs[index  ][grid->Nx-1] = coef_0 + coef_1*(tC+1) ;
		grid->grid_funcs[index+1][grid->Nx-1] = coef_0 + coef_1*(tC  ) ;
	}

	return ;
}
/*==========================================================================*/
/* We interpolate by computing the Taylor expansion about grid point nm1
   to second order (quadratic order) 
   We need at least three time levels for this to work */
/*==========================================================================*/
/*static*/ void set_interior_hyperbolic_boundary_quad_interp(
	amr_field* field, 
	amr_grid* parent,
	amr_grid* grid)
{
	double coef_0, coef_1, coef_2 ; 

	int perim ;

	int tC = (grid->tC)%REFINEMENT ;

	int index = field->index ;

	if (grid->perim_interior[0] == true) {
		perim  =  grid->perim_coords[0] ;
		coef_0 =  parent->grid_funcs[index][perim] ;
		coef_1 = (
			parent->grid_funcs[index  ][perim]
		- 	parent->grid_funcs[index+2][perim]
		)/(2*REFINEMENT) ;
		coef_2 = (
			   parent->grid_funcs[index  ][perim]
		- 	(2*parent->grid_funcs[index+1][perim])
		+    	   parent->grid_funcs[index+2][perim]
		)/pow(REFINEMENT,2) ;

		grid->grid_funcs[index  ][0] = coef_0 + coef_1*(tC+1) + 0.5*coef_2*pow(tC,2) ;
		grid->grid_funcs[index+1][0] = coef_0 + coef_1*(tC  ) + 0.5*coef_2*pow(tC,2) ;
	}
	if (grid->perim_interior[1] == true) {
		perim  =  grid->perim_coords[1] ;
		coef_0 =  parent->grid_funcs[index][perim] ;
		coef_1 = (
			parent->grid_funcs[index  ][perim]
		- 	parent->grid_funcs[index+2][perim]
		)/(2*REFINEMENT) ;
		coef_2 = (
			   parent->grid_funcs[index  ][perim]
		-	(2*parent->grid_funcs[index+1][perim])
		+	   parent->grid_funcs[index+2][perim]
		)/pow(REFINEMENT,2) ;

		grid->grid_funcs[index  ][grid->Nx-1] = coef_0 + coef_1*(tC+1) + 0.5*coef_2*pow(tC,2) ;
		grid->grid_funcs[index+1][grid->Nx-1] = coef_0 + coef_1*(tC  ) + 0.5*coef_2*pow(tC,2) ;
	}
	return ;
}
/*==========================================================================*/
/* using quadratic interpolation in time for boundaries */
/*==========================================================================*/
void amr_set_interior_hyperbolic_boundary(
	amr_field* fields, amr_grid* parent, amr_grid* grid)
{
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->pde_type,HYPERBOLIC) == 0) {
			set_interior_hyperbolic_boundary_quad_interp(
				field, parent, grid)
			;
		}
	}
	return ;
}
/*==========================================================================*/
/* ode methods */
/*==========================================================================*/
void amr_linear_extrapapolate_field(amr_field* field, amr_grid* grid) 
{
	int field_index = field->index ;
	int extrap_index = (field->index)+(field->time_levels) ;
	int tC = grid->tC ;
	double p_0, p_1 ;

	for (int jC=0; jC<(grid->Nx); jC++) {
		p_0 = grid->grid_funcs[extrap_index][jC] ;
		p_1 = (
			(grid->grid_funcs[extrap_index][jC])-(grid->grid_funcs[extrap_index+1][jC])
		) / REFINEMENT ;

		grid->grid_funcs[field_index][jC] = p_0 + p_1*(tC%REFINEMENT) ;
	}	
	return ;
}
/*==========================================================================*/
void amr_extrapolate_ode_fields(amr_field* fields, amr_grid* grid) 
{
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->pde_type,"ode") == 0) {
			amr_linear_extrapapolate_field(field,grid) ;
		}
	}
	return ;
}
/*==========================================================================*/
/* use value of coarser grid at finer grid point lower boundary to  set initial condition */
/*==========================================================================*/
void set_ode_initial_condition(amr_field* fields, amr_grid* grid) 
{
	int index = 0 ;
	int perim_coord = grid->perim_coords[0] ;

	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->pde_type,"ode") == 0) {
			index = field->index ;
			printf("%s\t%d\t%d\n", field->name, index, perim_coord) ;
			grid->grid_funcs[index][0] = grid->parent->grid_funcs[index][perim_coord] ;
		}
	}
	return ; 
}
/*==========================================================================*/
void amr_solve_ode_fields(
	amr_field* fields, amr_grid* grid, void (*solve_ode)(amr_grid*)) 
{
	for (amr_grid* iter=grid; iter!=NULL; iter=iter->child) {
		if (grid->perim_interior[0]==true) {
			set_ode_initial_condition(fields, iter) ;
		}
		solve_ode(iter) ;
	}
	return ;
}
/*==========================================================================*/
/* set extrap level whenever finer grids and this grid in sync */
/*==========================================================================*/
void amr_set_extrap_levels(amr_field* field, amr_grid* grid) 
{
	int field_index = field->index ;
	int extrap_levels = field->extrap_levels ;
	int extrap_index = (field->index) + (field->time_levels) ;
	for (int jC=0; jC<(grid->Nx); jC++) {
		for (int iC=extrap_index+extrap_levels-1; iC>extrap_index; iC--) {
			grid->grid_funcs[iC][jC] = grid->grid_funcs[iC-1][jC] ;
		}
		grid->grid_funcs[extrap_index][jC] = grid->grid_funcs[field_index][jC] ;
	}	
	return ;
}
/*==========================================================================*/
void amr_set_ode_extrap_levels(amr_field* fields, amr_grid* grid) 
{
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->pde_type,"ode") == 0) {
			amr_set_extrap_levels(field, grid) ;
		}	
	}
	return ;
}
/*==========================================================================*/
/* do twice as there are two extrapolation levels */
/*==========================================================================*/
void amr_set_all_grid_ode_extrap_levels(amr_grid_hierarchy* gh) 
{
	amr_field* fields = gh->fields ;
	for (amr_grid* grid=gh->grids; grid!=NULL; grid=grid->child) {
		for (int iC=0; iC<2; iC++) {
			amr_set_ode_extrap_levels(fields, grid) ;
			amr_set_ode_extrap_levels(fields, grid) ;
		}
	}
	return ;
}
/*==========================================================================*/
/* evolves all grids in hierarchy: from coarsest to finest.
 * We solve the ODE fields as outlined in gr-qc/0508110: we extrapapolate
 * from previous solutions, then resolve the ODE over the whole hierarchy
 * when all levels align. The previous ODE level used for extrapapoaltion
 * is then reset to make it agree with linear extrapapolation with resolved
 * ODE values. */
/*==========================================================================*/
static void amr_evolve_grid(
	amr_field* fields,
	amr_grid* grid,
	int num_t_steps,
	void (*evolve_hyperbolic_pde)(amr_grid*),
	void (*solve_ode)(amr_grid*))
{
	for (int tC=0; tC<num_t_steps; tC++) {
		shift_fields_one_time_level(fields, grid) ;
		grid->tC   += 1 ;
		grid->time += grid->dt ;
		if (grid->tC%REGRID == 0) {
			/* 
				TO DO: regrid all finer levels 
			*/
		}
		if (grid->parent != NULL) {
			amr_set_interior_hyperbolic_boundary(fields, grid->parent, grid) ;
		}
		amr_extrapolate_ode_fields(fields, grid) ;
		evolve_hyperbolic_pde(grid) ;
		if (grid->child != NULL) {
			amr_evolve_grid(
				fields,
				grid->child,
				REFINEMENT,
				evolve_hyperbolic_pde,
				solve_ode)
			;
		} 	
		amr_solve_ode_fields(fields, grid, solve_ode) ;
	}
	amr_set_ode_extrap_levels(fields, grid) ;
	if (grid->parent != NULL) {
		/* 
			TO DO: compute truncation error 
		*/
		inject_overlaping_fields(fields, grid->parent, grid) ;	
	}
	return ;
}
/*==========================================================================*/
/* sets initial data on all pre assigned grid levels */
/*==========================================================================*/
static void set_initial_data(
	amr_grid_hierarchy* gh,
	void (*initial_data)(amr_grid*))
{
	amr_grid* grid = gh->grids ;
	while (grid != NULL) {
		initial_data(grid) ;
		grid = grid->child ;
	}
	shift_grids_one_time_level(gh) ;
	shift_grids_one_time_level(gh) ;
	amr_set_all_grid_ode_extrap_levels(gh) ; 
	return ;
}
/*==========================================================================*/
/* evolves all grids in hierarchy */
/*==========================================================================*/
void amr_main(
	amr_grid_hierarchy* gh, 
	void (*initial_data)(amr_grid*),
	void (*evolve_hyperbolic_pde)(amr_grid*),
	void (*solve_ode)(amr_grid*),
	void (*compute_diagnostics)(amr_grid*),
	void (*save_to_file)(amr_grid*))
{
	add_self_similar_initial_grids(gh, 2) ;
	set_initial_data(gh, initial_data) ;
	for (amr_grid* grid=gh->grids; grid != NULL; grid=grid->child) {
		save_to_file(grid) ;
	}
	for (int tC=1; tC<(gh->Nt); tC++) {
		amr_evolve_grid(
			gh->fields, 
			gh->grids,
			1,
			evolve_hyperbolic_pde,
			solve_ode) 
		;
		if (tC%(gh->t_step_save)==0) {
			for (amr_grid* grid=gh->grids; grid != NULL; grid=grid->child) {
				compute_diagnostics(grid) ;
				save_to_file(grid) ;
			}
		}
	}
	return ;	
}
