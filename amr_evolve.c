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
static void copy_to_2nd_array(int Nx, double *field_1, double *field_2) 
{
	for (int iC=0; iC<Nx; iC++) {
		field_2[iC] = field_1[iC] ;
	}
	return ;
}
/*==========================================================================*/
/* shifting fields from farthest back in time to most recent in time */
/*==========================================================================*/
static void shift_field(int field_index, int time_levels, amr_grid *grid)
{
	int Nx = grid->Nx ;
	for (int iC=field_index+time_levels-1; iC>field_index; iC--) {
		copy_to_2nd_array(Nx, grid->grid_funcs[iC-1], grid->grid_funcs[iC]) ;
	}
}
/*==========================================================================*/
static void shift_fields_one_time_level(
	amr_field *fields,
	amr_grid *grid)
{
	int index, time_levels ;
	for (amr_field *field=fields; field!=NULL; field=field->next) {
		index = field->index ;
		time_levels = field->time_levels ;
		shift_field(index,time_levels,grid) ;
	}
	return ;
}
/*==========================================================================*/
static void shift_grids_one_time_level(amr_grid_hierarchy *gh)
{
	for (amr_grid *grid=gh->grids; grid!=NULL; grid=grid->child) {
		shift_fields_one_time_level(gh->fields, grid) ;
	}
}
/*==========================================================================*/
/* for now its linear prolongation */
/*==========================================================================*/
/*static*/ void prolong_grid_func(
	int Nx, int perim_coord_left, double *gf, double *gf_child)
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
	int Nx, int perim_coord_left, double *gf_parent, double *gf, int *trunc_err_flags)
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
	int Nx, int perim_coord_left, double *gf_parent, double *gf)
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
		amr_field *fields, amr_grid *parent, amr_grid *grid)
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
/* {ordered field_n, field_nm1, ...} */
/*==========================================================================*/
/*static*/ void set_interior_hyperbolic_boundary_linear_interp(
	amr_field* field, amr_grid* parent, amr_grid* grid)
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
static void set_interior_hyperbolic_boundary_quad_interp(
	amr_field* field, amr_grid* parent, amr_grid* grid)
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
static void set_interior_hyperbolic_boundary(
	amr_field* fields, amr_grid* parent, amr_grid* grid)
{
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->pde_type,HYPERBOLIC)==0) {
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
static void linear_extrapapolate_level_0_field(amr_field* field, amr_grid* grid) 
{
	int field_index = field->index ;
	int extrap_index = (field_index)+(field->time_levels) ;
	double p_0, p_1 ;

	for (int jC=0; jC<(grid->Nx); jC++) {
		p_0 = grid->grid_funcs[extrap_index][jC] ;
		p_1 = (
			(grid->grid_funcs[extrap_index][jC])-(grid->grid_funcs[extrap_index+1][jC])
		) ;
		grid->grid_funcs[field_index][jC] = p_0 + p_1 ;
	}	
	return ;
}
/*==========================================================================*/
/* extrapolation for level n>0 */
/*==========================================================================*/
static void linear_extrapapolate_level_n_field(amr_field* field, amr_grid* grid) 
{
	int field_index = field->index ;
	int extrap_index = (field_index)+(field->time_levels) ;
	double p_0, p_1 ;
	int step ;
	if ((grid->tC)==(grid->parent->tC)*REFINEMENT) {
		step = REFINEMENT ;
	} else {
		step = (grid->tC)%REFINEMENT ;
	}
	for (int jC=0; jC<(grid->Nx); jC++) {
		p_0 = grid->grid_funcs[extrap_index][jC] ;
		p_1 = (
			(grid->grid_funcs[extrap_index][jC])-(grid->grid_funcs[extrap_index+1][jC])
		) / REFINEMENT ;

		grid->grid_funcs[field_index][jC] = p_0 + p_1*step ;
	}	
	return ;
}
/*==========================================================================*/
static void amr_extrapolate_ode_fields(amr_field* fields, amr_grid* grid) 
{
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->pde_type,ODE) == 0) {
			if ((grid->level)==0) {
				linear_extrapapolate_level_0_field(field,grid) ;
			} else {
				linear_extrapapolate_level_n_field(field,grid) ;
			}
		}	
	}
	return ;
}
/*==========================================================================*/
/* for now: finer grid sets boundary condition for coarser to be integrated
 * outwards */
/*==========================================================================*/
static void set_ode_initial_condition(amr_field* fields, amr_grid* grid) 
{
	int excised_jC = grid->excised_jC ;
	int child_Nx = grid->child->Nx ;
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->pde_type,ODE) == 0) {
			int index = field->index ;
			grid->grid_funcs[index][excised_jC] = grid->child->grid_funcs[index][child_Nx-1] ;
		}
	}
	return ; 
}
/*==========================================================================*/
/* solve from leftmost grid to rightmost-for now it's finer to coarser
 * as all begin from r=0. The coarser level is integrated from the finer
 * level boundary outwards */
/*==========================================================================*/
static void solve_ode_fields(
	amr_field* fields, amr_grid* grid, void (*solve_ode)(amr_grid*)) 
{
	int excised_jC = grid->excised_jC ;
	int child_lower_jC = grid->child->perim_coords[0] ;
	int child_upper_jC = grid->child->perim_coords[1] ;
	int Nx = grid->Nx ;

	if (((grid->child)==NULL) 
	||  (child_upper_jC<excised_jC)
	) {
		solve_ode(grid) ;
	} else {
		if (child_lower_jC>excised_jC) {
			grid->Nx = child_lower_jC ;
	
			solve_ode(grid) ;

			grid->Nx = Nx ;
		}
		if (child_upper_jC<(Nx-1)) {
			(grid->excised_jC) = child_upper_jC ;
			(grid->excision_on) = false ;

			set_ode_initial_condition(fields, grid) ;
			solve_ode(grid) ;

			grid->excised_jC = excised_jC ;
			(grid->excision_on) = true ;
		}
	}
	return ;
}
/*==========================================================================*/
/* set extrap level whenever finer grids and this grid in sync */
/*==========================================================================*/
static void set_extrap_levels(amr_field* field, amr_grid* grid) 
{
	int field_index = field->index ;
	int extrap_levels = field->extrap_levels ;
	int extrap_index = (field->index) + (field->time_levels) ;
	for (int jC=0; jC<(grid->Nx); jC++) {
	/*	adjust first extrapolation level for linear extrapolation
	*/
		grid->grid_funcs[extrap_index][jC] = (
			grid->grid_funcs[extrap_index+1][jC] + grid->grid_funcs[field_index][jC] 
		) / 2 ;
	/*	shift extrapolation levels 
	*/
		for (int iC=extrap_index+extrap_levels-1; iC>extrap_index; iC--) {
			grid->grid_funcs[iC][jC] = grid->grid_funcs[iC-1][jC] ;
		}
		grid->grid_funcs[extrap_index][jC] = grid->grid_funcs[field_index][jC] ;
	}	
	return ;
}
/*==========================================================================*/
/* TO DO: rescale Al (the lapse) HERE for unit evolution at spatial infty */
/*==========================================================================*/
static void set_grid_ode_extrap_levels(amr_field* fields, amr_grid* grid) 
{
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		if (strcmp(field->pde_type,ODE) == 0) {
			set_extrap_levels(field, grid) ;
		}	
	}
	return ;
}
/*==========================================================================*/
/* shift extrapolation levels once and set first extrapolation level to
 * the most recent grid value */
/*==========================================================================*/
static void shift_extrap_field(int field_index, int time_levels, int extrap_levels, amr_grid* grid)
{
	int Nx = grid->Nx ;
	for (int iC=(field_index+time_levels+extrap_levels-1); iC>(field_index+time_levels); iC--) {
		copy_to_2nd_array(Nx, grid->grid_funcs[iC-1], grid->grid_funcs[iC]) ;
	}
	copy_to_2nd_array(Nx, grid->grid_funcs[field_index], grid->grid_funcs[field_index+time_levels]) ;
}
/*==========================================================================*/
static void shift_extrap_fields_one_level(
	amr_field* fields,
	amr_grid* grid)
{
	int index, time_levels, extrap_levels ;
	for (amr_field* field=fields; field!=NULL; field=field->next) {
		index = field->index ;
		time_levels = field->time_levels ;
		extrap_levels = field->extrap_levels ;
		if ((field->extrap_levels)!=0) {
			shift_extrap_field(index,time_levels,extrap_levels,grid) ;
		}
	}
	return ;
}
/*==========================================================================*/
static void shift_grids_ode_extrap_levels(amr_grid_hierarchy* gh)
{
	for (amr_grid* grid=gh->grids; grid!=NULL; grid=grid->child) {
		shift_extrap_fields_one_level(gh->fields, grid) ;
	}
}
/*==========================================================================*/
static void solve_ode_initial_data(
	amr_grid_hierarchy* gh, void (*solve_ode)(amr_grid*)) 
{
	amr_grid* grid = gh->grids ;
	amr_set_to_tail(&grid) ;
	do {
		solve_ode_fields(gh->fields, grid, solve_ode) ;
		inject_overlaping_fields(gh->fields, grid->parent, grid) ;
		grid = grid->parent ;
	} while ((grid->level)!=0) ;	

	return ;
}
/*==========================================================================*/
/* Recursive routine to evolve amr (or fmr if no regridding) hierarchy
 * 
 * Hyperbolics solved as in Berger&Oliger: coarse step onwards to finer
 * levels, and inject on overlapping grids.
 *
 * We solve the ODE fields as outlined in gr-qc/0508110: we extrapapolate
 * from previous solutions, then resolve the ODE over the whole hierarchy
 * when all levels align. The previous ODE level used for extrapapoaltion
 * is then reset to make it agree with linear extrapapolation with resolved
 * ODE values. 
 */
/*==========================================================================*/
static void evolve_grid(
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
		if ((grid->tC)%REGRID == 0) {
			/* 
				TO DO: regrid all finer levels 
			*/
		}
/* 	do not interpolate coarse grid (level 1) and  shadow (level 0),
	both which span the entire domain so boundary conditions are physical 
*/	
		if ((grid->level)>1) {
			set_interior_hyperbolic_boundary(fields, grid->parent, grid) ;
		}
		amr_extrapolate_ode_fields(fields, grid) ;
		evolve_hyperbolic_pde(grid) ;
		if ((grid->child)!=NULL) {
			evolve_grid(
				fields,
				grid->child,
				REFINEMENT,
				evolve_hyperbolic_pde,
				solve_ode)
			;
		} 
		if ((grid->level)>0) { /* all grids finer than shadow grid */
			solve_ode_fields(fields, grid, solve_ode) ;
		}
	}
	if ((grid->child)!=NULL) {
		/* 
			TO DO: compute truncation error 
		*/
		inject_overlaping_fields(fields, grid, grid->child) ;	
	}
	set_grid_ode_extrap_levels(fields, grid) ;
	return ;
}
/*==========================================================================*/
/* sets free initial data on all pre assigned grid levels - does NOT
 * solve the constraints */
/*==========================================================================*/
static void set_free_initial_data(
	amr_grid_hierarchy* gh,
	void (*free_initial_data)(amr_grid*))
{
	amr_grid* grid = gh->grids ;
	while (grid != NULL) {
		free_initial_data(grid) ;
		grid = grid->child ;
	}
	return ;
}
/*==========================================================================*/
static void save_all_grids(
	amr_grid_hierarchy* gh, void (*save_to_file)(amr_grid*))
{
	for (amr_grid* grid=gh->grids; grid != NULL; grid=grid->child) {
		save_to_file(grid) ;
	}
	return ;
}
/*==========================================================================*/
/* compute diagnostics on each grid, and syncronize the excision point
 * with the shadow grid */
/*==========================================================================*/
void set_excision_point(amr_grid* grid)
{
	if (((grid->level) > 0)
	&&  ((grid->parent->excised_jC) > (grid->perim_coords[0]))
	) {
		if ((grid->parent->excised_jC) < (grid->perim_coords[1])) {
			grid->excised_jC = REFINEMENT * (
				grid->parent->excised_jC - grid->perim_coords[0]
			) ;
		} else {
			grid->excised_jC = grid->Nx-1 ;
		}
	}
	return ;
}
/*==========================================================================*/
static void compute_all_grid_diagnostics(
	amr_grid_hierarchy* gh, void (*compute_diagnostics)(amr_grid*))
{
	for (amr_grid* grid=gh->grids; grid!=NULL; grid=grid->child) {
		compute_diagnostics(grid) ;
		set_excision_point(grid) ;
	}
	return ;
}
/*==========================================================================*/
static void set_past_t_data_first_order(amr_grid_hierarchy* gh) 
{
	shift_grids_one_time_level(gh) ;
	shift_grids_one_time_level(gh) ;
	shift_grids_ode_extrap_levels(gh) ; 
	shift_grids_ode_extrap_levels(gh) ; 
	return ;
}
/*==========================================================================*/
/* flip the extrapolation levels about the grid at level t */
/*==========================================================================*/
static void flip_earlier_time_levels(amr_field* field, amr_grid* grid)
{
	int time_levels= field->time_levels ; 
	int extrap_levels= field->extrap_levels ; 
	int field_index= field->index ;
	for (int jC=0; jC<(grid->Nx); jC++) {
		for (int iC=1; iC<(time_levels+extrap_levels); iC++) {
			
			grid->grid_funcs[field_index+iC][jC] 
			= 	2*(grid->grid_funcs[field_index][jC]) 
			- 	(grid->grid_funcs[field_index+iC][jC]) 
			;
		}
	}	
	return ;
}
/*==========================================================================*/
static void flip_dt(amr_field* fields, amr_grid* grids)
{
	for (amr_grid* grid=grids; grid!=NULL; grid=grid->child) {
		for (amr_field* field=fields; field!=NULL; field=field->next) {
			flip_earlier_time_levels(field, grid) ;
		}
		grid->dt *= -1 ;
	}
	return ;
}
/*==========================================================================*/
/* initial data: the ode are assumed to set the "constrained" degrees
 * of freedom. details of this algorithm can be found in the appendix of
 * gr-qc/0508110 */
/*==========================================================================*/
static void set_initial_data(
	amr_grid_hierarchy* gh, 
	void (*free_initial_data)(amr_grid*),
	void (*evolve_hyperbolic_pde)(amr_grid*),
	void (*solve_ode)(amr_grid*))
{
	set_free_initial_data(gh, free_initial_data) ;
	solve_ode_initial_data(gh, solve_ode) ; 	
	set_past_t_data_first_order(gh) ;	
	for (int iC=0; iC<1; iC++) {
		evolve_grid(
			gh->fields, 
			gh->grids,
			2,
			evolve_hyperbolic_pde,
			solve_ode) 
		;
		flip_dt(gh->fields,gh->grids) ;
		evolve_grid(
			gh->fields, 
			gh->grids,
			2,
			evolve_hyperbolic_pde,
			solve_ode) 
		;
		flip_dt(gh->fields,gh->grids) ;
	}
	set_free_initial_data(gh, free_initial_data) ;
	solve_ode_initial_data(gh, solve_ode) ; 	

	return ;
}
/*==========================================================================*/
/* evolves all grids in hierarchy */
/*==========================================================================*/
void amr_main(
	amr_grid_hierarchy* gh, 
	void (*free_initial_data)(amr_grid*),
	void (*evolve_hyperbolic_pde)(amr_grid*),
	void (*solve_ode)(amr_grid*),
	void (*compute_diagnostics)(amr_grid*),
	void (*save_to_file)(amr_grid*))
{
	int compute_diagnostics_tC = 1/(gh->cfl_num) ;
/* 
	add_initial...: for the fixed amr grid hierarchy 
*/
	int grid_size_ratio = 2 ;
	int grid_levels = 2 ;
	add_self_similar_initial_grids(gh, grid_size_ratio, grid_levels) ;

	set_initial_data(
		gh, 
		free_initial_data, 
		evolve_hyperbolic_pde,
		solve_ode)
	;
	compute_all_grid_diagnostics(gh, compute_diagnostics) ;
	save_all_grids(gh, save_to_file) ;

	for (int tC=1; tC<(gh->Nt); tC++) {
		evolve_grid(
			gh->fields, 
			gh->grids,
			1,
			evolve_hyperbolic_pde,
			solve_ode) 
		;
		if ((tC%(compute_diagnostics_tC)==0) 
		||  (tC%(gh->t_step_save)==0)
		) {
			compute_all_grid_diagnostics(gh, compute_diagnostics) ;
		}
		if (tC%(gh->t_step_save)==0) {
			save_all_grids(gh, save_to_file) ;
		}
	}
	return ;	
}
