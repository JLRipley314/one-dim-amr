/*===========================================================================*/
/*===========================================================================*/
#include "amr_grid_hierarchy.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*============================================================================*/
const int amr_max_levels = 5 ; 
const int refinement = 2 ; 
const int regrid = 80 ; 
const int buffer_coord = 40 ; 
const int min_grid_size = 40 ;

const double trunc_err_tolerance = 1e-7 ; 

const char HYPERBOLIC[] = "hyperbolic" ;
const char ELLIPTIC[] = "elliptic" ;
const char ODE[] = "ode" ;
const char DIAGNOSTIC[] = "diagnostic" ;

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
static inline int max_int(int val_1, int val_2) 
{
	return (val_1>val_2) ? val_1 : val_2 ;
}
static inline int min_int(int val_1, int val_2) 
{
	return (val_1<val_2) ? val_1 : val_2 ;
}
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
	assert(array[0]!=NULL) ;
        free(array[0]) ;    
        array[0] = NULL ;    
	assert(array!=NULL) ;
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
	assert((current->next)!=NULL) ;

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
	assert(grid!=NULL) ; 
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
	assert(grid!=NULL) ;
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
	assert(grid!=NULL) ;
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
inline static double cubic_polynomial(
	double p0, double p1, double p2, double p3, int kC)
{
	double step = (double)kC/(double)refinement ;
	return p0 + p1*step + (p2*pow(step,2)/2.) + (p3*pow(step,3)/6.) ;
}
/*==========================================================================*/
/*static*/ void interpolate_cubic_grid_func(
	int lower_coord, int Nx, double *parent_gf, double *child_gf)
{
	assert(parent_gf!=NULL) ;
	assert(child_gf!=NULL) ;
/* lower */
	int jC = 0 ;
	double vj   = parent_gf[jC+lower_coord] ;
	double vjp1 = parent_gf[jC+lower_coord+1] ;
	double vjp2 = parent_gf[jC+lower_coord+2] ;
	double vjp3 = parent_gf[jC+lower_coord+3] ;

	double p0 = vj ;
	double p1 = (-11*vj)/6. + 3*vjp1 - (3*vjp2)/2. + vjp3/3. ;
	double p2 = 2*vj - 5*vjp1 + 4*vjp2 - vjp3 ;
	double p3 = - vj + 3*vjp1 - 3*vjp2 + vjp3 ;
	for (int kC=0; kC<refinement; kC++) {
		child_gf[refinement*jC+kC] = cubic_polynomial(p0,p1,p2,p3,kC) ; 
		printf("%d\t%f\n", refinement*jC+kC, child_gf[refinement*jC+kC]) ;
	}
/* upper */
	jC = Nx-1 ;
	       vj   = parent_gf[jC+lower_coord] ;
	double vjm1 = parent_gf[jC+lower_coord-1] ;
	double vjm2 = parent_gf[jC+lower_coord-2] ;
	double vjm3 = parent_gf[jC+lower_coord-3] ;

	p0 = vj ;
	p1 = (11*vj)/6. - 3*vjm1 + (3*vjm2)/2. - vjm3/3. ;
	p2 = 2*vj - 5*vjm1 + 4*vjm2 - vjm3 ;
	p3 =   vj - 3*vjm1 + 3*vjm2 - vjm3 ;
	for (int kC=0; kC<refinement; kC++) {
		child_gf[refinement*jC-kC] = cubic_polynomial(p0,p1,p2,p3,-kC) ; 
		printf("%d\t%f\n", refinement*jC-kC, child_gf[refinement*jC-kC]) ;
	}
/* interior */
	for (jC=1; jC<Nx-2; jC++) {
		vjm1 = parent_gf[jC+lower_coord-1] ;
		vj   = parent_gf[jC+lower_coord] ;
		vjp1 = parent_gf[jC+lower_coord+1] ;
		vjp2 = parent_gf[jC+lower_coord+2] ;

		p0 = vj ;
		p1 = -vj/2. - vjm1/3. + vjp1 - vjp2/6. ;
		p2 = -2*vj + vjm1 + vjp1 ;
		p3 = 3*vj - vjm1 - 3*vjp1 + vjp2 ;

		for (int kC=0; kC<refinement; kC++) {
			child_gf[refinement*jC+kC] = cubic_polynomial(p0,p1,p2,p3,kC) ;
			printf("%d\t%f\n", refinement*jC+kC, child_gf[refinement*jC+kC]) ;
		}
	}
	return ;
}
/*==========================================================================*/
static void interpolate_grid_func(
	int lower_coord, int Nx, double *parent_gf, double *child_gf)
{
	assert(parent_gf!=NULL) ;
	assert(child_gf!=NULL) ;
	for (int jC=0; jC<(Nx-1); jC++) {
		double p0 = parent_gf[jC+lower_coord] ; 
		double p1 = (
			parent_gf[jC+lower_coord+1]
		-	parent_gf[jC+lower_coord]
		)/refinement ; 
		for (int kC=0; kC<refinement; kC++) {
			child_gf[refinement*jC+kC] = p0 + kC*p1 ;
		}
	}
	child_gf[refinement*(Nx-1)] = parent_gf[Nx-1+lower_coord] ;

	return ;
}
/*==========================================================================*/
static void interpolate_all_grid_funcs(amr_grid *parent, amr_grid *child)
{
	assert(parent!=NULL) ;
	assert(child!=NULL) ;
	int lower_coord = child->perim_coords[0] ;
	int upper_coord = child->perim_coords[1] ;
	int Nx = upper_coord - lower_coord + 1 ;
	for (int iC=0; iC<(parent->num_grid_funcs); iC++) {
//		interpolate_cubic_grid_func(
//			lower_coord, Nx, parent->grid_funcs[iC], child->grid_funcs[iC])
//		;
		interpolate_grid_func(
			lower_coord, Nx, parent->grid_funcs[iC], child->grid_funcs[iC])
		;
	}
	return ;
}
/*==========================================================================*/
/* Kreiss Oliger filter the grid functions */
/*==========================================================================*/
void smooth_all_grid_funcs(amr_field *fields, amr_grid *grid)
{
	double epsilon_ko = 1.0 ;
	int Nx = grid->Nx ;
	int exc_jC = grid->excised_jC ;
	for (amr_field *field=fields; field!=NULL; field=field->next) {	
		double *vals = grid->grid_funcs[field->index] ;
		for (int jC=exc_jC+2; jC<(Nx-2); jC++) {
			vals[jC] -= (epsilon_ko/16.) * (
				vals[jC+2] 
			+       (-4.*vals[jC+1])
			+       (6.*vals[jC]) 
			+       (-4.*vals[jC-1])
			+       vals[jC-2]
			)
			;
        	}
/* for outer excision boundary */
		vals[Nx-2] += (epsilon_ko/16.) * (
			vals[Nx-1] 
		+       (-4.*vals[Nx-2]) 
		+       (6.*vals[Nx-3]) 
		+       (-4.*vals[Nx-4]) 
		+       vals[Nx-5] 
		) ;
		vals[1] += (epsilon_ko/16.) * ( 
                        vals[0] 
                +       (-4.*vals[1]) 
                +       (6.*vals[2]) 
                +       (-4.*vals[3]) 
                +       vals[4] 
        	) ;
		vals = NULL ;
	}
	return ;
}
/*==========================================================================*/
/* 	if the parent grid moved, then the finer grid should have different
	bounding coordinates to not move itself in physical space. Only
	the immediate child needs to be shifted as it does not actually
	change in size and perim coords are defined with respect to
	immediate parent. */
/*==========================================================================*/
static void reset_child_grid_perim_coords(
	int old_parent_lower_coord, int new_parent_lower_coord, amr_grid *child)
{
	if (child==NULL) return ;

	int shift = refinement*(old_parent_lower_coord-new_parent_lower_coord) ;

	child->perim_coords[0] += shift ;
	child->perim_coords[1] += shift ;

	return ;
}
/*==========================================================================*/
/* inject grid functions if overlapping */
/*==========================================================================*/
static void inject_old_child_vals(amr_grid *old_child, amr_grid *new_child)
{
	assert(old_child!=NULL) ;
	assert(new_child!=NULL) ;
	int old_lower_coord = old_child->perim_coords[0] ;
	int new_lower_coord = new_child->perim_coords[0] ;
	int old_upper_coord = old_child->perim_coords[1] ;
	int new_upper_coord = new_child->perim_coords[1] ;

	int Nx = refinement*(
		min_int(old_upper_coord,new_upper_coord) - max_int(old_lower_coord,new_lower_coord)
	) + 1 
	;
	if (Nx<=0) return ;
	int offset = refinement * (old_lower_coord-new_lower_coord) ;

	if (offset>0) {
		for (int iC=0; iC<(old_child->num_grid_funcs); iC++) {
			for (int jC=0; jC<Nx; jC++) {
				new_child->grid_funcs[iC][jC+offset] = old_child->grid_funcs[iC][jC] ; 
			}
		}
	} else  {
		offset *= -1 ;
		for (int iC=0; iC<(old_child->num_grid_funcs); iC++) {
			for (int jC=0; jC<Nx; jC++) {
				new_child->grid_funcs[iC][jC] = old_child->grid_funcs[iC][jC+offset] ; 
			}
		}
	}
	return ;
}
/*==========================================================================*/
/* if cannot fit into parent grid + buffer and if grandchild (if exists)
   cannot fit into this grid then send back new/lower coords with -1 */
/*==========================================================================*/
static void make_coords_fit_into_parent_and_fit_grandchild(
		amr_grid *grid, int *new_lower_coord, int *new_upper_coord) 
{
	int Nx = grid->Nx ;
	int parent_bounding_lower_coord = 0 ;
	if (grid->perim_interior[0]==true) {
		parent_bounding_lower_coord += buffer_coord ;
	}
	int parent_bounding_upper_coord = Nx-1 ;
	if (grid->perim_interior[1]==true) {
		parent_bounding_lower_coord -= buffer_coord ;
	}
	if (parent_bounding_lower_coord>=parent_bounding_upper_coord) {
		*new_lower_coord = -1 ;
		*new_upper_coord = -1 ;
		return ;
	} 
/* if there is no grandchild then does not constrain the new grid */
	int grandchild_bounding_lower_coord = *new_lower_coord ;
	int grandchild_bounding_upper_coord = *new_upper_coord ;

	if ((grid->child!=NULL)
	&&  (grid->child->child!=NULL)
	) {
		amr_grid *grandchild = grid->child->child ;
		grandchild_bounding_lower_coord = grandchild->perim_coords[0]/refinement ;
		grandchild_bounding_upper_coord = grandchild->perim_coords[1]/refinement ;

		if (grandchild->perim_interior[0]==true) {
			grandchild_bounding_lower_coord -= buffer_coord/refinement ;
		}
		if (grandchild->perim_interior[1]==true) {
			grandchild_bounding_upper_coord += buffer_coord/refinement ;
		}
		if ((grandchild_bounding_lower_coord<parent_bounding_lower_coord) 
		||  (grandchild_bounding_upper_coord>parent_bounding_upper_coord) 
		) {
			*new_lower_coord = -1 ;
			*new_upper_coord = -1 ;
			return ; 
		}
		grandchild = NULL ;
	}
/* format:
	parent_lower_coord | new_lower_coord | grandchild_lower_coord
 and
	grandchild_upper_coord | new_upper_coord | parent_upper_coord	
*/	
	if (*new_lower_coord < parent_bounding_lower_coord) {
		*new_lower_coord = parent_bounding_lower_coord ;
	}
	if (*new_lower_coord > grandchild_bounding_lower_coord) {
		*new_lower_coord = grandchild_bounding_lower_coord ;
	}
	if (*new_upper_coord > parent_bounding_upper_coord) {
		*new_upper_coord = parent_bounding_upper_coord ;
	}
	if (*new_upper_coord < grandchild_bounding_upper_coord) {
		*new_upper_coord = grandchild_bounding_upper_coord ;
	}

	return ;
}
/*==========================================================================*/
static void add_flagged_child_grid(amr_grid *grid)
{
	assert(grid!=NULL) ;
	int new_lower_coord = grid->flagged_coords[0] ;
	new_lower_coord = 0 ;
	int new_upper_coord = grid->flagged_coords[1] ;

	amr_grid *old_child = grid->child ;
	int old_lower_coord = 0 ;
	int old_upper_coord = 0 ;
	if (old_child!=NULL) {
		old_lower_coord = old_child->perim_coords[0] ;
		old_upper_coord = old_child->perim_coords[1] ;
/* do not regrid if the new grid is too close to the old grid value 
*/
		if ((fabs(old_lower_coord-new_lower_coord)<8)
		&&  (fabs(old_upper_coord-new_upper_coord)<8)
		) {
			return ;
		}
	}
	make_coords_fit_into_parent_and_fit_grandchild(grid, &new_lower_coord, &new_upper_coord) ; 

	if (((new_upper_coord-new_lower_coord)<min_grid_size) 
	||  (new_lower_coord==-1)
	||  (new_upper_coord==-1)
	) {
		if ((old_child!=NULL)
		&& ((old_child->child)==NULL)
		) {
			amr_destroy_grid(old_child) ;
		}
		return ;
	}
	amr_grid *new_child = amr_make_finer_grid(
		new_lower_coord, new_upper_coord, grid
	) ;
	interpolate_all_grid_funcs(grid, new_child) ;

	if (old_child!=NULL) {
		inject_old_child_vals(old_child,new_child) ;
		amr_destroy_grid(old_child) ; 
	} 
	if (grid->child!=NULL) {
		printf("grid->child->level %d\n", grid->child->level) ;
	}
	amr_insert_grid(new_child, grid) ;
	if (new_child->child!=NULL) {
		printf("new_child->child->level %d\n", new_child->child->level) ;
		reset_child_grid_perim_coords(old_lower_coord, new_lower_coord, grid->child->child) ;
	}

	return ;
}
/*============================================================================*/
void amr_add_finer_grid(int left_coord, int right_coord, amr_grid* parent)
{
	assert(parent!=NULL) ;
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
	assert(gf_parent!=NULL) ;
	assert(gf!=NULL) ;
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
	assert(parent!=NULL) ;
	assert(grid!=NULL) ;
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
	assert(parent!=NULL) ;
	assert(grid!=NULL) ;
	double regrid_err_lim = trunc_err_tolerance ;
	for (amr_field *field=fields; field!=NULL; field=(field->next)) {
		if (strcmp(field->pde_type,HYPERBOLIC)!=0) {
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
	}
	return ;
}
/*==========================================================================*/
static void determine_grid_coords(
	amr_field *fields, amr_grid *grid)
{
	assert(grid!=NULL) ;
	int Nx = grid->Nx ;
	int lower_child_grid_coord = Nx-1 ; 
	int upper_child_grid_coord = 0 ; 
/* find min and max cords */
	for (amr_field *field=fields; field!=NULL; field=(field->next)) {
		int lower_coord = field->flagged_coords[0] ;
		int upper_coord = field->flagged_coords[1] ;
		if ((	(strcmp(field->pde_type,HYPERBOLIC)==0) 
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
/* if not close to physical boundary then add buffer */
	if ((grid->perim_interior[0]==true)
	&&  (lower_child_grid_coord - buffer_coord > 0)
	) {
		lower_child_grid_coord -= buffer_coord ;
	}
	if ((grid->perim_interior[1]==true)
	&&  (lower_child_grid_coord + buffer_coord < Nx-1)
	) {
		upper_child_grid_coord +=buffer_coord ;
	}
	grid->flagged_coords[0] = lower_child_grid_coord ;
	grid->flagged_coords[1] = upper_child_grid_coord ;

	return ;
}
/*==========================================================================*/
/* only regrid one grid as others have already be syncronised-so 
 * estimated truncation error would be zero. when evolving the finer
 * grid they will be regrided anyways so this is not a problem. */
/*==========================================================================*/
void regrid_all_finer_levels(amr_field *fields, amr_grid *grid)
{
	assert(grid!=NULL) ;
	flag_field_regridding_coords(fields, grid->parent, grid) ;
	determine_grid_coords(fields, grid) ;
	printf("level %d\n", grid->level) ;
	printf("lower %d\tupper %d\n", grid->flagged_coords[0], grid->flagged_coords[1]) ;
	fflush(NULL) ;
	if ((grid->level)<amr_max_levels-1) {
		add_flagged_child_grid(grid) ;
		smooth_all_grid_funcs(fields, grid) ;
	}
	return ;
}
/*============================================================================*/
/* starting from level 1 grid, add self similar grid hierarchy from origin */
/*============================================================================*/
void add_initial_grids(
	amr_grid_hierarchy* gh) 
{
	amr_grid* grid = gh->grids->child ; /* start at grid level=1 */
	int Nx = grid->Nx ; 

	int num_grids = 2 ;

	for (int iC=0; iC<num_grids; iC++) {
		switch (iC) {
			case 0:
				amr_add_finer_grid(0, (int)(Nx/3), grid) ;
				break ;
			case 1:
				amr_add_finer_grid(0, (int)(Nx-2*buffer_coord), grid) ;
				break ;	
			default:
				amr_add_finer_grid(0, (int)(Nx/2), grid) ;
				break ;
		}
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
	assert(grid_to_insert!=NULL) ;
	assert(grid!=NULL) ;
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
	assert(grid!=NULL) ;
	printf("destroying level %d\n", grid->level) ;
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
