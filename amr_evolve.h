#include "amr_grid_hierarchy.h"

/*==========================================================================*/
/* currently: evolves mixed hyperbolic-ode system in 1+1 dim with a
 * FIXED grid */
/*==========================================================================*/
void amr_main(
	amr_grid_hierarchy* gh, 
	void (*initial_data)(amr_grid*),
	void (*evolve_hyperbolic_pde)(amr_grid*),
	void (*solve_ode)(amr_grid*),
	void (*compute_diagnostics)(amr_grid*),
	void (*save_to_file)(amr_grid*))
;
