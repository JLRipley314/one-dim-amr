#include "amr_grid_hierarchy.h"

/*==========================================================================*/
void amr_main(
	amr_grid_hierarchy* gh, 
	void (*initial_data)(amr_grid*),
	void (*evolve_pde)(amr_grid*),
	void (*save_to_file)(amr_grid*))
;
