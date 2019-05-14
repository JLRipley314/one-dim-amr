#include "amr_grid_hierarchy.h"

/*==========================================================================*/
void amr_main(
	struct amr_grid_hierarchy* gh, 
	void (*initial_data)(struct amr_grid*),
	void (*evolve_pde)(struct amr_grid*),
	void (*save_to_file)(struct amr_grid*))
;
