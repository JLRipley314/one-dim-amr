#include "amr_grid_hierarchy.h"

/*==========================================================================*/
void amr_main(
	struct amr_grid_hierarchy* gh, 
	int num_t_steps,
	int save_time,
	void (*evolve_pde)(void),
	void (*save_to_file)(void))
;
