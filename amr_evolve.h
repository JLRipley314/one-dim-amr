#include "amr_grid_hierarchy.h"

#define REFINEMENT 4
#define REGRID 8

void amr_evolve(
	struct grid* grid,
	int num_t_steps,
	void (*evolve_pde)(void))
;
