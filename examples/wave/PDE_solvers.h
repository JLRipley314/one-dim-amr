/* we set Dirichlet boundary conditions for now */

void advance_tStep_wave(
	int Nx,
	double dt, 	double dx,
	double* P_n, 	double* P_nm1,
	double* Q_n,	double* Q_nm1)
;

void Kreiss_Oliger_Dissipation(
	int Nx,
	double* field)
;
