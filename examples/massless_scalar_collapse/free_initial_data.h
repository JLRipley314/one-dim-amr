#ifndef _FREE_INITIAL_DATA_H_
#define _FREE_INITIAL_DATA_H_

/*==========================================================================*/
/* Gaussian pulse */
/*==========================================================================*/
void initial_data_Gaussian(
	double s_L,
	int Nx, 	
	double dt, double dx,
	double bbox[2],
	char* direction,
	double amp, double width, double center,
	double* Al, double* Ze, 
	double*  P, double*  Q)
;
/*==========================================================================*/
/* Black hole with specified mass (in Painleve-Gullstrand coordinates). */ 
/*==========================================================================*/
void initial_data_black_hole(
	double s_L,
	int Nx, 	
	double dt, double dx,
	double bbox[2],
	double mass,
	int exc_jC,
	double* Al, double* Ze)
;
#endif /* _FREE_INITIAL_DATA_H_ */
