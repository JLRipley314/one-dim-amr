#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/*===========================================================================*/
/* basic matrix computations */
/*===========================================================================*/
double compute_2dMatrix_determinant(double matrix[2][2])
{
	return (matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1]) ;
}
/*===========================================================================*/
double compute_2dMatrix_trace(double matrix[2][2])
{
	return (matrix[0][0] + matrix[1][1]) ;
}
/*===========================================================================*/
void compute_2dMatrix_inverse(double matrix[2][2], double inverse[2][2])
{
	double determinant = compute_2dMatrix_determinant(matrix) ;

	if (fabs(determinant) <= 1e-16) {
		printf("ERROR: compute_2dMatrix_inverse: fabs(det)<1e-16\n") ;
		exit(EXIT_FAILURE) ;
	}

	inverse[0][0] =   matrix[1][1] / determinant ;
	inverse[1][0] = - matrix[1][0] / determinant ;	
	inverse[0][1] = - matrix[0][1] / determinant ; 
	inverse[1][1] =   matrix[0][0] / determinant ; 

	return ;
}
/*===========================================================================*/
void compute_2dMatrix_multiply(
	double matrix_1[2][2], double matrix_2[2][2], double matrix_multiple[2][2])
{
	matrix_multiple[0][0] = (matrix_1[0][0]*matrix_2[0][0]) + (matrix_1[0][1]*matrix_2[1][0]) ;
	matrix_multiple[1][0] = (matrix_1[1][0]*matrix_2[0][0]) + (matrix_1[1][1]*matrix_2[1][0]) ;
	matrix_multiple[0][1] = (matrix_1[0][0]*matrix_2[0][1]) + (matrix_1[0][1]*matrix_2[1][1]) ;
	matrix_multiple[1][1] = (matrix_1[1][0]*matrix_2[0][1]) + (matrix_1[1][1]*matrix_2[1][1]) ;

	return ;
}
/*===========================================================================*/
void compute_2dMatrix_subtract(
	double matrix_1[2][2], double matrix_2[2][2], double matrix_subtract[2][2])
{
	matrix_subtract[0][0] = matrix_1[0][0] - matrix_2[0][0] ;
	matrix_subtract[1][0] = matrix_1[1][0] - matrix_2[1][0] ;
	matrix_subtract[0][1] = matrix_1[0][1] - matrix_2[0][1] ;
	matrix_subtract[1][1] = matrix_1[1][1] - matrix_2[1][1] ;

	return ;
}
/*===========================================================================*/
void compute_2dMatrix_add(
	double matrix_1[2][2], double matrix_2[2][2], double matrix_add[2][2])
{
	matrix_add[0][0] = matrix_1[0][0] + matrix_2[0][0] ;
	matrix_add[1][0] = matrix_1[1][0] + matrix_2[1][0] ;
	matrix_add[0][1] = matrix_1[0][1] + matrix_2[0][1] ;
	matrix_add[1][1] = matrix_1[1][1] + matrix_2[1][1] ;

	return ;
}
/*===========================================================================*/
/* returns true, false if real. side effect: eigenvales_2[iC]  */
/*===========================================================================*/
bool compute_2dMatrix_eigenvalues(
	double matrix_2[2][2], double eigenvalues_2[2])
{
	double trace = compute_2dMatrix_trace(matrix_2) ;
	double det   = compute_2dMatrix_determinant(matrix_2) ;
	double discriminant = pow(trace,2) - (4.*det) ;

	if (discriminant < 0) {
		eigenvalues_2[0] = 0. ; 
		eigenvalues_2[1] = 0. ; 

		return true ;
	} 
	if (trace < -1e-16) {
		eigenvalues_2[1] = 0.5 * (trace - sqrt(discriminant)) ;
		eigenvalues_2[0] = det / eigenvalues_2[1] ;
	} else if (trace > 1e-16) {
		eigenvalues_2[0] = 0.5 * (trace + sqrt(discriminant)) ;
		eigenvalues_2[1] = det / eigenvalues_2[0] ;
	} else {
		eigenvalues_2[0] =   0.5 * sqrt(discriminant) ;
		eigenvalues_2[1] = - 0.5 * sqrt(discriminant) ;
	}	
	return false ;
}
