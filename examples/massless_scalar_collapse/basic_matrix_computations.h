#ifndef _BASIC_MATRIX_COMPUTATIONS_H_
#define _BASIC_MATRIX_COMPUTATIONS_H_

#include <math.h>
#include <stdbool.h>
/*===========================================================================*/
/* basic matrix computations */
/*===========================================================================*/
double compute_2dMatrix_determinant(
	double matrix[2][2])
;
/*===========================================================================*/
double compute_2dMatrix_trace(
	double matrix[2][2])
;
/*===========================================================================*/
void compute_2dMatrix_inverse(
	double matrix[2][2], double inverse[2][2])
;
/*===========================================================================*/
void compute_2dMatrix_multiply(
	double matrix_1[2][2], double matrix_2[2][2], double matrix_multiple[2][2])
;
/*===========================================================================*/
void compute_2dMatrix_subtract(
	double matrix_1[2][2], double matrix_2[2][2], double matrix_subtract[2][2])
;
/*===========================================================================*/
void compute_2dMatrix_add(
	double matrix_1[2][2], double matrix_2[2][2], double matrix_add[2][2])
;
/*===========================================================================*/
/* returns 0 if imaginary, returns 1 if real. side effect: eigenvales_2[iC]  */
/*===========================================================================*/
bool compute_2dMatrix_eigenvalues(
	double matrix_2[2][2], double eigenvalues_2[2])
;

#endif /* _BASIC_MATRIX_COMPUTATIONS_H_ */
