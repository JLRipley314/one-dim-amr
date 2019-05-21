#include <math.h>

#include "stereographic_routines.h"

/*
 * Dn: n = order of derivative (1: first, 2, second, etc derivative)
 * center: center difference about grid point n
 * forward: forward ' '
 * backward: backward ' '
 * CrankNicolson: about grid point n+1/2
 */
/*==========================================================================*/
extern double stereographic_r(double s_L, double x) ;

extern double stereographic_dr(double s_L, double x, double dx) ;
/*==========================================================================*/
/* 2nd order accurate stencils */
/*==========================================================================*/
extern double D1_center_2ndOrder(double vp1, double vm1, double dx) ;

extern double D2_center_2ndOrder(double vp1, double v0 , double vm1, double dx) ;

extern double D4_center_2ndOrder(double vp2, double vp1, double v0 , double vm1, double vm2, double dx) ;

extern double D1_forward_2ndOrder(double vp2, double vp1, double v0,  double dx) ;

extern double D1_backward_2ndOrder(double v0 , double vm1, double vm2, double dx) ;

extern double D2_forward_2ndOrder(double vp3, double vp2, double vp1, double v0,  double dx) ;

extern double D2_backward_2ndOrder(double v0, double vm1, double vm2, double vm3,  double dx) ;

extern double D1_CrankNicolson_2ndOrder(double vtp1, double vt, double dt) ;

extern double AVG_CrankNicolson_2ndOrder(double v1, double v2) ;
/*==========================================================================*/
/* 2nd order stencils
 * stereographic projection: r = x / (1-x/s_L) */
/*==========================================================================*/
extern double D1_stereographic_center_2ndOrder(double x, double s_L, double vp1, double vm1, double dx) ;

extern double D1_stereographic_forward_2ndOrder(double x, double s_L, double vp2, double vp1, double v0,  double dx) ;

extern double D1_stereographic_backward_2ndOrder(double x, double s_L, double v0 , double vm1, double vm2, double dx) ;

extern double D2_stereographic_center_2ndOrder(double x, double s_L, double vp1 , double v0, double vm1, double dx) ;
