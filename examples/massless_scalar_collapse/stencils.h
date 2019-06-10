#ifndef _STENCILS_H_
#define _STENCILS_H_

#include <math.h>

/*
 * center: centered difference about grid point n
 * forward: forward ' '
 * backward: backward ' '
 * CrankNicolson: about grid point n+1/2
 */

/*==========================================================================*/
/* 2nd order accurate stencils */
/*==========================================================================*/
inline double D1_center_2ndOrder(double vp1, double vm1, double dx)
{
        return 
		(vp1 - vm1) / (2.0 * dx) ;
        ;
}
inline double D2_center_2ndOrder(double vp1, double v0 , double vm1, double dx)
{
        return 
		(vp1 - (2.*v0) + vm1) / (dx * dx) ;
        ;
}
inline double D4_center_2ndOrder(
        double vp2, double vp1, double v0 , double vm1, double vm2, double dx)
{
        return
                (vp2 + (-4.*vp1) + (6.*v0) + (-4.*vm1) + vm2) * pow(dx,-4)
        ;
}
inline double D1_forward_2ndOrder(double vp2, double vp1, double v0,  double dx)
{
        return 
		((-(vp2) + (4.*(vp1)) - (3.*(v0)))/(2.*(dx)))
        ;
}
inline double D1_backward_2ndOrder(double v0 , double vm1, double vm2, double dx)
{
        return 
		((+(vm2) - (4.*(vm1)) + (3.*(v0)))/(2.*(dx)))
        ;
}
inline double D2_forward_2ndOrder(double vp3, double vp2, double vp1, double v0,  double dx)
{
        return 
		(-(1.*(vp3)) + (4.*(vp2)) - (5.*(vp1)) + (2.*(v0))) / (dx*dx)
        ;
}
inline double D2_backward_2ndOrder(double v0, double vm1, double vm2, double vm3,  double dx)
{
        return 
		((-(1.*(vm3)) + (4.*(vm2)) - (5.*(vm1)) + (2.*(v0))) / ((dx)*(dx)))
        ;
}
inline double D1_CrankNicolson_2ndOrder(double vtp1, double vt, double dt)
{
        return 
		(((vtp1) - (vt)) / (dt))
        ;
}
inline double AVG_CrankNicolson_2ndOrder(double v1, double v2)
{
        return 
		(((v1) + (v2)) / 2.)
        ;
}
/*==========================================================================*/
/* stereographic projection: r = x / (1-x/s_L) */
/*==========================================================================*/
inline double stereographic_r(double s_L, double x)
{
        return pow(1. - (x/s_L), -1) * x ;
}
inline double stereographic_dr(double s_L, double x, double dx)
{
        return pow(1. - (x/s_L), -2) * dx ;
}
/*==========================================================================*/
/* stereographic derivatives */
/*==========================================================================*/
inline double D1_stereographic_center_2ndOrder(double s_L, double x, double vp1, double vm1, double dx)
{
        return
                pow(1-x/s_L,2) * (vp1 - vm1) / (2.0 * dx) ;
        ;
}
inline double D1_stereographic_forward_2ndOrder(double s_L, double x, double vp2, double vp1, double v0,  double dx)
{
        return
                pow(1-x/s_L,2) * ((-(vp2) + (4.*(vp1)) - (3.*(v0)))/(2.*(dx)))
        ;
}
inline double D1_stereographic_backward_2ndOrder(double s_L, double x, double v0 , double vm1, double vm2, double dx)
{
        return
                pow(1-x/s_L,2) * ((+(vm2) - (4.*(vm1)) + (3.*(v0)))/(2.*(dx)))
        ;
}
inline double D2_stereographic_center_2ndOrder(double s_L, double x, double vp1, double v0 , double vm1, double dx)
{
        return
		pow(1-x/s_L,4) * (vp1 - (2.*v0) + vm1) / (dx*dx)
	-	(2./s_L) * pow(1-x/s_L,3) * (vp1 - vm1) / (2*dx)
        ;
}
/*==========================================================================*/
/*==========================================================================*/
void Kreiss_Oliger_filter(
	int Nx,
	double* field)
;
/*==========================================================================*/
void Kreiss_Oliger_filter_origin(
	double *field,
	char *parity)
;
/*==========================================================================*/

#endif  /* _STENCILS_H_ */
