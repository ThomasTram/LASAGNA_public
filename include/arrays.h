#include "common.h"

#define _SPLINE_NATURAL_ 0 /**< natural spline: ddy0=ddyn=0 */
#define _SPLINE_EST_DERIV_ 1 /**< spline with estimation of first derivative on both edges */

int arrays_spline_natural(double *x, double *y, double *d2y, int n);
int arrays_spline_interpolate(double *xa, 
			      double *ya, 
			      double *d2ya, 
			      int n,
			      int xdir,
			      double x,
			      double *y,
			      int *ptr_klo,
			      int *ptr_khi);
