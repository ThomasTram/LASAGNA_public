#ifndef __NWT__
#define __NWT__
#include "common.h"
#include "evolver_ndf15.h"
/**************************************************************/

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
int Newton(int (*vecfun)(double * y, double * Fy, void *param),
	   double *y0,
	   void *param,
	   double *maxstep,
	   double rtol,
	   int *iter,
	   int maxiter,
	   int neq,
	   ErrorMsg error_message);
int jacobian_for_Newton(int (*vecfun)(double * y, double * Fy, void *param),
			double *y0,
			double *Fval,
			int neq,
			void *param,
			double **jac);


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
