#ifndef __NWT__
#define __NWT__
#include "common.h"
#include "linalg_wrapper_dense_NR.h"
/**************************************************************/

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
int Newton(int (*vecfun)(double * y, double * Fy, void *param),
	   int (*jacfun)(double * y, double **jac, void *param),
	   double *y0,
	   void *param,
	   double *maxstep,
	   double rtol,
	   int *iter,
	   int maxiter,
	   size_t neq,
	   ErrorMsg error_message);
int jacobian_for_Newton(int (*vecfun)(double * y, double * Fy, void *param),
			double *y0,
			double *Fval,
			size_t neq,
			void *param,
			double **jac);


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
