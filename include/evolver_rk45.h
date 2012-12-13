#ifndef __RK__
#define __RK__
#include "common.h"
#include "evolver_common.h"
/**************************************************************/

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
  
int evolver_rk45(int (*derivs)(double x,double * y,double * dy,
			       void * parameters_and_workspace, ErrorMsg error_message),
		 void * ppaw,
		 double x_ini,
		 double x_final,
		 double * y_inout, 
		 size_t neq, 
		 EvolverOptions *options,
		 ErrorMsg error_message);

int evolver_rkdp45(int (*derivs)(double x,
				 double * y,
				 double * dy,
				 void * parameters_and_workspace, 
				 ErrorMsg error_message),
		   void *ppaw,
		   double t_ini,
		   double t_final,
		   double * y_inout, 
		   size_t neq, 
		   EvolverOptions *options,
		   ErrorMsg error_message);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
