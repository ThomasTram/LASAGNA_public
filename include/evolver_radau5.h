#ifndef __RADAU__
#define __RADAU__
#include <complex.h>
#include "common.h"
#include "evolver_common.h"
/**************************************************************/


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
  int evolver_radau5(int (*derivs)(double x,double * y,double * dy,
				   void * parameters_and_workspace, ErrorMsg error_message),
		     void * parameters_and_workspace_for_derivs,
		     double t0,
		     double tfinal,
		     double * y_inout, 
		     size_t neq, 
		     EvolverOptions *options,
		     ErrorMsg error_message);
  int update_linear_system_radau5(MultiMatrix *J,
				  MultiMatrix *A,
				  MultiMatrix *Z,
				  double hnew);
  int transform_C_tensor_I(double *C, 
			   int s,
			   double *vec_in,
			   double *vec_out,
			   double *vec_init,
			   int n);
  int dense_output_radau5(double tinterp, 
			  double *yi, 
			  double t0, 
			  double h,
			  double *y0, 
			  double *Fi,
			  int *interpidx,
			  size_t neq);
  double norm_inf(double *y, 
		  double *err_y, 
		  double threshold,
		  size_t neq);
  double norm_L2(double *y, 
		 double *err_y,
		 double threshold,
		 size_t neq);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
