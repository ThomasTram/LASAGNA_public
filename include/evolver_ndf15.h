#ifndef __EVONDF15__
#define __EVONDF15__
#include "common.h"
#include "evolver_common.h"
/**************************************************************/


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
  int interp_from_dif(double tinterp,double tnew,double *ynew,double h,double **dif,int k, double *yinterp,
		      double *ypinterp, double *yppinterp, int* index, size_t neq, int output);
  int update_linear_system_ndf15(MultiMatrix *J, MultiMatrix *A, double hinvGak);
  int adjust_stepsize(double **dif, double abshdivabshlast, size_t neq,int k);
  void eqvec(double *datavec,double *emptyvec, int n);  
  int evolver_ndf15(int (*derivs)(double x,double * y,double * dy,
				  void * parameters_and_workspace, ErrorMsg error_message),
		    void * parameters_and_workspace_for_derivs,
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
