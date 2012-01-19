#ifndef __RK__
#define __RK__
#include "common.h"
/**************************************************************/

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
  
int evolver_rk45(
	int (*derivs)(double x,double * y,double * dy,
		void * parameters_and_workspace, ErrorMsg error_message),
	double x_ini,
	double x_final,
	double * y_inout, 
 	int * used_in_output,
	int neq, 
	void * parameters_and_workspace_for_derivs,
	double rtol, 
	double minimum_variation, 
	double * t_vec, 
	int t_res,
	int (*output)(double x,double y[],double dy[],int index_x,void * parameters_and_workspace,
		ErrorMsg error_message),
	int (*print_variables)(double x, double y[], double dy[], void *parameters_and_workspace,
		ErrorMsg error_message),
	ErrorMsg error_message);

int evolver_rkdp45(
		   int (*derivs)(double x,
				 double * y,
				 double * dy,
				 void * parameters_and_workspace, 
				 ErrorMsg error_message),
		   double t_ini,
		   double t_final,
		   double * y_inout, 
		   int * used_in_output,
		   int neq, 
		   void * ppaw,
		   double rtol, 
		   double abstol, 
		   double * t_vec, 
		   int tres,
		   int *ignore1,
		   int *ignore2,
		   int (*output)(double x,
				 double y[],
				 double dy[],
				 int index_x,
				 void * parameters_and_workspace,
				 ErrorMsg error_message),
		   int (*print_variables)(double x, 
					  double y[], 
					  double dy[], 
					  void *ppaw,
					  ErrorMsg error_message),
		   int (*stop_function)(double x, 
					double y[], 
					double dy[], 
					void *parameters_and_workspace,
					ErrorMsg error_message),
		   ErrorMsg error_message);


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
