#ifndef __RADAU__
#define __RADAU__
#include <complex.h>
#include "common.h"
#include "sparse.h"
#include "evolver_ndf15.h"
#ifdef _SUPERLU
/* Define my integer type int_t */
#include "slu_dcomplex.h"
#endif
/**************************************************************/

struct jacobian_plus{
  //For normal method:
  double complex **LU_cx;
  int *luidx;
  //For sparse method:
  int sparse_stuff_initialised;
  sp_mat_cx *spJ_cx;
  sp_num_cx *Numerical_cx;
#ifdef _SUPERLU
  doublecomplex *Axz;
  SuperMatrix A;
  SuperMatrix AC;
  SuperMatrix L;
  SuperMatrix U;
  superlumt_options_t superlumt_options;
  Gstat_t  Gstat;
  int SLU_info;
#endif
};


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
#ifdef _SUPERLU
extern void zCreate_CompCol_Matrix(SuperMatrix *, 
				   int, int, int, 
				   doublecomplex *,
				   int *, int *, 
				   Stype_t, Dtype_t, Mtype_t);

extern void zCreate_Dense_Matrix(SuperMatrix *, 
				 int, int, 
				 doublecomplex *, 
				 int,
				 Stype_t, Dtype_t, Mtype_t);
#endif
  int initialize_jacobian_plus(struct jacobian *jac, 
			       struct jacobian_plus *jac_plus, 
			       int neq,
			       ErrorMsg error_message);
  int uninitialize_jacobian_plus(struct jacobian_plus *jac_plus);
  int evolver_radau5(
		     int (*derivs)(double x,double * y,double * dy,
				   void * parameters_and_workspace, 
				   ErrorMsg error_message),
		     double x_ini,
		     double x_final,
		     double * y_inout, 
		     int * used_in_output,
		     int neq, 
		     void * parameters_and_workspace_for_derivs,
		     double rtol, 
		     double abstol, 
		     double * t_vec, 
		     int t_res,
		     int *Ap,
		     int *Ai,
		     int (*output)(double x,double y[],double dy[],
				   int index_x,void * parameters_and_workspace,
				   ErrorMsg error_message),
		     int (*print_variables)(double x, 
					    double y[], 
					    double dy[], 
					    void *parameters_and_workspace,
					    ErrorMsg error_message),
		     int (*stop_function)(double x, 
					  double y[], 
					  double dy[], 
					  void *parameters_and_workspace,
					  ErrorMsg error_message),
		     ErrorMsg error_message);
  int new_linearisation_radau5(struct jacobian *jac, 
			       struct jacobian_plus *jac_plus,
			       double hnew,
			       int neq,
			       ErrorMsg error_message);
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
			  int neq);
  double norm_inf(double *y, 
		  double *err_y, 
		  double threshold,
		  int neq);
  double norm_L2(double *y, 
		 double *err_y,
		 double threshold,
		 int neq);
  int lubksb_cx(double complex **a, int n, int *indx, double complex b[]);
  int ludcmp_cx(double complex **a, int n, int *indx, double *d, double *vv);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
