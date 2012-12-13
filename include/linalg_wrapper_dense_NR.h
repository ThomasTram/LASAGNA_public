#ifndef __WRAPPER_DENSE_NR__ /* allow multiple inclusions */
#define __WRAPPER_DENSE_NR__

#include "common.h"
#include <complex.h>
#include "evolver_common.h"

typedef struct {
  MultiMatrix *A;  //Pointer to MultiMatrix A
  MultiMatrix *LU; //The LU decomposition.
  double *LUw;     //Workspace for LU
  int *luidx;      //Permutation vector for LU.
  size_t neq;
} DNR_structure;


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int linalg_initialise_dense_NR(MultiMatrix *A, 
				 EvolverOptions *options,
				 void **linalg_workspace,
				 ErrorMsg error_message);
  int linalg_finalise_dense_NR(void *linalg_workspace,
			       ErrorMsg error_message);
  int linalg_factorise_dense_NR(void *linalg_workspace,
				int has_changed_significantly,
				ErrorMsg error_message);
  int linalg_solve_dense_NR(MultiMatrix *B, 
			    MultiMatrix *X,
			    void *linalg_workspace,
			    ErrorMsg error_message);
  
  int ludcmp(double **a, int n, int *indx, double *d, double *vv);
  int ludcmp_cx(double complex **a, int n, int *indx, double *d, double *vv);
  int lubksb(double **a, int n, int *indx, double b[]);
  int lubksb_cx(double complex **a, int n, int *indx, double complex b[]);

#ifdef __cplusplus
}
#endif

#endif
