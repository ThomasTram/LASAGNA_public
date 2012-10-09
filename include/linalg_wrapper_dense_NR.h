#include "common.h"
#include <complex.h>
#include "multimatrix.h"

typedef struct {
  MultiMatrix *LU; //The LU decomposition.
  double *LUw;     //Workspace for LU
  int *luidx;      //Permutation vector for LU.
  int neq;
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
  int linalg_factorise_dense_NR(MultiMatrix *A, 
				void *linalg_workspace,
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
