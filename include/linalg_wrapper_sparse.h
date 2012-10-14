#ifndef __WRAPPER_SPARSE__ /* allow multiple inclusions */
#define __WRAPPER_SPARSE__
#include "common.h"
#include <complex.h>
#include "evolver_common.h"
#include "sparse.h"

typedef struct {
  void *SparseNumerical;
  void *A;
  DataType Dtype;
  double PivotTolerance;
  int Factorised;
  int RefactorCount;
  int RefactorMax;
} SP_structure;


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int linalg_initialise_sparse(MultiMatrix *A,
			       EvolverOptions *options,
			       void **linalg_workspace,
			       ErrorMsg error_message);
  int linalg_finalise_sparse(void *linalg_workspace,
			     ErrorMsg error_message);
  int linalg_factorise_sparse(void *linalg_workspace,
			      int has_changed_significantly,
			      ErrorMsg error_message);
  int linalg_solve_sparse(MultiMatrix *B, 
			  MultiMatrix *X,
			  void *linalg_workspace,
			  ErrorMsg error_message);
  

#ifdef __cplusplus
}
#endif

#endif
