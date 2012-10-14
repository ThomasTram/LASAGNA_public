#ifndef __WRAPPER_SUPERLU__ /* allow multiple inclusions */
#define __WRAPPER_SUPERLU__

#include "common.h"
#include <complex.h>
#include "evolver_common.h"
typedef int int_t; /* default */
#include "supermatrix.h"
#include "slu_mt_util.h"


typedef struct {
  SuperMatrix A;
  SuperMatrix AC;
  SuperMatrix L;
  SuperMatrix U;
  SuperMatrix X;
  superlumt_options_t superlumt_options;
  Gstat_t  Gstat;
  int SLU_info;
  int neq;
  int *perm_r;
  int *perm_c;
} SLU_structure;


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int linalg_initialise_SuperLU(MultiMatrix *A,
				EvolverOptions *options,
				void **linalg_workspace,
				ErrorMsg error_message);
  int linalg_finalise_SuperLU(void *linalg_workspace,
			      ErrorMsg error_message);
  int linalg_factorise_SuperLU(void *linalg_workspace,
			       int has_changed_significantly,
			       ErrorMsg error_message);
  int linalg_solve_SuperLU(MultiMatrix *B, 
			   MultiMatrix *X,
			   void *linalg_workspace,
			   ErrorMsg error_message);
  

#ifdef __cplusplus
}
#endif



#endif
