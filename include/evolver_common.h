#ifndef __EVO_COMMON__
#define __EVO_COMMON__

#include "common.h"
#define TINY 1e-50
#include "multimatrix.h"

typedef struct _EvolverOptions{
  int tres;             /**If t_vec!=NULL, length of t_vec.
			   If t_vec==NULL, refinement factor.*/
  double AbsTol;
  double RelTol;
  int * used_in_output; /**Which indices are required in output? */
  double * t_vec;       /**Output at specified points */
  int Stats[10]; /** Evolver statistics */
  int Flags[10]; /** Can for instance be used for communication between evolver
		     and linalg wrapper, or different instances of the wrapper. */
  int Cores;     /** Number of cores available for linalg_wrapper */
  int EvolverVerbose;  /** Level of output from evolver. */
  int LinAlgVerbose;   /** Level of output from linalg wrapper. */
  /** Pointers to the 4 linear algebra wrapper functions: */
  int (*linalg_initialise)(MultiMatrix *, struct _EvolverOptions *, void **, ErrorMsg);
  int (*linalg_finalise)(void *, ErrorMsg);
  int (*linalg_factorise)(void *, int, ErrorMsg);
  int (*linalg_solve)(MultiMatrix *, MultiMatrix *, void *, ErrorMsg);
  /** Pointers to the evolver utility functions:*/
  int (*output)(double t, double *y, double *dy, int i, void *p, ErrorMsg err);
  int (*print_variables)(double t, double *y, double *dy, void *p, ErrorMsg err);
  int (*stop_function)(double t, double *y, double *dy, void *p, ErrorMsg err);
  /** Jacobian specific stuff:*/
  int use_sparse;
  int *Ai;
  int *Ap;
  /** This flags makes the evolver return a pointer to the jacobian and 
      calculate numjac twice each time it is called.**/
  int J_pointer_flag;
  MultiMatrix *J_pointer;
} EvolverOptions;

struct numjac_workspace{
  /* Allocate vectors and matrices: */
  double *jacvec;
  int *col_group;
  int max_group; /*Number of columngroups -1 */
  double *yscale;
  double *del;
  double * Difmax;
  double * absFdelRm;
  double * absFvalue;
  double * absFvalueRm;
  double * Fscale;
  double * ffdel;
  double * yydel;
  double * tmp;

  double **ydel_Fdel;

  int * logj;
  int * Rowmax;
};

  int initialize_numjac_workspace(MultiMatrix *J, void ** numjac_workspace, ErrorMsg error_message);
  int uninitialize_numjac_workspace(void * numjac_workspace);
  int numjac(int (*derivs)(double x,double * y,double * dy,void * parameters_and_workspace,ErrorMsg error_message),
	     double t, double *y, double *fval, MultiMatrix *J, void* numjac_workspace,
	     double thresh, int neq, int *nfe,
	     void * parameters_and_workspace_for_derivs, ErrorMsg error_message);



/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif



#ifdef __cplusplus
}
#endif



#endif
