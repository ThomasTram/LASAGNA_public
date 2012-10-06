#ifndef __EVO__
#define __EVO__
#include "common.h"
#include "sparse.h"
#define TINY 1e-50
#ifdef _SUPERLU
typedef int int_t; /* default */
#include "supermatrix.h"
#include "slu_mt_util.h"
#define COLAMD_KNOBS 20
#endif
/**************************************************************/

enum lu_package{
  dense=0,
  sparse=1,
  SuperLU=2,
  superLU=2, //synonyms to catch typos.
  Superlu=2,
  superlu=2
};

struct jacobian{
  /*Stuff for normal method: */
  double **dfdy;
  double *jacvec; /*Stores experience gained from subsequent calls */
  double **LU;
  double *LUw;
  int *luidx;
  /*Sparse stuff:*/
  enum lu_package lu_pack;
  int sparse_stuff_initialized;
  int max_nonzero;     /*Maximal number of non-zero entries to be considered sparse */
  int repeated_pattern;
  int trust_sparse; /* Number of times a pattern is repeated (actually included) before we trust it. */
  int pattern_supplied;
  int refactor_count;
  int refactor_max;
  int has_grouping;
  int has_pattern;
  int new_jacobian; /* True if sp_ludcmp has not been run on the current jacobian. */
  int cnzmax;
  int *col_group; /* Column grouping. Groups go from 0 to max_group*/
  int *col_wi; /* Workarray for column grouping*/
  int max_group; /*Number of columngroups -1 */
  sp_mat *spJ; /* Stores the matrix we want to decompose */
  double *xjac; /*Stores the values of the sparse jacobian. (Same pattern as spJ) */
  sp_num *Numerical; /*Stores the LU decomposition.*/
  int *Cp; /* Stores the column pointers of the spJ+spJ' sparsity pattern. */
  int *Ci; /* Stores the row indices of the  spJ+spJ' sparsity pattern. */
#ifdef _SUPERLU
  SuperMatrix A;
  SuperMatrix AC;
  SuperMatrix L;
  SuperMatrix U;
  superlumt_options_t superlumt_options;
  Gstat_t  Gstat;
  int SLU_info;
#endif
};

struct numjac_workspace{
  /* Allocate vectors and matrices: */
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

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
#ifdef _SUPERLU
  extern void dCreate_CompCol_Matrix(SuperMatrix *, 
				     int, int, int, 
				     double *,
				     int *, int *, 
				     Stype_t, Dtype_t, Mtype_t);
  
  extern void dCreate_Dense_Matrix(SuperMatrix *, 
				   int, int, 
				   double *, 
				   int,
				   Stype_t, Dtype_t, Mtype_t);
  extern int colamd_recommended(int nnz, int n_row, int n_col) ;
  extern int colamd(int n_row, int n_col,int Alen, int A [],
		    int p [], double knobs [COLAMD_KNOBS]) ;
  extern int     sp_ienv(int);
  extern void pdgstrf (superlumt_options_t *, SuperMatrix *, int *, 
		       SuperMatrix *, SuperMatrix *, Gstat_t *, int *);
  extern void pdgstrf_init (int, fact_t, trans_t, yes_no_t, int, int, 
			    double, yes_no_t, double,
			    int *, int *, void *, int, SuperMatrix *,
			    SuperMatrix *, superlumt_options_t *, Gstat_t *);
  extern void dgstrs (trans_t, SuperMatrix *, SuperMatrix*, 
		      int*, int*, SuperMatrix*, Gstat_t *, int *);
#endif

  int initialize_jacobian(struct jacobian *jac, int neq, ErrorMsg error_message);
  int uninitialize_jacobian(struct jacobian *jac);
  int initialize_numjac_workspace(struct numjac_workspace * nj_ws,int neq, ErrorMsg error_message);
  int uninitialize_numjac_workspace(struct numjac_workspace * nj_ws);
  int calc_C(struct jacobian *jac);
  int interp_from_dif(double tinterp,double tnew,double *ynew,double h,double **dif,int k, double *yinterp,
		      double *ypinterp, double *yppinterp, int* index, int neq, int output);
  int new_linearisation(struct jacobian *jac,double hinvGak,int neq, ErrorMsg error_message);
  int adjust_stepsize(double **dif, double abshdivabshlast, int neq,int k);
  void eqvec(double *datavec,double *emptyvec, int n);
  int lubksb(double **a, int n, int *indx, double b[]);
  int ludcmp(double **a, int n, int *indx, double *d, double *vv);
  
  int numjac(int (*derivs)(double x,double * y,double * dy,void * parameters_and_workspace,ErrorMsg error_message),
	     double t, double *y, double *fval, struct jacobian *jac, struct numjac_workspace *nj_ws,
	     double thresh, int neq, int *nfe,
	     void * parameters_and_workspace_for_derivs, ErrorMsg error_message);
  
  
int evolver_ndf15(
	int (*derivs)(double x,double * y,double * dy,
		void * parameters_and_workspace, ErrorMsg error_message),
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
	int (*output)(double x,
		      double y[],
		      double dy[],
		      int index_x,
		      void * parameters_and_workspace,
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


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
