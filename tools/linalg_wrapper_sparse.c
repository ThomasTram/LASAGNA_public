#include "linalg_wrapper_sparse.h"
#include "complex.h"


/** Wrapper for linear algebra methods used by evolver_ndf15 and
    evolver_radau5. The following subroutines are called from
    the evolver:

    int linalg_initialise(MultiMatrix *A, 
                          void **linalg_workspace,
                          EvolverOptions *options,
			  ErrorMsg error_message)
    -- This is called once in the beginning for each matrix. 
       Use for allocating memory, do symbolic manipulations 
       of sparse pattern (if any), etc. The include file should 
       define a structure containing all arrays, and this structure 
       should be allocated and pointer to it returned.

       For sparse methods, the sparsity pattern should be in
       A, and a column-preorderings should be made. In the case of
       RADAU5, the sparsity patterns of the two matrices are identical,
       so I do not have to calculate the preordering again. However,
       it is not costly, and this way I keep it more general.

    int linalg_finalise(void *linalg_workspace, ErrorMsg error_message)
    -- This is called once, in the end of the computation.

    int linalg_factorise(MultiMatrix *A, 
                         void *linalg_workspace,
                         ErrorMsg error_message)
    -- Factorise A and store the factors in linalg_workspace. If
       the method can take advantage of subsequent factorisations,
       there should be a flag in linalg_workspace indicating if
       linalg_factorise has been called before.

    int linalg_solve(MultiMatrix *B, 
                     MultiMatrix *X,
                     void *linalg_workspace,
		     ErrorMsg error_message)
    -- Solve A X^T = B^T. X and B are dense, and we using row-major
       storage which explains the transpositions, ^T.

    The names of these 4 functions could in principle be anything, but
    I suggest to use the names above and append _"NAME_OF_WRAPPER".
*/
int linalg_initialise_sparse(MultiMatrix *A, 
			     EvolverOptions *options,
			     void **linalg_workspace,
			     ErrorMsg error_message){
  SCCformat *Store=A->Store;
  SP_structure *ws;
  sp_mat *spmat_dbl;
  sp_mat_cx *spmat_dbl_cx;
  int nnz;
  int *Cp, *Ci;
  int *perm_c, *wamd;

  ncol = A->ncol; nrow = A->nrow; Dtype = A->Dtype;
  //Test input:
  lasagna_test(A->ncol != A->nrow, 
	       error_message, 
	       "Matrix not square!");
  lasagna_test((A->Dtype!=L_DBL)&&(A->Dtype!=L_DBL_CX), 
	       error_message,
	       "Unknown datatype in A.");
  lasagna_test(A->Stype!=L_SCC, 
	       error_message,
	       "This wrapper only supports sparse input matrix.");
  
  nnz = Store->nnz;
  lasagna_alloc(ws,sizeof(SP_structure),error_message);
  switch (A->Dtype){
  case (L_DBL):
    lasagna_call(sp_num_alloc(&(ws->SparseNumerical), 
			      ncol, 
			      error_message),
		 error_message,error_message);
    perm_c = ((sp_num *) ws->SparseNumerical)->q;
    wamd = ((sp_num *) ws->SparseNumerical)->wamd;
    lasagna_alloc(ws->A,sizeof(sp_mat *),error_message);
    sp_mat = (sp_mat *) ws->A;
    sp_mat->ncols = ncol;
    sp_mat->nrows = nrow;
    sp_mat->maxnz = nnz;
    sp_mat->Ai = Store->Ai;
    sp_mat->Ap = Store->Ap;
    sp_mat->Ax = (double *) Store->Ax;
    break;
  case (L_DBL_CX):
    lasagna_call(sp_num_alloc_cx(&(ws->SparseNumerical), 
				 ncol, 
				 error_message),
		 error_message,error_message);
    perm_c = ((sp_num_cx *) ws->SparseNumerical)->q;
    wamd = ((sp_num_cx *) ws->SparseNumerical)->wamd;
    lasagna_alloc(ws->A,sizeof(sp_mat_cx *),error_message);
    sp_mat = (sp_mat_cx *) ws->A;
    sp_mat->ncols = ncol;
    sp_mat->nrows = nrow;
    sp_mat->maxnz = nnz;
    sp_mat->Ai = Store->Ai;
    sp_mat->Ap = Store->Ap;
    sp_mat->Ax = (double complex*) Store->Ax;
    break;
  }
  //Calculate sparsity pattern of C = A + A^T for use with AMD:
  lasagna_call(get_pattern_A_plus_AT(Store->Ap, 
				     Store->Ai, 
				     ncol, 
				     &(Cp), 
				     &(Ci), 
				     error_message), 
	       error_message,error_message);
  
  /* Calculate the optimal ordering: */
  sp_amd(Cp, Ci, ncol, Cp[ncol],perm_c,wamd);
  free(Cp);
  free(Ci);
  /** Set options for solver: */
  ws->PivotTolerance = 0.1;
  ws->RefactorCount = 0;
  ws->RefactorMax = 10;
  ws->Factorised = _FALSE_;
  *linalg_workspace = (void *) ws;
  
  return _SUCCESS_;
}

int linalg_finalise_sparse(void *linalg_workspace,
			   ErrorMsg error_message){
  
  SP_structure *ws=linalg_workspace;  
  switch (A->Dtype){
  case (L_DBL):
    sp_num_free((sp_num *) ws->SparseNumerical);
    break;
  case (L_DBL_CX):
    sp_num_free_cx((sp_num_cx *) ws->SparseNumerical);
    break;
  }
  free(ws->A);
  free(ws);
}

int linalg_factorise_sparse(void *linalg_workspace,
			    ErrorMsg error_message){
  SP_structure *ws= linalg_workspace;
  int fr;
  /** The actual data in ws->A is the same as in the MultiMatrix
      A which was passed to linalg_initialise_sparse.
  */

  switch(A->Dtype){
  case (L_DBL):
    if ((ws->Factorised==_TRUE_)&&(ws->RefactorCount < ws->RefactorMax)){
      fr = sp_refactor((sp_num *) ws->SparseNumerical, 
		       (sp_mat *) ws->A);
      ws->RefactorCount++;
    }
    else{
      fr = sp_ludcmp((sp_num *) ws->SparseNumerical, 
		     (sp_mat *) ws->A, 
		     ws->PivotTolerance);
      ws->Factorised = _TRUE_;
      ws->RefactorCount = 0;
    }
    break;
  case (L_DBL_CX):
    if ((ws->Factorised==_TRUE_)&&(ws->RefactorCount < ws->RefactorMax)){
      fr = sp_refactor_cx((sp_num_cx *) ws->SparseNumerical, 
			  (sp_mat_cx *) ws->A);
      ws->RefactorCount++;
    }
    else{
      fr = sp_ludcmp_cx((sp_num_cx *) ws->SparseNumerical, 
			(sp_mat_cx *) ws->A, 
			ws->PivotTolerance);
      ws->Factorised = _TRUE_;
      ws->RefactorCount = 0;
    }
    break;
  }
  return fr;
}

int linalg_solve_sparse(MultiMatrix *B, 
			MultiMatrix *X,
			void *linalg_workspace,
			ErrorMsg error_message){
  SP_structure *ws= linalg_workspace;
  DNRformat *StoreB=B->Store;
  DNRformat *StoreX=X->Store;
  double **MatB_dbl;
  double **MatX_dbl;
  double complex **MatB_dbl_cx;
  double complex **MatX_dbl_cx;
  int i,fr;

  switch(B->Dtype){
  case (L_DBL):
    for(i=1; i<=B->nrow; i++){ 
      MatX_dbl = (double **) StoreX->Matrix;
      MatB_dbl = (double **) StoreB->Matrix;
      //Compensate for zero indexing scheme and solve:
      fr = sp_lusolve((sp_num *) ws->SparseNumerical, 
		      MatB_dbl[i]+1, MatX_dbl[i]+1);
    }	
    break;
  case (L_DBL_CX):
    for(i=1; i<=B->nrow; i++){ 
      MatX_dbl_cx = (double complex **) StoreX->Matrix;
      MatB_dbl_cx = (double complex **) StoreB->Matrix;
      //Compensate for zero indexing scheme and solve:
      fr = sp_lusolve_cx((sp_num_cx *) ws->SparseNumerical, 
			 MatB_dbl_cx[i]+1, MatX_dbl_cx[i]+1);
    }      
    break;
  }
  return fr;
}
