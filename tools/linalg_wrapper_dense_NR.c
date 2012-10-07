#include "linalg_wrapper_dense_NR"
/** Wrapper for linear algebra methods used by evolver_ndf15 and
    evolver_radau5. The following subroutines are called from
    the evolver:

    int linalg_initialise(MultiMatrix *A, 
                          void **linalg_workspace,
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
                         void **linalg_workspace,
                         ErrorMsg error_message)
    -- Factorise A and store the factors in linalg_workspace. If
       the method can take advantage of subsequent factorisations,
       there should be a flag in linalg_workspace indicating if
       linalg_factorise has been called before.

    int linalg_solve(MultiMatrix *B, 
                     MultiMatrix *X,
                     void **linalg_workspace,
		     ErrorMsg error_message)
    -- Solve A X^T = B^T. X and B are dense, and we using row-major
       storage which explains the transpositions, ^T.

    The names of these 4 functions could in principle be anything, but
    I suggest to use the names above and append _"NAME_OF_WRAPPER".
*/
int linalg_initialise_dense_NR(MultiMatrix *A, 
			       void **linalg_workspace,
			       ErrorMsg error_message){
  DNR_structure *workspace;
  void *data;
  DataType Dtype;
  int ncol, nrow;
  ncol = A->ncol; nrow = A->nrow; Dtype = A->Dtype;
  //Test input:
  lasagna_test(A->ncol == A->nrow, error_message);
  lasagna_test((A->Dtype!=L_DBL)&&(A->Dtype!=L_DBL_CX), error_message);
  lasagna_test(A->Stype!=L_DNR, error_message);
  
  //Allocate stuff:
  lasagna_alloc(workspace,sizeof(DNR_structure),error_message);
  lasagna_alloc(workspace->luidx,sizeof(int)*(nrow+1),error_message);
  lasagna_alloc(data,GetByteSize(Dtype)*(ncol*nrow+1),error_message);
  lasagna_alloc(workspace->LUw,GetByteSize(Dtype)*(nrow+1),error_message);
  
  lasagna_call(CreateMatrix_DNR(workspace->LU,
				Dtype, 
				nrow,
				ncol,
				data,
				error_message),
	       error_message,error_message);
  workspace->neq = nrow;
  *linalg_workspace = (void *) workspace;

  return _SUCCESS_;
}

int linalg_finalise_dense_NR(void *linalg_workspace,
			     ErrorMsg error_message){
  DNR_structure *ws;
  ws = (DNR_structure *)linalg_workspace;
  free(ws->LU->data);
  //Data freed, now destroy multimatrix LU:
  DestroyMultiMatrix(ws->LU);
  free(ws->LUw);
  free(ws->luidx);
  free(ws);
  return _SUCCESS_;
}

