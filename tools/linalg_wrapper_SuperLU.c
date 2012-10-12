#include "linalg_wrapper_SuperLU.h"
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

    int linalg_factorise(void *linalg_workspace,
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
int linalg_initialise_SuperLU(MultiMatrix *A, 
			      EvolverOptions *options,
			      void **linalg_workspace,
			      ErrorMsg error_message){
  SCCformat *Store=A->Store;
  SLU_structure *workspace;
  void *data;
  DataType Dtype;
  int ncol, nrow, j;
  int *Aw;
  int Alen;
  //SuperLU parameters:
  int panel_size, relax;
  trans_t trans;
  fact_t   fact;
  yes_no_t refact, usepr;
  double diag_pivot_thresh, drop_tol;
  void *work;
  int lwork;


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
  
  //Allocate stuff:
  lasagna_alloc(workspace,sizeof(SLU_structure),error_message);
  lasagna_alloc(workspace->perm_c,sizeof(int)*(nrow+1),error_message);
  lasagna_alloc(workspace->perm_r,sizeof(int)*(nrow),error_message);
  switch (A->Dtype){
  case (L_DBL):
    dCreate_CompCol_Matrix(&(workspace->A),
			   ncol,nrow,Store->nnz,
			   Store->Ax,
			   Store->Ai,
			   Store->Ap,
			   SLU_NC,
			   SLU_D,
			   SLU_GE);
    dCreate_Dense_Matrix(&(workspace->X), nrow, 1, NULL, nrow,
			 SLU_DN, SLU_D, SLU_GE);
    break;
  case (L_DBL_CX):
    zCreate_CompCol_Matrix(&(workspace->A),
			   ncol,nrow,Store->nnz,
			   Store->Ax,
			   Store->Ai,
			   Store->Ap,
			   SLU_NC,
			   SLU_Z,
			   SLU_GE);
    zCreate_Dense_Matrix(&(workspace->X), nrow, 1, NULL, nrow,
			 SLU_DN, SLU_Z, SLU_GE);
  }
  /**Handle case when no factorisations are done, 
     so L and U has not been allocated :
  */
  workspace->L.Store = NULL;
  workspace->U.Store = NULL;

  Alen = colamd_recommended (Store->nnz, ncol, ncol);
  lasagna_alloc(Aw,sizeof(int)*Alen,error_message);
  for (j=0; j<Store->nnz; j++)
    Aw[j] = Store->Ai[j];
  for (j=0; j<=ncol; j++)
    workspace->perm_c[j] = Store->Ap[j];

  lasagna_test(colamd (ncol, ncol, Alen, Aw, workspace->perm_c, NULL) != TRUE,
	       error_message,
	       "colamd failed!");
  free(Aw);

  /** Set SuperLU parameters: 
      ------------------------*/
  fact               = DOFACT;
  refact             = NO;
  trans              = NOTRANS;
  panel_size         = sp_ienv(1);
  relax              = sp_ienv(2);
  diag_pivot_thresh  = 1.0;
  usepr              = NO;
  drop_tol           = 0.0;
  work               = NULL;
  lwork              = 0;
  /** --------------------- */

  StatAlloc(ncol, options->Cores, panel_size, relax, &workspace->Gstat);
  StatInit(ncol, options->Cores, &workspace->Gstat);
  switch (A->Dtype){
  case (L_DBL):
    pdgstrf_init(options->Cores, fact, trans, refact, panel_size, relax,
		 diag_pivot_thresh, usepr, drop_tol, 
		 workspace->perm_c, workspace->perm_r,
		 work, lwork, &(workspace->A), &(workspace->AC), 
		 &(workspace->superlumt_options), &(workspace->Gstat));
    break;
  case (L_DBL_CX):
    pzgstrf_init(options->Cores, fact, trans, refact, panel_size, relax,
		 diag_pivot_thresh, usepr, drop_tol, 
		 workspace->perm_c, workspace->perm_r,
		 work, lwork, &(workspace->A), &(workspace->AC), 
		 &(workspace->superlumt_options), &(workspace->Gstat));
    break;
  }
  
  workspace->neq = nrow;
  *linalg_workspace = (void *) workspace;
  
  return _SUCCESS_;
}

int linalg_finalise_SuperLU(void *linalg_workspace,
			    ErrorMsg error_message){
  
  SLU_structure *ws=linalg_workspace;  
  pxgstrf_finalize(&(ws->superlumt_options), &(ws->AC));
  StatFree(&(ws->Gstat));
  
  if (ws->L.Store != NULL) Destroy_SuperNode_SCP(&(ws->L));
  if (ws->U.Store != NULL) Destroy_CompCol_NCP(&(ws->U));
  
  Destroy_SuperMatrix_Store(&(ws->A)); 
  Destroy_SuperMatrix_Store(&(ws->X));
  free(ws->perm_c);
  free(ws->perm_r);
  free(ws);
  return _SUCCESS_;
}


int linalg_factorise_SuperLU(void *linalg_workspace,
			     int has_changed_significantly,
			     ErrorMsg error_message){
  SLU_structure *ws= linalg_workspace;
  /** The actual data in ws->A is the same as in the MultiMatrix
      A which was passed to linalg_initialise_SuperLU. Thus, since
      AC has the same data as in A, we can just call the factorisation
      routine now, given that A has been modified accordingly outside.
  */

  switch(ws->A.Dtype){
  case (SLU_D):
    pdgstrf(&(ws->superlumt_options), 
	    &(ws->AC), 
	    ws->perm_r, 
	    &(ws->L), 
	    &(ws->U), 
	    &(ws->Gstat), 
	    &(ws->SLU_info));
    break;
  case (SLU_Z):
    pzgstrf(&(ws->superlumt_options), 
	    &(ws->AC), 
	    ws->perm_r, 
	    &(ws->L), 
	    &(ws->U), 
	    &(ws->Gstat), 
	    &(ws->SLU_info));
    break;
  }
  //Update superlu_mt_options structure:
  ws->superlumt_options.fact = FACTORED;
  ws->superlumt_options.refact = YES;
  return _SUCCESS_;
}

int linalg_solve_SuperLU(MultiMatrix *B, 
			 MultiMatrix *X,
			 void *linalg_workspace,
			 ErrorMsg error_message){
  SLU_structure *ws= linalg_workspace;
  DNRformat *StoreB=B->Store;
  DNRformat *StoreX=X->Store;
  int i;
  DNformat *ws_StoreX = ws->X.Store;



  /** Copy data from MultiMatrix B into MultiMatrix X. 
      This can be done by memcpy. */
  memcpy(StoreX->Data, 
	 StoreB->Data, 
	 GetByteSize(B->Dtype)*((B->nrow)*(B->ncol)+1));

  /** Copy information from MultiMatrix X into SuperMatrix X:
      Remember that it is the rows in the input matrix B we are
      solving for. (Because of the DNR format scheme).
      We do the transposition implicitly by assigning lda=ncol. */
  ws_StoreX->lda = X->ncol;
  ws->X.ncol = X->nrow;
  ws->X.nrow = X->ncol;
  //Solve system:
  switch(B->Dtype){
  case (L_DBL):
    //Compensate for indexing scheme:
    ws_StoreX->nzval = (void *) (( (double*) StoreX->Data) + 1);
    dgstrs(NOTRANS, &(ws->L),&(ws->U), ws->perm_r, ws->perm_c, 
	   &(ws->X), &(ws->Gstat), &(ws->SLU_info));
    break;
  case (L_DBL_CX):
    //Compensate for indexing scheme:
    ws_StoreX->nzval = (void *) (( (double complex*) StoreX->Data) + 1);
    /** The elements in StoreX->Data is of type double complex as defined
	in complex.h, but SuperLU has its own complex numbers defined as
	typedef struct{double r; double i} doublecomplex. However, the two
	are actually equivalent. */
    zgstrs(NOTRANS, &(ws->L),&(ws->U), ws->perm_r, ws->perm_c, 
	   &(ws->X), &(ws->Gstat), &(ws->SLU_info));
    break;
  }
  
  return _SUCCESS_;
}
