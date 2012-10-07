#include "multimatrix.h"
#include "complex.h"

int CreateMatrix_SCC(MultiMatrix *A,
		     StoreType Stype, 
		     int nrow,
		     int ncol,
		     int nnz;
		     int Ai*,
		     int Ap*,
		     void Ax*,
		     ErrorMsg error_message){
  SCCformat Store*;
  /** Note: The Ai, Ap and Ax arrays are assumed to be 
      allocated, but not neccessarily initialised. Thus,
      we can not use nnz = Ap[ncol]. */
  //Test inputs:
  lasagna_test(nrow<1, error_message);
  lasagna_test(ncol<1, error_message);
  lasagna_test((Dtype!=L_DBL)&&(Dtype!=L_DBL_CX), error_message);

  //Set values:
  A->Stype = L_SCC;
  A->Dtype = Dtype;
  A->nrow = nrow;
  A->ncol = ncol;
  
  //Allocate storage structure:
  lasagna_alloc(A->Store,sizeof(SCCformat),error_message);
  Store = (SCCformat *) A->Store;
  Store->nnz = nnz;
  Store->Ai = Ai;
  Store->Ap = Ap;
  Store->Ax = Ax;

  return _SUCCESS_;
}

int CreateMatrix_DNR(MultiMatrix *A,
		     DataType Dtype, 
		     int nrow,
		     int ncol,
		     void data*,
		     ErrorMsg error_message){
  DNRformat Store*;
  double *x_dbl,**p2p_dbl;
  double complex *x_dbl_cx,**p2p_dbl_cx;
  
  /** Important: data is assumed to be allocated and of size nrow*ncol+1.
      The matrix is stored in row-major order from Ax+1. Given a matrix B
      as ordinary, 0-indexed row-major, one can pass B-1 for Ax. This is not
      in compliance with standards but should work nonetheless.*/
  //Test inputs:
  lasagna_test(nrow<1, error_message);
  lasagna_test(ncol<1, error_message);
  lasagna_test((Dtype!=L_DBL)&&(Dtype!=L_DBL_CX), error_message);
  //Set values:
  A->Stype = DNR;
  A->Dtype = Dtype;
  A->nrow = nrow;
  A->ncol = ncol;
  
  //Allocate storage structure:
  lasagna_alloc(A->Store,sizeof(DNRformat),error_message);
  Store = (DNRformat *) A->Store;
  Store->Data = data;

  switch (A->Dtype){
  case (L_DBL):
    lasagna_alloc(Store->pointer_to_pointers,
		  sizeof(double *)*(nrow+1),
		  error_message);
    p2p_dbl = (double **) Store->Matrix;
    p2p_dbl[0] = NULL; //For definiteness.
    p2p_dbl[1] = (double *) Store->Data;
    for (i=2; i<=nrow; i++)
      p2p_dbl[i] = p2p_dbl[i-1]+ncol;
    break;
  case (L_DBL_CX):
    lasagna_alloc(Store->pointer_to_pointers,
		  sizeof(double complex *)*(nrow+1),
		  error_message);
    p2p_dbl_cx = (double complex **) Store->Matrix;
    p2p_dbl_cx[0] = NULL; //For definiteness.
    p2p_dbl_cx[1] = (double complex *) Store->Data;
    for (i=2; i<=nrow; i++)
      p2p_dbl_cx[i] = p2p_dbl_cx[i-1]+ncol;
    break;
  }
  return _SUCCESS_;
}

int DestroyMultiMatrix(MultiMatrix *A){
  /** Deallocates everything that was allocated by 
      CreateMatrix_xxx. If the pointers to the actual
      storage Ax, Ai, Ap, data,... is not available,
      they should be freed before to avoid memory leak. */
  int i;
  DNRformat *StoreDNR;
  if (A->Stype == L_DNR){
    StoreDNR = (DNRformat *) A->Store;
    free(StoreDNR->Matrix);
  }
  free(A->Store);
  return _SUCCESS_;
}

size_t GetByteSize(DataType Dtype){
  switch (A->Dtype){
  case (L_DBL):
    return sizeof(double);
    break;
  case (L_DBL_CX):
    return sizeof(double complex);
    break;
  }
  return -1;
}
