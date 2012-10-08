#include "linalg_wrapper_dense_NR.h"
#include "complex.h"
#define TINY 1e-50
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
int linalg_initialise_dense_NR(MultiMatrix *A, 
			       void **linalg_workspace,
			       ErrorMsg error_message){
  DNR_structure *workspace;
  void *data;
  DataType Dtype;
  int ncol, nrow;
  ncol = A->ncol; nrow = A->nrow; Dtype = A->Dtype;
  //Test input:
  lasagna_test(A->ncol == A->nrow, 
	       error_message, 
	       "Matrix not square!");
  lasagna_test((A->Dtype!=L_DBL)&&(A->Dtype!=L_DBL_CX), 
	       error_message,
	       "Unknown datatype in A.");
  lasagna_test(A->Stype!=L_DNR, 
	       error_message,
	       "This wrapper only supports dense input matrix.");
  
  //Allocate stuff:
  lasagna_alloc(workspace,sizeof(DNR_structure),error_message);
  lasagna_alloc(workspace->luidx,sizeof(int)*(nrow+1),error_message);
  lasagna_alloc(data,GetByteSize(Dtype)*(ncol*nrow+1),error_message);
  lasagna_alloc(workspace->LUw,sizeof(double)*(nrow+1),error_message);
  
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
  DNRformat *Store;
  ws = (DNR_structure *)linalg_workspace;
  Store = (DNRformat *) ws->LU->Store;
  free(Store->Data);
  //Data freed, now destroy multimatrix LU:
  DestroyMultiMatrix(ws->LU);
  free(ws->LUw);
  free(ws->luidx);
  free(ws);
  return _SUCCESS_;
}


int linalg_factorise_dense_NR(MultiMatrix *A, 
			      void *linalg_workspace,
			      ErrorMsg error_message){
  DNR_structure *ws= linalg_workspace;
  DNRformat *StoreA=A->Store;
  DNRformat *StoreLU=ws->LU->Store;
  double luparity;
  //Copy data from A into ws->LU. This can be done by memcpy.
  memcpy(StoreLU->Data, StoreA->Data, GetByteSize(ws->LU->Dtype)*(ws->neq+1));

  switch(A->Dtype){
  case (L_DBL):
    ludcmp((double **) StoreLU->Matrix, 
	   ws->neq, 
	   ws->luidx, 
	   &luparity, 
	   ws->LUw);
    break;
  case (L_DBL_CX):
    ludcmp_cx((double complex **) StoreLU->Matrix,
	      ws->neq, 
	      ws->luidx, 
	      &luparity, 
	      ws->LUw);
    break;
  }
  return _SUCCESS_;
}


int linalg_solve_dense_NR(MultiMatrix *B, 
			  MultiMatrix *X,
			  void *linalg_workspace,
			  ErrorMsg error_message){
  DNR_structure *ws= linalg_workspace;
  DNRformat *StoreLU=ws->LU->Store;
  DNRformat *StoreB=B->Store;
  DNRformat *StoreX=X->Store;
  int i;
  double **xmat;
  double complex **xmat_cx;

  //Copy data from B into X. This can be done by memcpy.
  memcpy(StoreX->Data, 
	 StoreB->Data, 
	 GetByteSize(B->Dtype)*((B->nrow)*(B->ncol)+1));

  //Loop over right hand sides:
  for(i=1; i<=B->nrow; i++){ 
    switch(B->Dtype){
    case (L_DBL):
      xmat = (double **) StoreX->Matrix;
      lubksb((double **) StoreLU->Matrix,
	     B->ncol,
	     ws->luidx,
	     xmat[i]);
      break;
    case (L_DBL_CX):
      xmat_cx = (double complex **) StoreX->Matrix;;
      lubksb_cx((double complex **) StoreLU->Matrix,
		B->ncol,
		ws->luidx,
		xmat_cx[i]);
      break;
    }
  }
  return _SUCCESS_;
}


















int ludcmp(double **a, int n, int *indx, double *d, double *vv){
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++) {
      if ((temp=fabs(a[i][j])) > big) big=temp;
    }
    if (big == 0.0) return _FAILURE_;
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  return _SUCCESS_;
}

int ludcmp_cx(double complex **a, int n, int *indx, double *d, double *vv){
  int i,imax=0,j,k;
  double big,dum,temp;
  double complex sum,dum2;
  *d=1.0;
  //Loop over rows to find largest element of each column
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++) {
      if ((temp=cabs(a[i][j])) > big) big=temp;
    }
    if (big == 0.0) return _FAILURE_;
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*cabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum2=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum2;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (cabs(a[j][j]) == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum2=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum2;
    }
  }
  return _SUCCESS_;
}


int lubksb(double **a, int n, int *indx, double b[]){
  int i,ii=0,ip,j;
  double sum;
  for (i=1;i<=n;i++) {
    ip=indx[i];
    if((ip<1)||(ip>n)) printf("Error in LU backsubstitution. (index is %d, n=%d)\n",ip,n);
    sum=b[ip];
    b[ip]=b[i];
    if (ii) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
  return _SUCCESS_;
}


int lubksb_cx(double complex **a, int n, int *indx, double complex b[]){
  int i,ii=0,ip,j;
  double complex sum;
  for (i=1;i<=n;i++) {
    ip=indx[i];
    if((ip<1)||(ip>n)) printf("WTFWTFWTF!!! i=%d, n=%d.\n",ip,n);
    sum=b[ip];
    b[ip]=b[i];
    if (ii) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (cabs(sum)!=0.0) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
  return _SUCCESS_;
}
