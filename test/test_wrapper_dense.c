#include "common.h"
#include "linalg_wrapper_dense_NR.h"

int main(){
  MultiMatrix A,B,X;
  void *linalg_workspace;
  int i,j, n=479, nsys = 2;
  double *Adata, *Bdata, *Xdata;
  ErrorMsg error_message;
  DNRformat *StoreA;
  DNRformat *StoreB;
  DNRformat *StoreX;
  int funcreturn, readcount;
  FILE *myfile;
  EvolverOptions options;
  Adata = malloc(sizeof(double)*(n*n+1));
  Bdata = malloc(sizeof(double)*(n*nsys+1));
  Xdata = malloc(sizeof(double)*(n*nsys+1));

  CreateMatrix_DNR(&A,
		   L_DBL, 
		   n,
		   n,
		   Adata,
		   error_message);
  CreateMatrix_DNR(&B,
		   L_DBL, 
		   nsys,
		   n,
		   Bdata,
		   error_message);
  CreateMatrix_DNR(&X,
		   L_DBL, 
		   nsys,
		   n,
		   Xdata,
		   error_message);

  //Initialise solver:
  funcreturn = linalg_initialise_dense_NR(&A,
					  &options,
					  &linalg_workspace,
					  error_message);
  if (funcreturn==_FAILURE_){
    printf("Error: %s\n",error_message);
    return _FAILURE_;
  }

  /** Solve equation */
  //StoreA = (DNRformat *) A.Store;
  //StoreB = (DNRformat *) B.Store;
  /**Setup test problem:
  Amat = (double **) StoreA->Matrix;
  StoreB = (DNRformat *)B.Store;
  Bmat = (double **) StoreB->Matrix;
  StoreX = (DNRformat *)X.Store;
  Xmat = (double **) StoreX->Matrix;
  for (j=1; j<=n; j++){
    for (i=1; i<=n; i++){
      Amat[j][i] = sqrt(i+j*i)-i+2.0*j*j;
    }
    Bmat[1][j] = j;
  }
  */
  myfile = fopen("testlin_Adense.dat","r");
  readcount=fread(Adata,sizeof(double),n*n+1,myfile);
  printf("Read %d entries.\n",readcount);
  fclose(myfile);
  
  myfile = fopen("testlin_bdense.dat","r");
  readcount=fread(Bdata,sizeof(double),n*nsys+1,myfile);
  printf("Read %d entries.\n",readcount);
  fclose(myfile);
  

  PrintMultiMatrix(&A,"Matrix A before factorisation.");
  PrintMultiMatrix(&B,"Right hand side B before factorisation.");
  
  funcreturn = linalg_factorise_dense_NR(&A, 
					 linalg_workspace,
					 error_message);
  if (funcreturn==_FAILURE_){
    printf("Error: %s\n",error_message);
    return _FAILURE_;
  }
  funcreturn = linalg_solve_dense_NR(&B, 
				     &X,
				     linalg_workspace,
				     error_message);
  if (funcreturn==_FAILURE_){
    printf("Error: %s\n",error_message);
    return _FAILURE_;
  }
  
  PrintMultiMatrix(&X,"Solution vector:");


  funcreturn = linalg_finalise_dense_NR(linalg_workspace,
					error_message);
  if (funcreturn==_FAILURE_){
    printf("Error: %s\n",error_message);
    return _FAILURE_;
  }
  

  DestroyMultiMatrix(&A);
  DestroyMultiMatrix(&B);
  DestroyMultiMatrix(&X);
  free(Xdata);
  free(Bdata);
  free(Adata);


  return 0;
}
