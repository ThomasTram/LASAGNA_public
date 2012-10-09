#include "multimatrix.h"
#include "linalg_wrapper_SuperLU.h"

int main(){
  int nnz,readcount;
  int cols=479;
  int nsys=2;
  int *Ai;
  int *Ap;
  double *Ax, *b, *x;
  FILE *myfile;
  MultiMatrix A, B, X;
  EvolverOptions options= {.Cores=4, .flag=1};
  ErrorMsg error_message;
  int funcreturn;
  void *linalg_workspace;

  //Construct A:
  //------------------------------------
  Ap = malloc(sizeof(int)*(cols+1));
  myfile = fopen("testlin_Ap.dat","r");
  readcount=fread(Ap,sizeof(int),(cols+1),myfile);
  printf("Read %d entries.\n",readcount);
  fclose(myfile);
  nnz = Ap[cols];

  Ai = malloc(sizeof(int)*nnz);
  myfile = fopen("testlin_Ai.dat","r");
  readcount=fread(Ai,sizeof(int),nnz,myfile);
  printf("Read %d entries.\n",readcount);
  fclose(myfile);

  Ax = malloc(sizeof(double)*nnz);
  b = malloc(sizeof(double)*(cols*nsys+1));
  x = malloc(sizeof(double)*(cols*nsys+1));
  CreateMatrix_DNR(&(B),
		   L_DBL, 
		   nsys,
		   cols,
		   b,
		   error_message);
  CreateMatrix_DNR(&(X),
		   L_DBL, 
		   nsys,
		   cols,
		   x,
		   error_message);

  CreateMatrix_SCC(&(A),
		   L_DBL,
		   cols,
		   cols,
		   nnz,
		   Ai,
		   Ap,
		   Ax,
		   error_message);

  
    //Initialise solver:
  funcreturn = linalg_initialise_SuperLU(&A,
					 &options,
					 &linalg_workspace,
					 error_message);

  myfile = fopen("testlin_Ax.dat","r");
  readcount=fread(Ax,sizeof(double),nnz,myfile);
  printf("Read %d entries.\n",readcount);
  fclose(myfile);

  myfile = fopen("testlin_bdense.dat","r");
  readcount=fread(b,sizeof(double),(cols*nsys+1),myfile);
  printf("Read %d entries from testlin_bdense.\n",readcount);
  fclose(myfile);

  funcreturn = linalg_factorise_SuperLU(&A, 
					linalg_workspace,
					error_message);
  if (funcreturn==_FAILURE_){
    printf("Error: %s\n",error_message);
    return _FAILURE_;
  }
  funcreturn = linalg_solve_SuperLU(&B, 
				    &X,
				    linalg_workspace,
				    error_message);
  if (funcreturn==_FAILURE_){
    printf("Error: %s\n",error_message);
    return _FAILURE_;
  }
  
  PrintMultiMatrix(&X,"Solution vector:");


  funcreturn = linalg_finalise_SuperLU(linalg_workspace,
					error_message);
  if (funcreturn==_FAILURE_){
    printf("Error: %s\n",error_message);
    return _FAILURE_;
  }
  

  DestroyMultiMatrix(&A);
  DestroyMultiMatrix(&B);
  DestroyMultiMatrix(&X);
  free(Ax);
  free(b);
  free(x);
  free(Ap);
  free(Ai);

  return 0;
}
