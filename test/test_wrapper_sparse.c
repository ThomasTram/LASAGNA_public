#include "multimatrix.h"
#include "linalg_wrapper_SuperLU.h"
#include "linalg_wrapper_sparse.h"

int main(){
  int nnz,readcount;
  int cols=479;
  int nsys=2;
  int *Ai;
  int *Ap;
  double *Ax, *b, *x, *Axim;
  double complex *Axz, *bz, *xz;
  void *Adat, *Xdat, *Bdat;
  FILE *myfile;
  MultiMatrix A, B, X;
  EvolverOptions options;
  ErrorMsg error_message;
  int funcreturn;
  DataType Dtype;
  void *linalg_workspace;
  int i, loop, maxloop=4;
  //Function pointers:
  int (*linalg_initialise)(MultiMatrix *, EvolverOptions *, void **, ErrorMsg);
  int (*linalg_finalise)(void *, ErrorMsg);
  int (*linalg_factorise)(void *, ErrorMsg);
  int (*linalg_solve)(MultiMatrix *, MultiMatrix *, void *, ErrorMsg);

  DefaultEvolverOptions(&options);
  options.Cores=4;

  for (loop=0; loop<maxloop; loop++){
    Dtype = L_DBL_CX;
    switch (loop%4){
    case (0):
      Dtype = L_DBL;
    case (1):
      fprintf(stdout,"Testing SuperLU wrapper:\n");
      linalg_initialise = &(linalg_initialise_SuperLU);
      linalg_finalise = &(linalg_finalise_SuperLU);
      linalg_factorise = &(linalg_factorise_SuperLU);
      linalg_solve = &(linalg_solve_SuperLU);
      break;      
    case (2):
      Dtype = L_DBL;
    case (3):
      fprintf(stdout,"Testing sparse wrapper:\n");
      linalg_initialise = &(linalg_initialise_sparse);
      linalg_finalise = &(linalg_finalise_sparse);
      linalg_factorise = &(linalg_factorise_sparse);
      linalg_solve = &(linalg_solve_sparse);
      break;      
    }

    switch(Dtype){
    case (L_DBL):
      printf("Testing double routines:\n");
      break;
    case (L_DBL_CX):
      printf("Testing double complex routines:\n");
      break;
    }
    //Construct A:
    //------------------------------------
    Ap = malloc(sizeof(int)*(cols+1));
    myfile = fopen("testlin_Ap.dat","r");
    readcount=fread(Ap,sizeof(int),(cols+1),myfile);
    printf("Read %d entries.\n",readcount);
    fclose(myfile);
    nnz = Ap[cols];


    printf("nnz=%d\n",nnz);
    Ai = malloc(sizeof(int)*nnz);
    myfile = fopen("testlin_Ai.dat","r");
    readcount=fread(Ai,sizeof(int),nnz,myfile);
    printf("Read %d entries.\n",readcount);
    fclose(myfile);

    /** Debug region 
    int *Cp, *Ci, i,j;
    FILE *fI, *fJ;
    get_pattern_A_plus_AT(Ap, Ai, cols, &Cp, &Ci, error_message);
    fI = fopen("I.dat","w");
    fJ = fopen("J.dat","w");

    for (j=0; j<cols; j++){
      for (i=Cp[j]; i<Cp[j+1]; i++){
	fprintf(fI,"%d ",Ci[i]+1);
	fprintf(fJ, "%d ",j+1);
      }
    }
    fclose(fI);
    fclose(fJ);
    return 0;
    */

    Ax = malloc(sizeof(double)*nnz);
    Axim = malloc(sizeof(double)*nnz);
    Axz = malloc(sizeof(double complex)*nnz);
    b = malloc(sizeof(double)*(cols*nsys+1));
    x = malloc(sizeof(double)*(cols*nsys+1));
    bz = malloc(sizeof(double complex)*(cols*nsys+1));
    xz = malloc(sizeof(double complex)*(cols*nsys+1));
   
    switch (Dtype){
    case (L_DBL):
      Adat = (void*) Ax;
      Bdat = (void*) b;
      Xdat = (void*) x;
      break;
    case (L_DBL_CX):
      Adat = (void*) Axz;
      Bdat = (void*) bz;
      Xdat = (void*) xz;
      break;
    }
    
    CreateMatrix_DNR(&(B),
		     Dtype, 
		     nsys,
		     cols,
		     Bdat,
		     error_message);
    CreateMatrix_DNR(&(X),
		     Dtype, 
		     nsys,
		     cols,
		     Xdat,
		     error_message);
    CreateMatrix_SCC(&(A),
		     Dtype,
		     cols,
		     cols,
		     nnz,
		     Ai,
		     Ap,
		     Adat,
		     error_message);
    
    //Initialise solver:
    funcreturn = linalg_initialise(&A,
				   &options,
				   &linalg_workspace,
				   error_message);

    myfile = fopen("testlin_Ax.dat","r");
    readcount=fread(Ax,sizeof(double),nnz,myfile);
    printf("Read %d entries.\n",readcount);
    fclose(myfile);
    if (Dtype==L_DBL_CX){
      for (i=0; i<nnz; i++)
	Axz[i] = Ax[i];
      myfile = fopen("testlin_Azim.dat","r");
      readcount=fread(Ax,sizeof(double),nnz,myfile);
      printf("Read %d entries.\n",readcount);
      fclose(myfile);
      for (i=0; i<nnz; i++)
	Axz[i] += I*Ax[i];
    }


    switch (Dtype){
    case (L_DBL):
      myfile = fopen("testlin_bdense.dat","r");
      readcount=fread(b,sizeof(double),(cols*nsys+1),myfile);
      printf("Read %d entries from testlin_bdense.\n",readcount);
      fclose(myfile);
      break;
    case (L_DBL_CX):
      myfile = fopen("testlin_bzrdense.dat","r");
      readcount=fread(b,sizeof(double),(cols*nsys+1),myfile);
      printf("Read %d entries from testlin_bdense.\n",readcount);
      fclose(myfile);
      for (i=0; i<(cols*nsys+1); i++)
	bz[i] = b[i];
      myfile = fopen("testlin_bzidense.dat","r");
      readcount=fread(b,sizeof(double),(cols*nsys+1),myfile);
      printf("Read %d entries from testlin_bdense.\n",readcount);
      fclose(myfile);
      for (i=0; i<(cols*nsys+1); i++)
	bz[i] += I*b[i];
      break;
    }

    funcreturn = linalg_factorise(linalg_workspace,
				  error_message);
    if (funcreturn==_FAILURE_){
      printf("Error: %s\n",error_message);
      return _FAILURE_;
    }
    funcreturn = linalg_solve(&B, 
			      &X,
			      linalg_workspace,
			      error_message);
    if (funcreturn==_FAILURE_){
      printf("Error: %s\n",error_message);
      return _FAILURE_;
    }
  
    PrintMultiMatrix(&X,"Solution vector:");


    funcreturn = linalg_finalise(linalg_workspace,
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
    free(Axz);
    free(bz);
    free(xz);
  }
  return 0;
}
