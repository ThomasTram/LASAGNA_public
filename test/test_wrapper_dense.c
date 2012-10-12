#include "common.h"
#include "linalg_wrapper_dense_NR.h"

int main(){
  MultiMatrix A,B,X;
  void *linalg_workspace;
  int i,j, n=479, nsys = 2;
  double *Adata, *Bdata, *Xdata;
  double complex *Axzdata, *Bzdata, *Xzdata;
  void *Adat, *Bdat, *Xdat;
  ErrorMsg error_message;
  DNRformat *StoreA;
  DNRformat *StoreB;
  DNRformat *StoreX;
  int funcreturn, readcount;
  FILE *myfile;
  EvolverOptions options;
  DataType Dtype;
  int loop, maxloop=10;
  
  for (loop=0; loop<maxloop; loop++){
 

    Adata = malloc(sizeof(double)*(n*n+1));
    Bdata = malloc(sizeof(double)*(n*nsys+1));
    Xdata = malloc(sizeof(double)*(n*nsys+1));
    Axzdata = malloc(sizeof(double complex)*(n*n+1));
    Bzdata = malloc(sizeof(double complex)*(n*nsys+1));
    Xzdata = malloc(sizeof(double complex)*(n*nsys+1));

    switch(loop%2){
    case (0):
      printf("Testing real double routine:\n");
      Dtype = L_DBL;
      Adat = Adata;
      Bdat = Bdata;
      Xdat = Xdata;
      break;
    case(1):
      printf("Testing complex double routine:\n");
      Dtype = L_DBL_CX;
      Adat = Axzdata;
      Bdat = Bzdata;
      Xdat = Xzdata;
      break;
    }
  

    CreateMatrix_DNR(&A,
		     Dtype, 
		     n,
		     n,
		     Adat,
		     error_message);
    CreateMatrix_DNR(&B,
		     Dtype, 
		     nsys,
		     n,
		     Bdat,
		     error_message);
    CreateMatrix_DNR(&X,
		     Dtype, 
		     nsys,
		     n,
		     Xdat,
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
    myfile = fopen("testlin_Adense.dat","r");
    readcount=fread(Adata,sizeof(double),n*n+1,myfile);
    printf("Read %d entries.\n",readcount);
    fclose(myfile);
    
    switch (Dtype){
    case (L_DBL):
      myfile = fopen("testlin_bdense.dat","r");
      readcount=fread(Bdata,sizeof(double),n*nsys+1,myfile);
      printf("Read %d entries.\n",readcount);
      fclose(myfile);
      break;
    case (L_DBL_CX):
      for (i=0; i<(n*n+1); i++)
	Axzdata[i] = Adata[i];
      myfile = fopen("testlin_Adense_imag.dat","r");
      readcount=fread(Adata,sizeof(double),n*n+1,myfile);
      printf("Read %d entries.\n",readcount);
      fclose(myfile);
      for (i=0; i<(n*n+1); i++)
	Axzdata[i] += I*Adata[i];
      myfile = fopen("testlin_bzrdense.dat","r");
      readcount=fread(Bdata,sizeof(double),n*nsys+1,myfile);
      printf("Read %d entries.\n",readcount);
      fclose(myfile);
      for (i=0; i<(n*nsys+1); i++)
	Bzdata[i] = Bdata[i];
      myfile = fopen("testlin_bzidense.dat","r");
      readcount=fread(Bdata,sizeof(double),n*nsys+1,myfile);
      printf("Read %d entries.\n",readcount);
      fclose(myfile);
      for (i=0; i<(n*nsys+1); i++)
	Bzdata[i] += I*Bdata[i];
      break;
    }
    //PrintMultiMatrix(&A,"Matrix A before factorisation.");
    PrintMultiMatrix(&B,"Right hand side B before factorisation.");
  
    funcreturn = linalg_factorise_dense_NR(linalg_workspace,
					   _TRUE_;
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
    free(Xzdata);
    free(Bzdata);
    free(Axzdata);
  
  }
  return 0;
}
