#include "evolver_common.h"

int DefaultEvolverOptions(EvolverOptions *opt, LinAlgWrapper linalg){
  int i;
  extern int linalg_initialise_dense_NR();
  extern int linalg_finalise_dense_NR();
  extern int linalg_factorise_dense_NR();
  extern int linalg_solve_dense_NR();

  extern int linalg_initialise_sparse();
  extern int linalg_finalise_sparse();
  extern int linalg_factorise_sparse();
  extern int linalg_solve_sparse();

  extern int linalg_initialise_SuperLU();
  extern int linalg_finalise_SuperLU();
  extern int linalg_factorise_SuperLU();
  extern int linalg_solve_SuperLU();

  opt->AbsTol=1e-6;
  opt->RelTol=1e-3;  
  opt->used_in_output=NULL;
  opt->t_vec=NULL;
  opt->tres=1;
  opt->Cores=1;
  opt->EvolverVerbose=1; 
  opt->output=NULL;
  opt->print_variables=NULL;
  opt->stop_function=NULL;
  for (i=0; i<10; i++){
    opt->Stats[i]= 0;
    opt->Flags[i]= 0;
  }
  opt->LinAlgVerbose=1;
  opt->Ap = NULL;
  opt->Ai = NULL;
  opt->J_pointer_flag = _FALSE_;
  switch (linalg){
  case (LINALG_WRAPPER_DENSE_NR):
    opt->linalg_initialise=linalg_initialise_dense_NR;
    opt->linalg_finalise=linalg_finalise_dense_NR;
    opt->linalg_factorise=linalg_factorise_dense_NR;
    opt->linalg_solve=linalg_solve_dense_NR;
    opt->use_sparse = _FALSE_;
    break;
  case (LINALG_WRAPPER_SPARSE):
    opt->linalg_initialise=linalg_initialise_sparse;
    opt->linalg_finalise=linalg_finalise_sparse;
    opt->linalg_factorise=linalg_factorise_sparse;
    opt->linalg_solve=linalg_solve_sparse;
    opt->use_sparse = _TRUE_;
    break;
  case (LINALG_WRAPPER_SUPERLU):
    opt->linalg_initialise=linalg_initialise_SuperLU;
    opt->linalg_finalise=linalg_finalise_SuperLU;
    opt->linalg_factorise=linalg_factorise_SuperLU;
    opt->linalg_solve=linalg_solve_SuperLU;
      opt->use_sparse = _TRUE_;
    break;
  default:
    opt->linalg_initialise=NULL;
    opt->linalg_finalise=NULL;
    opt->linalg_factorise=NULL;
    opt->linalg_solve=NULL;
    opt->use_sparse = _FALSE_;
    break;
  }

  return _SUCCESS_;
}


/**********************************************************************/
/* Here are some routines related to the calculation of the jacobian: */
/* "numjac", "initialize_jacobian", "uninitialize_jacobian",					*/
/* "initialize_numjac_workspace", "uninitialize_numjac_workspace".		*/
/**********************************************************************/
int numjac(int (*derivs)(double x, 
			 double * y,
			 double * dy, 
			 void * parameters_and_workspace, 
			 ErrorMsg error_message),
	   double t, 
	   double *y, 
	   double *fval,
	   MultiMatrix* J, 
	   void *numjac_workspace,
	   double thresh, 
	   int neq, 
	   int *nfe, 
	   void * parameters_and_workspace_for_derivs,
	   ErrorMsg error_message){
  /*	Routine that computes the jacobian numerically. It is based on the numjac 
	implementation in MATLAB, but a feature for recognising sparsity in the
	jacobian and taking advantage of that has been added.
  */
  struct numjac_workspace *nj_ws=numjac_workspace;
  double eps=DBL_EPSILON, br=pow(eps,0.875),bl=pow(eps,0.75),bu=pow(eps,0.25);
  double facmin=pow(eps,0.78),facmax=0.1;
  int logjpos, pattern_broken;
  double tmpfac,difmax2=0.,del2,ffscale;
  int i,j,rowmax2;
  double maxval1,maxval2;
  int colmax,group,row,nz,nz2;
  double Fdiff_absrm,Fdiff_new;
  double **dFdy,*fac, *Ax;
  int *Ap=NULL, *Ai=NULL;
  FILE *jacfile;
  SCCformat *StoreSCC;
  DNRformat *StoreDNR;

  fac = nj_ws->jacvec;
  for(j=1;j<=neq;j++){
    nj_ws->yscale[j] = max(fabs(y[j]),thresh);
    nj_ws->del[j] = (y[j] + fac[j] * nj_ws->yscale[j]) - y[j];
  }

  /*Select an increment del for a difference approximation to
    column j of dFdy.  The vector fac accounts for experience
    gained in previous calls to numjac.
  */
  for(j=1;j<=neq;j++){
    if (nj_ws->del[j]==0.0){
      for(;;){
	if (fac[j] < facmax){
	  fac[j] = min(100*fac[j],facmax);
	  nj_ws->del[j] = (y[j] + fac[j]*nj_ws->yscale[j]) - y[j];
	  if(nj_ws->del[j]==0.0){
	    break;
	  }
	}
	else{
	  nj_ws->del[j] = thresh;
	  break;
	}
      }
    }
  }

  /* keep del pointing into region: */
  for(j=1;j<=neq;j++){
    if (fval[j]>=0.0){
      nj_ws->del[j] = fabs(nj_ws->del[j]);
    }
    else{
      nj_ws->del[j] = -fabs(nj_ws->del[j]);
    }
  }

  switch(J->Stype){
  case(L_SCC):
    StoreSCC = J->Store;
    Ap = StoreSCC->Ap;
    Ai = StoreSCC->Ai;
    Ax = (double*) StoreSCC->Ax;
    //printf("Sparse calculation..neq=%d, has grouping=%d\n",neq,jac->has_grouping);
    //printf("jac->pattern_supplied = %d\n",jac->pattern_supplied);
    colmax = nj_ws->max_group+1;
    //printf("->groups=%d/%d.\n",colmax,neq);
    for(j=1;j<=colmax;j++){
      /*loop over groups */
      group = j-1;
      for(i=1;i<=neq;i++){
	/*Add y-vector.. */
	nj_ws->ydel_Fdel[i][j] = y[i];
	/*Add del of all groupmembers:*/
	if(nj_ws->col_group[i-1]==group) 
	  nj_ws->ydel_Fdel[i][j] +=nj_ws->del[i];
      }
    }
    break;
  case(L_DNR):
    //Ordinary dense method
    StoreDNR = J->Store;
    dFdy = (double **)StoreDNR->Matrix; /* Assign pointer to dfdy directly for easier notation. */
    //printf("Normal calculation...\n");
    /*Normal calculation: */
    colmax = neq;
    for(j=1;j<=neq;j++){
      for(i=1;i<=neq;i++){
	nj_ws->ydel_Fdel[i][j] = y[i];
      }
      nj_ws->ydel_Fdel[j][j] += nj_ws->del[j];
    }

    break;
  }

  /* The next section should work regardless of sparse...*/
  /* Evaluate the function at y+delta vectors:*/
  for(j=1;j<=colmax;j++){
    for(i=1;i<=neq;i++){
      nj_ws->yydel[i] = nj_ws->ydel_Fdel[i][j];
    }
    lasagna_call((*derivs)(t,
			   nj_ws->yydel+1,
			   nj_ws->ffdel+1,
			   parameters_and_workspace_for_derivs,
			   error_message),
		 error_message,error_message);

    *nfe+=1;
    for(i=1;i<=neq;i++) 
      nj_ws->ydel_Fdel[i][j] = nj_ws->ffdel[i];
  }

  switch(J->Stype){
  case(L_SCC):
    /*Using the Fdel array, form the jacobian and construct max-value arrays.
      First we do it for the sparse case, then for the normal case:*/
    /* Sparse case:*/
    for(j=0;j<neq;j++){
      /*Loop over columns, and assign corresponding group:*/
      group = nj_ws->col_group[j];
      Fdiff_new = 0.0;
      Fdiff_absrm = 0.0;
      
      for(i=Ap[j];i<Ap[j+1];i++){
	/* Loop over rows in the sparse matrix */
	row = Ai[i]+1;
	/* Do I want to construct the full jacobian? No, that is ugly..*/
	Fdiff_absrm = max(Fdiff_absrm,fabs(Fdiff_new));
	Fdiff_new = nj_ws->ydel_Fdel[row][group+1]-fval[row]; /*Remember to access the column of the corresponding group */
	if (fabs(Fdiff_new)>=Fdiff_absrm){
	  nj_ws->Rowmax[j+1] = row;
	  nj_ws->Difmax[j+1] = Fdiff_new;
	}
	/* Assign value to sparse rep of jacobian: */
	Ax[i] = Fdiff_new/nj_ws->del[j+1];
      }
      /* The maximum numerical value of Fdel in true column j+1*/
      nj_ws->absFdelRm[j+1] = fabs(nj_ws->ydel_Fdel[nj_ws->Rowmax[j+1]][group+1]);
    }
    break;
  case(L_DNR):
    /*Normal case:*/
    for(j=1;j<=neq;j++){
      Fdiff_new = 0.0;
      Fdiff_absrm = 0.0;
      for(i=1;i<=neq;i++){
	Fdiff_absrm = max(fabs(Fdiff_new),Fdiff_absrm);
	Fdiff_new = nj_ws->ydel_Fdel[i][j] - fval[i];
	dFdy[i][j] = Fdiff_new/nj_ws->del[j];
	/*Find row maximums:*/
	if(fabs(Fdiff_new)>=Fdiff_absrm){
	  /* Found new max location in column */
	  nj_ws->Rowmax[j] = i;
	  nj_ws->Difmax[j] = fabs(Fdiff_new);
	}
      }
      nj_ws->absFdelRm[j] = fabs(nj_ws->ydel_Fdel[nj_ws->Rowmax[j]][j]);
    }
  }
  
  /* Adjust fac for next call to numjac. */
  for(i=1;i<=neq;i++){
    nj_ws->absFvalue[i] = fabs(fval[i]);
  }
  for(j=1;j<=neq;j++){
    nj_ws->absFvalueRm[j] = nj_ws->absFvalue[nj_ws->Rowmax[j]];
  }

  logjpos = 0;
  for(j=1;j<=neq;j++){
    if (((nj_ws->absFdelRm[j]<TINY)&&(nj_ws->absFvalueRm[j] < TINY))||(fabs(nj_ws->Difmax[j])<TINY)){
      nj_ws->logj[j] = 1;/*.true.*/
      logjpos = 1;
    }
    else{
      nj_ws->logj[j] = 0;
    }
  }

  if (logjpos ==1){
    for(i=1;i<=neq;i++){
      nj_ws->yydel[i] = y[i];
      nj_ws->Fscale[i] = max(nj_ws->absFdelRm[i],nj_ws->absFvalueRm[i]);
    }
    /* If the difference in f values is so small that the column might be just
       ! roundoff error, try a bigger increment. */
    for(j=1;j<=neq;j++){
      if ((nj_ws->logj[j]==1)&&(nj_ws->Difmax[j]<=(br*nj_ws->Fscale[j]))){
	tmpfac = min(sqrt(fac[j]),facmax);
	del2 = (y[j] + tmpfac*nj_ws->yscale[j]) - y[j];
	if ((tmpfac!=fac[j])&&(del2!=0.0)){
	  if (fval[j] >= 0.0){
	    /*! keep del pointing into region */
	    del2 = fabs(del2);
	  }
	  else{
	    del2 = -fabs(del2);
	  }
	  nj_ws->yydel[j] = y[j] + del2;
	  lasagna_call((*derivs)(t,nj_ws->yydel+1,nj_ws->ffdel+1,
			       parameters_and_workspace_for_derivs,error_message),
		     error_message,error_message);
	  *nfe+=1;
	  nj_ws->yydel[j] = y[j];
	  rowmax2 = 1;
	  Fdiff_new=0.0;
	  Fdiff_absrm = 0.0;
	  for(i=1;i<=neq;i++){
	    Fdiff_absrm = max(Fdiff_absrm,fabs(Fdiff_new));
	    Fdiff_new = nj_ws->ffdel[i]-fval[i];
	    nj_ws->tmp[i] = Fdiff_new/del2;
	    if(fabs(Fdiff_new)>=Fdiff_absrm){
	      rowmax2 = i;
	      difmax2 = fabs(Fdiff_new);
	    }
	  }
	  maxval1 = difmax2*fabs(del2)*tmpfac;
	  maxval2 = nj_ws->Difmax[j]*fabs(nj_ws->del[j]);
	  if(maxval1>=maxval2){
	    /* The new difference is more significant, so
	       use the column computed with this increment.
	       This depends on wether we are in sparse mode or not: */
	    switch(J->Stype){
	    case(L_SCC):
	      for(i=Ap[j-1];i<Ap[j];i++) 
		Ax[i]=nj_ws->tmp[Ai[i]];
	      break;
	    case(L_DNR):
	      for(i=1;i<=neq;i++) 
		dFdy[i][j]=nj_ws->tmp[i];
	      break;
	    }
	    /* Adjust fac for the next call to numjac. */
	    ffscale = max(fabs(nj_ws->ffdel[rowmax2]),nj_ws->absFvalue[rowmax2]);
	    if (difmax2 <= bl*ffscale){
	      /* The difference is small, so increase the increment. */
	      fac[j] = min(10*tmpfac, facmax);
	    }
	    else if(difmax2 > bu*ffscale){
	      /* The difference is large, so reduce the increment. */
	      fac[j] = max(0.1*tmpfac, facmin);
	    }
	    else{
	      fac[j] = tmpfac;
	    }
	  }
	}
      }
    }
  }
  return _SUCCESS_;
} /* End of numjac */


int initialize_numjac_workspace(MultiMatrix *J,
				void **numjac_workspace,
				ErrorMsg error_message){
  struct numjac_workspace * nj_ws;
  SCCformat *StoreSCC; 
  int i,neq = J->ncol, neqp=neq+1;
  
  lasagna_alloc(nj_ws, sizeof(struct numjac_workspace), error_message);
  /* Allocate vectors and matrices: */
  lasagna_alloc(nj_ws->jacvec,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->yscale,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->del,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->Difmax,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->absFdelRm,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->absFvalue,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->absFvalueRm,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->Fscale,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->ffdel,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->yydel,sizeof(double)*neqp,error_message);
  lasagna_alloc(nj_ws->tmp,sizeof(double)*neqp,error_message);
	
  /* Allocate vector of pointers to rows of matrix.*/
  lasagna_alloc(nj_ws->ydel_Fdel,sizeof(double*)*(neq+1),error_message); 
  lasagna_alloc(nj_ws->ydel_Fdel[1],sizeof(double)*(neq*neq+1),error_message);
  nj_ws->ydel_Fdel[0] = NULL;
  for(i=2;i<=neq;i++) nj_ws->ydel_Fdel[i] = nj_ws->ydel_Fdel[i-1]+neq; /* Set row pointers... */ 
	
  lasagna_alloc(nj_ws->logj,sizeof(int)*neqp,error_message);
  lasagna_alloc(nj_ws->Rowmax,sizeof(int)*neqp,error_message);
  /* Done allocating stuff */
  /* Initialize jacvec to sqrt(eps):*/
  for (i=1;i<=neq;i++) nj_ws->jacvec[i]=1.490116119384765597872e-8;
  nj_ws->col_group = NULL;
  if (J->Stype==L_SCC){
    StoreSCC = J->Store;
    lasagna_alloc(nj_ws->col_group, sizeof(int)*neq, error_message);
    nj_ws->max_group = get_column_grouping(StoreSCC->Ap, StoreSCC->Ai, 
					  neq, nj_ws->col_group, nj_ws->Rowmax);
  }
  *numjac_workspace = (void *) nj_ws;
  return _SUCCESS_;
}

int uninitialize_numjac_workspace(void *numjac_workspace){
  struct numjac_workspace * nj_ws = numjac_workspace;
  /* Deallocate vectors and matrices: */
  free(nj_ws->jacvec);
  free(nj_ws->yscale);
  free(nj_ws->del);
  free(nj_ws->Difmax);
  free(nj_ws->absFdelRm);
  free(nj_ws->absFvalue);
  free(nj_ws->absFvalueRm);
  free(nj_ws->Fscale);
  free(nj_ws->ffdel);
  free(nj_ws->yydel);
  free(nj_ws->tmp);

  free(nj_ws->ydel_Fdel[1]);
  free(nj_ws->ydel_Fdel);
  free(nj_ws->logj);
  free(nj_ws->Rowmax);

  if (nj_ws->col_group != NULL)
    free(nj_ws->col_group);

  free(nj_ws);
  return _SUCCESS_;
}
