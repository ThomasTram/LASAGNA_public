/****************************************/
/* Stiff ODE evolver for CLASS					*/
/* 19/11 2010 													*/
/* Thomas Tram					       							*/
/****************************************/
/*	This is an variable order, adaptive stepsize ODE evolver for lasagna.
	It is based on the Numerical Differentiaion Formulas	of order 1-5,
	and can be used to solve stiff problems. The algorithm is described in 
	[The MATLAB ODE Suite, L. F. Shampine and M. W. Reichelt, SIAM Journal 
	on Scientific Computing, 18-1, 1997].
	
	The code will call the (*output) routine at x-values in t_vec[]. It
	will interpolate to get the y-value aswell as (optionally) the first and
	second derivative at points returned by the routine relevant_indices. 
	Since the transfer function only depends on a few of these values, there
	is no need to interpolate all of them.
		
	Every time there is a stepsize change, we need to rebuild **dif, the matrix
	which holds the backward differences, to reflect the new stepsize. This is done
	in "adjust_stepsize" and is essentially a matrix multiplication. Every times 
	there is either a stepsize change, the method order k is changed, or we compute
	a new Jacobian, we must call "new_linearisation" to calculate a new LU 
	decomposition of a matrix.
	
	The method will not recompute the Jacobian in every step, but only when the
	Newton iterations fail to converge fast enough. This feature makes the
	solver competitive even for non-stiff problems.
	
	Statistics is saved in the stepstat[6] vector. The entries are:
	stepstat[0] = Successful steps.
	stepstat[1] = Failed steps.
	stepstat[2] = Total number of function evaluations.
	stepstat[3] = Number of Jacobians computed.
	stepstat[4] = Number of LU decompositions.
	stepstat[5] = Number of linear solves.
	If ppt->perturbations_verbose > 2, this statistic is printed at the end of
	each call to evolver.
	
	Sparsity:
	When the number of equations becomes high, too much times is spent on solving
	linear equations. Since the ODE-functions are not very coupled, one would normally
	supply a sparsity pattern of the jacobian, which would encode the couplings of the
	equations. However, this would be of some inconvenience to the users who just want
	to add some new physics without thinking too much about how the code work. So we
	pay a few more function evaluations, and calculate the full jacobian every time.
	
	Then, if jac->use_sparse==_TRUE_, numjac will try to construct a sparse matrix from 
	the dense matrix. If there are too many nonzero elements in the dense matrix, numjac
	will stop constructing the sparse matrix and set jac->use_sparse=_FALSE_. The sparse
	matrix is stored in the compressed column format. (See sparse.h).
	
	In the sparse case, we also do partial pivoting, but with diagonal preference. The
	structure of the equations are nearly optimal for the LU decomposition, so we don't
	want to mess it up by too many row permutations if we can avoid it. This is also why
	do not use any column permutation to pre-order the matrix.
*/
#include "common.h"
#include "evolver_ndf15.h"
#include "sparse.h"

int evolver_ndf15(int (*derivs)(double x,double * y,double * dy,
				void * parameters_and_workspace, ErrorMsg error_message),
		  void * parameters_and_workspace_for_derivs,
		  double t0,
		  double tfinal,
		  double * y_inout, 
		  size_t neq, 
		  EvolverOptions *options,
		  ErrorMsg error_message){
	
  /** Handle options: */
  int *interpidx, *used_in_output, *stepstat, verbose, tres; 
  double abstol, rtol, *t_vec;
  int (*linalg_initialise)(MultiMatrix *, EvolverOptions *, void **, ErrorMsg);
  int (*linalg_finalise)(void *, ErrorMsg);
  int (*linalg_factorise)(void *, int, ErrorMsg);
  int (*linalg_solve)(MultiMatrix *, MultiMatrix *, void *, ErrorMsg);
  int (*output)(double t, double *y, double *dy, int i, void *p, ErrorMsg err);
  int (*print_variables)(double t, double *y, double *dy, void *p, ErrorMsg err);
  int (*stop_function)(double t, double *y, double *dy, void *p, ErrorMsg err);
  used_in_output = options->used_in_output; stepstat = &(options->Stats[0]);
  verbose = options->EvolverVerbose; tres = options->tres; abstol = options->AbsTol; 
  rtol = options->RelTol; t_vec = options->t_vec; 
  linalg_initialise = options->linalg_initialise; linalg_finalise = options->linalg_finalise;
  linalg_factorise = options->linalg_factorise; linalg_solve = options->linalg_solve;
  output = options->output; print_variables=options->print_variables; 
  stop_function = options->stop_function;
  int *Ai, *Ap, nnz;

  /* Constants: */
  double G[5]={1.0,3.0/2.0,11.0/6.0,25.0/12.0,137.0/60.0};
  double alpha[5]={-37.0/200,-1.0/9.0,-8.23e-2,-4.15e-2, 0};
  double invGa[5],erconst[5];
  int maxit=4, maxk=2;
  double threshold = abstol/rtol;
	
  /* Logicals: */
  int Jcurrent,new_jacobian, havrate,done,at_hmin,nofailed,gotynew,tooslow;
	
  /* Storage: */
  double *f0,*y,*wt,*ddfddt,*pred,*ynew,*invwt,*rhs,*psi,*difkp1,*del,*yinterp;
  double *tempvec1,*tempvec2,*ypinterp,*yppinterp;
  double **dif;
	
  /* Method variables: */
  double t,ti,tnew=0;
  double rh,htspan,absh,hmin,hmax,h,tdel;
  double abshlast,hinvGak,minnrm,oldnrm=0.,newnrm;
  double err,hopt,errkm1,hkm1,errit,rate=0.,temp,errkp1,hkp1,maxtmp;
  int k,klast,nconhk,iter,next=0,kopt,tdir;
	
  /* Misc: */
  int nfenj,j,ii,jj, numidx;
  size_t neqp=neq+1;

  /* Matrices for jacobian and linearisation: */
  MultiMatrix J, A, RHS, DEL;
  double *Jval, *Aval;
  void *linalg_workspace_A, *nj_ws;
  double **Matrix;
  DNRformat *StoreDNR;

  /** Allocate memory . */
  void * buffer;
  int buffer_size;

  buffer_size=
    15*neqp*sizeof(double)
    +neqp*sizeof(int)
    +neqp*sizeof(double*)
    +(7*neq+1)*sizeof(double);

  lasagna_alloc(buffer,
	      buffer_size,
	      error_message);

  f0       =(double*)buffer;
  wt       =f0+neqp;
  ddfddt   =wt+neqp;
  pred     =ddfddt+neqp;
  y        =pred+neqp;
  invwt    =y+neqp;
  rhs      =invwt+neqp;
  psi      =rhs+neqp;
  difkp1   =psi+neqp;
  del      =difkp1+neqp;
  yinterp  =del+neqp;
  ypinterp =yinterp+neqp;
  yppinterp=ypinterp+neqp;
  tempvec1 =yppinterp+neqp;
  tempvec2 =tempvec1+neqp;

  interpidx=(int*)(tempvec2+neqp);

  dif      =(double**)(interpidx+neqp);
  dif[1]   =(double*)(dif+neqp);
  for(j=2;j<=neq;j++) dif[j] = dif[j-1]+7; /* Set row pointers... */
  dif[0] = NULL;
  /* for (ii=0;ii<(7*neq+1);ii++) dif[1][ii]=0.; */
  for (j=1; j<=neq; j++) {
    for (ii=1;ii<=7;ii++) {
      dif[j][ii]=0.;
    }
  }

  /* 	lasagna_alloc(f0,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(wt,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(ddfddt,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(pred,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(y,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(invwt,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(rhs,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(psi,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(difkp1,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(del,sizeof(double)*neqp,error_message); */ 
  /* 	lasagna_alloc(yinterp,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(ypinterp,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(yppinterp,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(tempvec1,sizeof(double)*neqp,error_message); */
  /* 	lasagna_alloc(tempvec2,sizeof(double)*neqp,error_message); */
	
  /* 	lasagna_alloc(interpidx,sizeof(int)*neqp,error_message); */
	
  /* Allocate vector of pointers to rows of dif:*/
  /* 	lasagna_alloc(dif,sizeof(double*)*neqp,error_message);  */
  /* 	lasagna_calloc(dif[1],(7*neq+1),sizeof(double),error_message); */
  /* 	dif[0] = NULL; */
  /* 	for(j=2;j<=neq;j++) dif[j] = dif[j-1]+7; */ /* Set row pointers... */ 
 
  /*Set pointers:*/
  ynew = y_inout-1; /* This way y_inout is always up to date. */

  /** Initialise MultiMatrix J, A and the linear method: */
  if (options->use_sparse == _TRUE_){
    printf("Use Sparse\n");
    Ai=options->Ai; Ap=options->Ap;
    nnz = Ap[neq];

    lasagna_calloc(Jval, nnz, sizeof(double), error_message);
    lasagna_alloc(Aval, sizeof(double)*nnz, error_message);
    lasagna_call(CreateMatrix_SCC(&(J),L_DBL,neq, neq, nnz, Ai, Ap, Jval, error_message),
		 error_message, error_message);
    lasagna_call(CreateMatrix_SCC(&(A),L_DBL,neq, neq, nnz, Ai, Ap, Aval, error_message),
		 error_message, error_message);
  }
  else{
    printf("Use dense\n");
    lasagna_calloc(Jval, (neq*neq+1), sizeof(double), error_message);
    lasagna_alloc(Aval, sizeof(double)*(neq*neq+1), error_message);
    lasagna_call(CreateMatrix_DNR(&(J),L_DBL,neq, neq, Jval, error_message),
		 error_message, error_message);
    lasagna_call(CreateMatrix_DNR(&(A),L_DBL,neq, neq, Aval, error_message),
		 error_message, error_message);
  }
  if(options->J_pointer_flag ==  _TRUE_) options->J_pointer = &(J);
  lasagna_call(linalg_initialise(&A, options, &linalg_workspace_A,error_message),
	       error_message, error_message);
  lasagna_call(CreateMatrix_DNR(&(RHS), L_DBL, 1, neq, rhs, error_message),
	       error_message, error_message);
  lasagna_call(CreateMatrix_DNR(&(DEL), L_DBL, 1, neq, del, error_message),
	       error_message, error_message);

  /* Initialize workspace for numjac: */
  lasagna_call(initialize_numjac_workspace(&J, &nj_ws,error_message),
	       error_message,error_message);
	
  /* Initialize some method parameters:*/
  for(ii=0;ii<5;ii++){
    invGa[ii] = 1.0/(G[ii]*(1.0 - alpha[ii]));
    erconst[ii] = alpha[ii]*G[ii] + 1.0/(2.0+ii);
  }

  /* Set the relevant indices which needs to be found by interpolation. */
  /* But if we want to print variables for testing purposes, just interpolate everything.. */
  for(ii=1;ii<=neq;ii++){
    y[ii] = y_inout[ii-1];
    if ((used_in_output != NULL)&&(print_variables==NULL)){
      interpidx[ii]=used_in_output[ii-1];
    }
    else{
      interpidx[ii]=1;
    }
  }

  /* Some lasagna-specific stuff:*/
  if ((tfinal-t0)<0.0){
    tdir = -1;
  }
  else{
    tdir = 1;
  }

  if (t_vec!=NULL){
    //Do output at specified locations
    for(next=0; (t_vec[next]-t0)*tdir<0.0; next++);
  }
 	
  if (verbose > 3){
    numidx=0;
    for(ii=1;ii<=neq;ii++){
      if (interpidx[ii]==_TRUE_) numidx++;
    }
    printf("%d/%d\n",numidx,(int) neq);
  }
  
  htspan = fabs(tfinal-t0);
  hmax = htspan/10.0;

  for(ii=0;ii<6;ii++) stepstat[ii] = 0;
  
  lasagna_call((*derivs)(t0,
			 y+1,
			 f0+1,
			 parameters_and_workspace_for_derivs,
			 error_message),
	       error_message,error_message);
  stepstat[2] +=1;

  t = t0;
  nfenj=0;
  lasagna_call(numjac((*derivs),
		      t,
		      y,
		      f0,
		      &J,
		      nj_ws,
		      abstol,
		      neq,
		      &nfenj,
		      parameters_and_workspace_for_derivs,
		      error_message),
	       error_message,error_message);  
  if(options->J_pointer_flag == _TRUE_){
    // Setting jacvec to default value.
    for(j=1;j<=neq;j++) ((struct numjac_workspace*) nj_ws)->jacvec[j]=1.490116119384765597872e-8;
    // Calling derivs and numjac to ensure a updated Jacobian.
    lasagna_call((*derivs)(t0,
			   y+1,
			   f0+1,
			   parameters_and_workspace_for_derivs,
			   error_message),
		 error_message,error_message);
    stepstat[2] +=1;
    lasagna_call(numjac((*derivs),
			t,
			y,
			f0,
			&J,
			nj_ws,
			abstol,
			neq,
			&nfenj,
			parameters_and_workspace_for_derivs,
			error_message),
		 error_message,error_message);
    stepstat[3] += 1;
  }
  stepstat[3] += 1;
  stepstat[2] += nfenj;
  Jcurrent = _TRUE_; 
  new_jacobian = _TRUE_;
	
  hmin = 16.0*DBL_EPSILON*fabs(t);
  /*Calculate initial step */
  rh = 0.0;

  for(jj=1;jj<=neq;jj++){
    wt[jj] = max(fabs(y[jj]),threshold);
    /*printf("wt: %4.8f \n",wt[jj]);*/
    rh = max(rh,1.25/sqrt(rtol)*fabs(f0[jj]/wt[jj]));
    //printf("Index: %d, rh=%e\n",jj,1.25/sqrt(rtol)*fabs(f0[jj]/wt[jj]));
  }
  
  absh = min(hmax, htspan);
  if (absh * rh > 1.0) absh = 1.0 / rh;
	
  absh = max(absh, hmin);
  h = tdir * absh;
  tdel = (t + tdir*min(sqrt(DBL_EPSILON)*max(fabs(t),fabs(t+h)),absh)) - t;

  lasagna_call((*derivs)(t+tdel,y+1,tempvec1+1,parameters_and_workspace_for_derivs,error_message),
	     error_message,error_message);
  stepstat[2] += 1;

  /*I assume that a full jacobi matrix is always calculated in the beginning...*/
  //Must do something here:
  switch(J.Stype){
  case(L_DNR):
    StoreDNR = (DNRformat *) J.Store;
    Matrix = (double **) StoreDNR->Matrix;
    for(ii=1;ii<=neq;ii++){
      ddfddt[ii]=0.0;
      for(jj=1;jj<=neq;jj++){
	ddfddt[ii]+=(Matrix[ii][jj])*f0[jj];
      }
    }
    rh = 0.0;
    for(ii=1;ii<=neq;ii++){
      ddfddt[ii] += (tempvec1[ii] - f0[ii]) / tdel;
      rh = max(rh,1.25*sqrt(0.5*fabs(ddfddt[ii]/wt[ii])/rtol));
    }
    absh = min(hmax, htspan);
    if (absh * rh > 1.0) absh = 1.0 / rh;
    absh = max(absh, hmin);
    h = tdir * absh;
    break;
  }  
  /* Done calculating initial step
     Get ready to do the loop:*/
  k = 1;			/*start at order 1 with BDF1	*/
  klast = k;
  abshlast = absh;

  for(ii=1;ii<=neq;ii++) dif[ii][1] = h*f0[ii];
	
  hinvGak = h*invGa[k-1];
  nconhk = 0; 	/*steps taken with current h and k*/
 
  update_linear_system_ndf15(&J, &A, hinvGak);
  lasagna_call(linalg_factorise(linalg_workspace_A, new_jacobian, error_message),
	       error_message, error_message);
  stepstat[4] += 1;
  new_jacobian = _FALSE_;
  havrate = _FALSE_; /*false*/

  /* Doing main loop: */
  done = _FALSE_;
  at_hmin = _FALSE_;
  while (done==_FALSE_){
    hmin = 16*DBL_EPSILON*fabs(t);
    maxtmp = max(hmin,absh);
    absh = min(hmax, maxtmp);
    if (fabs(absh-hmin)<100*DBL_EPSILON){
      /* If the stepsize has not changed */
      if (at_hmin==_TRUE_){
	absh = abshlast;	/*required by stepsize recovery */
      }
      at_hmin = _TRUE_;
    }
    else{
      at_hmin = _FALSE_;
    }
    h = tdir * absh;
    /* Stretch the step if within 10% of tfinal-t. */
    if (1.1*absh >= fabs(tfinal - t)){
      h = tfinal - t;
      absh = fabs(h);
      done = _TRUE_;
    }
    if (((fabs(absh-abshlast)/absh)>1e-6)||(k!=klast)){
      adjust_stepsize(dif,(absh/abshlast),neq,k);
      hinvGak = h * invGa[k-1];
      nconhk = 0;
      update_linear_system_ndf15(&J, &A, hinvGak);
      lasagna_call(linalg_factorise(linalg_workspace_A, new_jacobian, error_message),
		   error_message, error_message);
      stepstat[4] += 1;
      new_jacobian = _FALSE_;
      havrate = _FALSE_;
    }
    /*		Loop for advancing one step */
    nofailed = _TRUE_;
    for( ; ; ){
      gotynew = _FALSE_;	/* is ynew evaluated yet?*/
      while(gotynew==_FALSE_){
	/*Compute the constant terms in the equation for ynew.
	  Next FOR lop is just: psi = matmul(dif(:,1:k),(G(1:k) * invGa(k)))*/
	for(ii=1;ii<=neq;ii++){
	  psi[ii] = 0.0;
	  for(jj=1;jj<=k;jj++){
	    psi[ii] += dif[ii][jj]*G[jj-1]*invGa[k-1];
	  }
	}
	/* Predict a solution at t+h. */
	tnew = t + h;
	if (done==_TRUE_){
	  tnew = tfinal; /*Hit end point exactly. */
	}
	h = tnew - t; 		 /* Purify h. */
	for(ii=1;ii<=neq;ii++){
	  pred[ii] = y[ii];
	  for(jj=1;jj<=k;jj++){
	    pred[ii] +=dif[ii][jj];
	  }
	}
	eqvec(pred,ynew,neq);
							
	/*The difference, difkp1, between pred and the final accepted
	  ynew is equal to the backward difference of ynew of order
	  k+1. Initialize to zero for the iteration to compute ynew.
	*/
							
	minnrm = 0.0;
	for(j=1;j<=neq;j++){
	  difkp1[j] = 0.0;
	  maxtmp = max(fabs(ynew[j]),fabs(y[j]));
	  invwt[j] = 1.0 / max(maxtmp,threshold);
	  maxtmp = 100*DBL_EPSILON*fabs(ynew[j]*invwt[j]);
	  minnrm = max(minnrm,maxtmp);
	}
	/* Iterate with simplified Newton method. */
	tooslow = _FALSE_;
	for(iter=1;iter<=maxit;iter++){
	  for (ii=1;ii<=neq;ii++){
	    tempvec1[ii]=(psi[ii]+difkp1[ii]);
	  }
	  lasagna_call((*derivs)(tnew,ynew+1,f0+1,parameters_and_workspace_for_derivs,error_message),
		       error_message,error_message);
	  stepstat[2] += 1;
	  for(j=1;j<=neq;j++){
	    rhs[j] = hinvGak*f0[j]-tempvec1[j];
	  }
								
	  /*Solve the linear system A*x=del by using the LU decomposition stored in linalg_workspace.*/
	  lasagna_call(linalg_solve(&RHS, &DEL, linalg_workspace_A, error_message),
		       error_message, error_message);
	  stepstat[5]+=1;
	  newnrm = 0.0;
	  for(j=1;j<=neq;j++){
	    maxtmp = fabs(del[j]*invwt[j]);
	    newnrm = max(newnrm,maxtmp);
	  }
	  for(j=1;j<=neq;j++){
	    difkp1[j] += del[j];
	    ynew[j] = pred[j] + difkp1[j];
	  }
	  if (newnrm <= minnrm){
	    gotynew = _TRUE_;
	    break; /* Break Newton loop */
	  }
	  else if(iter == 1){
	    if (havrate==_TRUE_){
	      errit = newnrm * rate / (1.0 - rate);
	      if (errit <= 0.05*rtol){
		gotynew = _TRUE_;
		break; /* Break Newton Loop*/
	      }
	    }
	    else {
	      rate = 0.0;
	    }
	  }
	  else if(newnrm > 0.9*oldnrm){
	    tooslow = _TRUE_;
	    break; /*Break Newton lop */
	  }
	  else{
	    rate = max(0.9*rate, newnrm / oldnrm);
	    havrate = _TRUE_;
	    errit = newnrm * rate / (1.0 - rate);
	    if (errit <= 0.5*rtol){
	      gotynew = _TRUE_;
	      break; /* exit newton */
	    }
	    else if (iter == maxit){
	      tooslow = _TRUE_;
	      break; /*exit newton */
	    }
	    else if (0.5*rtol < errit*pow(rate,(maxit-iter))){
	      tooslow = _TRUE_;
	      break; /*exit Newton */
	    }
	  }
	  oldnrm = newnrm;
	}
	if (tooslow==_TRUE_){
	  stepstat[1] += 1;
	  /*	! Speed up the iteration by forming new linearization or reducing h. */
	  if (Jcurrent==_FALSE_){
	    lasagna_call((*derivs)(t,y+1,f0+1,parameters_and_workspace_for_derivs,error_message),
		       error_message,error_message);
	    nfenj=0;
	    lasagna_call(numjac((*derivs),t,y,f0,&J,nj_ws,abstol,neq,
			      &nfenj,parameters_and_workspace_for_derivs,error_message),
			 error_message,error_message);
	    if(options->J_pointer_flag == _TRUE_){
	      // Calling derivs and numjac to ensure a updated Jacobian.
	      lasagna_call((*derivs)(t,y+1,f0+1, parameters_and_workspace_for_derivs,error_message),
			   error_message,error_message);
	      stepstat[2] +=1;
	      lasagna_call(numjac((*derivs),t,y,f0,&J,nj_ws,abstol,neq,
				  &nfenj,parameters_and_workspace_for_derivs,error_message),
			   error_message,error_message);
	      stepstat[3] += 1;
	    }
	    stepstat[3] += 1;
	    stepstat[2] += (nfenj + 1);
	    Jcurrent = _TRUE_;
	    new_jacobian = _TRUE_;
	  }
	  else if (absh <= hmin){
	    lasagna_test(absh <= hmin, error_message,
		       "Step size too small: step:%g, minimum:%g, in interval: [%g:%g]\n",
		       absh,hmin,t0,tfinal);
	  }
	  else{
	    abshlast = absh;
	    absh = max(0.3 * absh, hmin);
	    h = tdir * absh;
	    done = _FALSE_;
	    adjust_stepsize(dif,(absh/abshlast),neq,k);
	    hinvGak = h * invGa[k-1];
	    nconhk = 0;
	  }
	  /* A new linearisation is needed in both cases */
	  update_linear_system_ndf15(&J, &A, hinvGak);
	  lasagna_call(linalg_factorise(linalg_workspace_A, new_jacobian, error_message),
		       error_message, error_message);
	  stepstat[4] += 1;
	  new_jacobian = _FALSE_;
	  havrate = _FALSE_;
	}
      }
      /*end of while loop for getting ynew
	difkp1 is now the backward difference of ynew of order k+1. */
      err = 0.0;
      int maxerr=0;
      for(jj=1;jj<=neq;jj++){
	if (err<fabs(difkp1[jj]*invwt[jj])){
	  maxerr = jj-1;
	}
	err = max(err,fabs(difkp1[jj]*invwt[jj]));
      }
      //printf("Max Err index: %d\n",maxerr);
      err = err * erconst[k-1];
      if (verbose>3)
	printf("%e %d %e %e .. v = [%g, %g]\n",
	       t,maxerr,err,absh,ynew[maxerr+1],ynew[1]);
      if (err>rtol){
	/*Step failed */
	stepstat[1]+= 1;
	if (absh <= hmin){
	  lasagna_test(absh <= hmin, error_message,
		     "Step size too small: step:%g, minimum:%g, in interval: [%g:%g]\n",
		     absh,hmin,t0,tfinal);
	}
	abshlast = absh;
	if (nofailed==_TRUE_){
	  nofailed = _FALSE_;
	  hopt = absh * max(0.1, 0.833*pow((rtol/err),(1.0/(k+1))));
	  if (k > 1){
	    errkm1 = 0.0;
	    for(jj=1;jj<=neq;jj++){
	      errkm1 = max(errkm1,fabs((dif[jj][k]+difkp1[jj])*invwt[jj]));
	    }
	    errkm1 = errkm1*erconst[k-2];
	    hkm1 = absh * max(0.1, 0.769*pow((rtol/errkm1),(1.0/k)));
	    if (hkm1 > hopt){
	      hopt = min(absh,hkm1); 		/* don't allow step size increase */
	      k = k - 1;
	    }
	  }
	  absh = max(hmin, hopt);
	}
	else{
	  absh = max(hmin, 0.5 * absh);
	}
	h = tdir * absh;
	if (absh < abshlast){
	  done = _FALSE_;
	}
	adjust_stepsize(dif,(absh/abshlast),neq,k);
	hinvGak = h * invGa[k-1];
	nconhk = 0;
	update_linear_system_ndf15(&J, &A, hinvGak);
	lasagna_call(linalg_factorise(linalg_workspace_A, new_jacobian, error_message),
		     error_message, error_message);
	stepstat[4] += 1;
	new_jacobian = _FALSE_;
	havrate = _FALSE_;
      }
      else {
	break; /* Succesfull step */
      }
    }
    /* End of conditionless FOR loop */
    if (print_variables != NULL){
      lasagna_call((*print_variables)(t,
				      ynew+1,
				      f0+1,
				      parameters_and_workspace_for_derivs,
				      error_message),
		   error_message,error_message);
    }
    stepstat[0] += 1;
		 
    /* Update dif: */
    for(jj=1;jj<=neq;jj++){
      dif[jj][k+2] = difkp1[jj] - dif[jj][k+1];
      dif[jj][k+1] = difkp1[jj];
    }
    for(j=k;j>=1;j--){
      for(ii=1;ii<=neq;ii++){
	dif[ii][j] += dif[ii][j+1];
      }
    }
    /** Output **/
    if (t_vec==NULL){
      //Refinement output:
      for (jj=1; jj<tres; jj++){
	//Interpolated outputs:
	ti = tnew-(1.0-jj/((double) tres))*h;
	interp_from_dif(ti,
			tnew,
			ynew,
			h,
			dif,
			k,
			yinterp,
			ypinterp,
			yppinterp,
			interpidx,
			neq,
			2);				
	lasagna_call((*output)(ti,
			       yinterp+1,
			       ypinterp+1,
			       jj,
			       parameters_and_workspace_for_derivs,
			       error_message),error_message,error_message);
      }
      lasagna_call((*output)(tnew,
			     ynew+1,
			     f0+1,
			     tres,
			     parameters_and_workspace_for_derivs,
			     error_message),error_message,error_message);
    }
    else {
      //Output at Tvec grid:
      while ((next<tres)&&(tdir * (tnew - t_vec[next]) >= 0.0)){
	/* Do we need to write output? */
	if (tnew==t_vec[next]){
	  lasagna_call((*output)(t_vec[next],
				 ynew+1,
				 f0+1,
				 next,
				 parameters_and_workspace_for_derivs,
				 error_message),
		       error_message,error_message);
	  if (print_variables != NULL){
	    lasagna_call((*print_variables)(t_vec[next],
					    ynew+1,
					    f0+1,
					    parameters_and_workspace_for_derivs,
					    error_message),
			 error_message,error_message);
	  }
	}
	else {
	  /*Interpolate if we have overshot sample values*/
	  interp_from_dif(t_vec[next],
			  tnew,
			  ynew,
			  h,
			  dif,
			  k,
			  yinterp,
			  ypinterp,
			  yppinterp,
			  interpidx,
			  neq,
			  2);				
	  lasagna_call((*output)(t_vec[next],
				 yinterp+1,
				 ypinterp+1,
				 next,
				 parameters_and_workspace_for_derivs,
				 error_message),error_message,error_message);
	}
	next++;	
      }
    }
    /** End of output **/
    if (done==_TRUE_) {
      break;
    }
    klast = k;
    abshlast = absh;
    nconhk = min(nconhk+1,maxk+2);
    if (nconhk >= k + 2){
      temp = 1.2*pow((err/rtol),(1.0/(k+1.0)));
      if (temp > 0.1){
	hopt = absh / temp;
      }
      else {
	hopt = 10*absh;
      }
      kopt = k;
      if (k > 1){
	errkm1 = 0.0;
	for(jj=1;jj<=neq;jj++){
	  errkm1 = max(errkm1,fabs(dif[jj][k]*invwt[jj]));
	}
	errkm1 = errkm1*erconst[k-2];
	temp = 1.3*pow((errkm1/rtol),(1.0/k));
	if (temp > 0.1){
	  hkm1 = absh / temp;
	}
	else {
	  hkm1 = 10*absh;
	}
	if (hkm1 > hopt){
	  hopt = hkm1;
	  kopt = k - 1;
	}
      }
      if (k < maxk){
	errkp1 = 0.0;
	for(jj=1;jj<=neq;jj++){
	  errkp1 = max(errkp1,fabs(dif[jj][k+2]*invwt[jj]));
	}
	errkp1 = errkp1*erconst[k];
	temp = 1.4*pow((errkp1/rtol),(1.0/(k+2.0)));
	if (temp > 0.1){
	  hkp1 = absh / temp;
	}
	else {
	  hkp1 = 10*absh;
	}
	if (hkp1 > hopt){
	  hopt = hkp1;
	  kopt = k + 1;
	}
      }
      if (hopt > absh){
	absh = hopt;
	if (k!=kopt){
	  k = kopt;
	}
      }
    }
    /* Advance the integration one step. */
    t = tnew;
    eqvec(ynew,y,neq);
    Jcurrent = _FALSE_;

    /* Perhaps use stop function: */
    if (stop_function != NULL){
      if ((stepstat[0]>500000000)||
	  (stop_function(t,y+1,f0+1,parameters_and_workspace_for_derivs,
			   error_message) == _TRUE_)){      //Stop condition
	lasagna_call((*output)(t,y+1,f0+1,next,
			       parameters_and_workspace_for_derivs,
			       error_message),error_message,error_message);
	printf("Stop condition met...\n");
	break;
      }
    }
  }

  /* a last call is compulsory to ensure that all quantitites in
     y,dy,parameters_and_workspace_for_derivs are updated to the
     last point in the covered range */
  printf("Last call to derivs at t=%.16e.\n",tnew);
  lasagna_call(
	     (*derivs)(tnew,
		       ynew+1,
		       f0+1,
		       parameters_and_workspace_for_derivs,error_message),
	     error_message,
	     error_message);
	
  if (verbose > 0){
    printf("\n End of evolver. Next=%d, t=%e and tnew=%e.",next,t,tnew);
    printf("\n Statistics: [%d %d %d %d %d %d] \n",stepstat[0],stepstat[1],
	   stepstat[2],stepstat[3],stepstat[4],stepstat[5]);
  }
	
  /** Deallocate memory */

  free(buffer);
  DestroyMultiMatrix(&J);
  DestroyMultiMatrix(&A);
  DestroyMultiMatrix(&DEL);
  DestroyMultiMatrix(&RHS);
  lasagna_call(linalg_finalise(linalg_workspace_A, error_message),
	       error_message, error_message);
  free(Jval);
  free(Aval);
  /* 	free(f0); */
  /* 	free(wt); */
  /* 	free(ddfddt); */
  /* 	free(pred); */
  /* 	free(y); */
  /* 	free(invwt); */
  /* 	free(rhs); */
  /* 	free(psi); */
  /* 	free(difkp1); */
  /* 	free(del); */
  /* 	free(yinterp); */
  /* 	free(ypinterp); */
  /* 	free(yppinterp); */
  /* 	free(tempvec1); */
  /* 	free(tempvec2); */

  /* 	free(interpidx); */
  /* 	free(dif[1]); */
  /* 	free(dif); */
	
  uninitialize_numjac_workspace(nj_ws);
  return _SUCCESS_;

} /*End of program*/

/**********************************************************************/
/* Here are some small routines used in evolver_ndf15:                */
/* "interp_from_dif", "eqvec", "adjust_stepsize", "calc_C",           */
/* "new_linearisation", "relevant_indices", "ludcmp", "lubksb".       */
/**********************************************************************/

void eqvec(double *datavec,double *emptyvec, int n){
  int i;
  for(i=1;i<=n;i++){
    emptyvec[i] = datavec[i];
  }
}

/* Subroutine that interpolates from information stored in dif */
int interp_from_dif(double tinterp,double tnew,double *ynew,double h,double **dif,int k, double *yinterp,
		    double *ypinterp, double *yppinterp, int* index, size_t neq, int output){
  /* Output=1: only y_vector. Output=2: y and y prime. Output=3: y, yp and ypp*/
  int i,j,m,l,p,factor;
  double sumj,suml,sump,prodm,s;
  s = (tinterp - tnew)/h;
  if (k==1){
    for(i=1;i<=neq;i++){
      if(index[i]==_TRUE_){
	yinterp[i] = ynew[i] + dif[i][1] * s;
	if (output>1) ypinterp[i] = dif[i][1]/h;
	if (output>2) yppinterp[i] = 0; /* No good estimate can be made of the second derivative */
      }
    }
  }
  else{
    /*This gets tricky */
    for(i=1;i<=neq;i++){
      if(index[i]==_TRUE_){
	/*First the value of the function:	*/
	sumj=0.0;
	factor=1;
	for(j=1;j<=k;j++){
	  prodm=1.0;
	  factor*=j;
	  for(m=0;m<j;m++) prodm*=(m+s);
	  sumj+=dif[i][j]/factor*prodm;
	}
	yinterp[i] = ynew[i]+sumj;
	/* Now the first derivative: */
	if (output>1){
	  factor = 1;
	  sumj=0.0;
	  for(j=1;j<=k;j++){
	    suml = 0.0;
	    factor *=j;
	    for(l=0;l<j;l++){
	      prodm=1.0;
	      for(m=0;m<j;m++){
		if(m!=l) prodm*=(m+s);
	      }
	      suml+=prodm;
	    }
	    sumj+=dif[i][j]/factor*suml;
	  }
	  ypinterp[i] = sumj/h;
	}
	/* The second derivative: */
	if (output>2){
	  factor=1;
	  sumj=0.0;
	  for(j=1;j<=k;j++){
	    suml=0.0;
	    factor*=j;
	    for(l=0;l<j;l++){
	      sump=0.0;
	      for(p=0;p<j;p++){
		if(p!=l){
		  prodm=1.0;
		  for(m=0;m<j;m++){
		    if((m!=l)&&(m!=p)){
		      prodm*=(m+s);
		    }
		  }
		  sump+=prodm;
		}
	      }
	      suml+=sump;
	    }
	    sumj+=dif[i][j]/factor*suml;
	  }
	  yppinterp[i] = sumj/(h*h);
	}
      }
    }
  }
  return _SUCCESS_;
}

int adjust_stepsize(double **dif, double abshdivabshlast, size_t neq,int k){
  double mydifU[5][5]={{-1,-2,-3,-4,-5},{0,1,3,6,10},{0,0,-1,-4,-10},{0,0,0,1,5},{0,0,0,0,-1}};
  double tempvec[5];
  double mydifRU[5][5];
  int ii,jj,kk;
	
  for(ii=1;ii<=5;ii++) mydifRU[0][ii-1] = -ii*abshdivabshlast;
  for(jj=2;jj<=5;jj++){
    for(ii=1;ii<=5;ii++){
      mydifRU[jj-1][ii-1] = mydifRU[jj-2][ii-1]*(1.0-(1.0+ii*abshdivabshlast)/jj);
    }
  }
  for(ii=0;ii<5;ii++){
    for(kk=0;kk<5;kk++){
      /* Save the i'th row of mydifRU */
      tempvec[kk] = mydifRU[ii][kk];
    }
    for(jj=0;jj<5;jj++){
      /* Now do the matrix multiplication: */
      mydifRU[ii][jj] = 0.0;
      for(kk=0;kk<5;kk++)	mydifRU[ii][jj] += tempvec[kk]*mydifU[kk][jj];
    }
  }

  for(ii=0;ii<neq;ii++){
    for(kk=0;kk<k;kk++){
      /* Save the k first values of the i'th row of dif */
      tempvec[kk] = dif[ii+1][kk+1];
    }
    for(jj=0;jj<k;jj++){
      /* Now do the matrix multiplication: */
      dif[ii+1][jj+1] = 0.0;
      for(kk=0;kk<k;kk++) dif[ii+1][jj+1] += tempvec[kk]*mydifRU[kk][jj];
    }
  }
  return _SUCCESS_;
}

int update_linear_system_ndf15(MultiMatrix *J, 
			       MultiMatrix *A, 
			       double hinvGak){
  size_t neq=J->ncol;
  double luparity, *Ax, *Jx;
  int i,j,*Ap,*Ai,funcreturn;
  SCCformat *JStoreSCC,*AStoreSCC;
  DNRformat *JStoreDNR,*AStoreDNR;
  double **Jmat, **Amat;
  switch(J->Stype){
  case(L_SCC):
    JStoreSCC = J->Store; AStoreSCC = A->Store;
    Ap = AStoreSCC->Ap; Ai = AStoreSCC->Ai; 
    Ax = AStoreSCC->Ax; Jx = JStoreSCC->Ax;
    /* Construct Ax from Jx, the jacobian:*/
    for(j=0;j<neq;j++){
      for(i=Ap[j];i<Ap[j+1];i++){
	if(Ai[i]==j){
	  /* I'm at the diagonal */
	  Ax[i] = 1.0-hinvGak*Jx[i];
	}
	else{
	  Ax[i] = -hinvGak*Jx[i];
	}
      }
    }
    break;
  case (L_DNR):
    /* Normal calculation: */
    JStoreDNR = J->Store; AStoreDNR = A->Store;
    Jmat = (double **) JStoreDNR->Matrix;
    Amat = (double **) AStoreDNR->Matrix;    
    for(i=1;i<=neq;i++){
      for(j=1;j<=neq;j++){
	Amat[i][j] = - hinvGak * Jmat[i][j];
	if(i==j) Amat[i][j] +=1.0;
      }
    }
  }
  return _SUCCESS_;
}

/** Helper functions */

