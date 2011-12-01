#include "common.h"
#include "evolver_radau5.h"
/**
   Statistics is saved in the stepstat[6] vector. The entries are:
   stepstat[0] = Successful steps.
   stepstat[1] = Failed steps.
   stepstat[2] = Total number of function evaluations.
   stepstat[3] = Number of Jacobians computed.
   stepstat[4] = Number of LU decompositions.
   stepstat[5] = Number of linear solves.
*/
int evolver_radau5(
		  int (*derivs)(double x,
				double * y,
				double * dy,
				void * parameters_and_workspace, 
				ErrorMsg error_message),
		  double t_ini,
		  double t_final,
		  double * y0, 
		  int * interpidx,
		  int neq, 
		  void * parameters_and_workspace_for_derivs,
		  double rtol, 
		  double abstol, 
		  double * t_vec, 
		  int tres,
		  int *Ap,
		  int *Ai,
		  int (*output)(double x,
				double y[],
				double dy[],
				int index_x,
				void * parameters_and_workspace,
				ErrorMsg error_message),
		  int (*print_variables)(double x, 
					 double y[], 
					 double dy[], 
					 void *parameters_and_workspace,
					 ErrorMsg error_message),
		  ErrorMsg error_message){
	
  /* Constants: */
  double Tinv[9]= 
    {4.178718591551904727,   0.327682820761062387,   0.523376445499449548,
     0.502872634945786876,   -2.571926949855605429,   0.596039204828224925,
     -4.17871859155190473,   -0.327682820761062387,   0.476623554500550452};
  double T[9] = 
    {0.09443876248897524149, 0.03002919410514742449, -0.1412552950209542084,
     0.250213122965333311, -0.382942112757261938, 0.204129352293799932,
     1.0,   0.0,   1.0};
  double ci[3] = 
    {(4.0-sqrt(6.0))/10.0, (4.0+sqrt(6.0))/10.0, 1.0};
  double ei[3] = 
    {(-13.0-7.0*sqrt(6.0))/3.0, (-13.0+7.0*sqrt(6.0))/3.0, -1.0/3.0};

  double B[9], C[9], Ainv[9] = 
    {3.22474487139158905,   1.167840084690405495,   -0.253197264742180826,
     -3.56784008469040549,   0.775255128608410951,   1.053197264742180826,
     5.53197264742180826,   -7.53197264742180826,   5.0};

  double gamma_hat = 3.0-pow(3.0,1.0/3.0)+pow(3.0,2.0/3.0);
  double alpha_hat = 0.5*(6.0+pow(3.0,1.0/3.0)-pow(3.0,2.0/3.0));
  double beta_hat = 0.5*pow(3.0,1.0/6.0)*(3.0+pow(3.0,2.0/3.0));
  double gamma, alpha, beta, sr6=sqrt(6.0);

  double safety_base = 0.9, step_fac_newt_failed = 0.5;
  double theta_reuse_jacobian=0.1, same_step_lower=0.99, same_step_upper=2.0;
  double threshold = abstol/rtol;
  double tol_newton = rtol*max(100*DBL_EPSILON,min(0.03,sqrt(rtol)));
  int newt_iter_max = 10;

  int stepstat[6] = {0, 0, 0, 0, 0, 0};
  int newt_iter, tdir, next, nfenj;
  double abshlast, absh, abshnew; 
  double theta, theta_k, theta_k_old, eta=1.0, conv_est, conv_cor;
  double t, tdel, h, abshmin, rh;
  double norm_dW, norm_dW_old, norm_err;
  double z1,z2,z3,tau;

  int first_step = _TRUE_, Newton_converged, got_ynew, J_current;
  int reuse_stepsize, reuse_jacobian, last_failed=_TRUE_, break_next;

  double (*error_norm)(double *y, double *err_y, double threshold, int neq);

  int i, j, debug_stop;


  double *W, *dW, *Y0pZ, *rhs, *Fi, *Zlast;
  double *err, *diff, *xtemp, *ytemp, *ynew, *f0, *ylast, *ftmp, *dfdt;
  double *delta_w;
  double *dW_buf, *rhs_buf, *err_buf;
  double complex *rhs_cx, *delta_w_cx, *delta_w_cx_buf;
  struct jacobian jac;
  struct jacobian_plus jac_plus;
  struct numjac_workspace nj_ws;


  /** Debug area: 
   
  int k;
  for (j=0; j<3; j++){
    for (i=0; i<3; i++){
      B[j*3+i] = 0.0;
      for (k=0; k<3; k++){
	B[j*3+i] += Tinv[j*3+k]*Ainv[k*3+i];
      }
    }
  }
  for (j=0; j<3; j++){
    for (i=0; i<3; i++){
      C[j*3+i] = 0.0;
      for (k=0; k<3; k++){
	C[j*3+i] += B[j*3+k]*T[k*3+i];
      }
      printf("%.18e ",C[j*3+i]);
    }
    printf("\n");
  }
  printf("Gamma = %.17e, Alpha = %.17e, beta = %.17e\n",
	 gamma_hat, alpha_hat, beta_hat);
  printf("ci = [%g %g %g]\n",ci[0],ci[1],ci[2]);
  */

  W = malloc(3*neq*sizeof(double));
  dW_buf = malloc((3*neq+1)*sizeof(double));
  dW = dW_buf+1;
  Y0pZ = malloc(3*neq*sizeof(double));
  Fi = malloc(3*neq*sizeof(double));
  Zlast = malloc(3*neq*sizeof(double));
  rhs_buf = malloc((3*neq+1)*sizeof(double));
  rhs = rhs_buf +1;
  err_buf = malloc((neq+1)*sizeof(double));
  err = err_buf+1;
  diff = malloc(neq*sizeof(double));
  xtemp = malloc(neq*sizeof(double));
  ytemp = malloc(neq*sizeof(double));
  ynew = malloc(neq*sizeof(double));
  f0 = malloc(neq*sizeof(double));
  ylast = malloc(neq*sizeof(double));
  ftmp = malloc(neq*sizeof(double));
  dfdt = malloc(neq*sizeof(double));
  delta_w = malloc(neq*sizeof(double));
  rhs_cx = malloc(neq*sizeof(double complex));
  delta_w_cx_buf = malloc((neq+1)*sizeof(double complex));
  delta_w_cx = delta_w_cx_buf + 1;

  error_norm = norm_inf;

  if ((t_final-t_ini)>0)
    tdir = 1;
  else
    tdir = -1;

  for (next=0; next<tres; next++){
    if ((t_vec[next]-t_ini)*tdir>0) break;
  }

  /*Initialize the jacobian:*/
  lasagna_call(initialize_jacobian(&jac,neq,error_message),
	       error_message,error_message);

  /*Initialise extra storage for the complex linear algebra:*/
  initialize_jacobian_plus(&jac, &jac_plus, neq, error_message);
  
  /* Initialize workspace for numjac: */
  lasagna_call(initialize_numjac_workspace(&nj_ws,neq,error_message),
	       error_message,error_message);

  t = t_ini;
  
  /** Find the initial step: */
  lasagna_call((*derivs)(t,
			 y0,
			 f0,
			 parameters_and_workspace_for_derivs,error_message),
	       error_message,
	       error_message);
  stepstat[2]++;
  
  nfenj=0;
  lasagna_call(numjac((*derivs),
		      t,
		      y0-1,
		      f0-1,
		      &jac,
		      &nj_ws,
		      abstol,
		      neq,
		      &nfenj,
		      parameters_and_workspace_for_derivs,
		      error_message),
	       error_message,
	       error_message);
  stepstat[3] += 1;
  stepstat[2] += nfenj;
  J_current = _TRUE_;
  rh = error_norm(f0, y0, threshold, neq);
  rh *=1.25/rtol;

  abshmin = 16.0*fabs(t)*DBL_EPSILON;
  absh = 0.1*fabs(t_final-t_ini);
  if (absh * rh > 1.0) 
    absh = 1.0 / rh;
  absh = max(absh, abshmin);

  h = tdir * absh;
  tdel = (t + tdir*min(sqrt(DBL_EPSILON)*max(fabs(t),fabs(t+h)),absh)) - t;

  lasagna_call((*derivs)(t+tdel,
		       y0,
		       ftmp,
		       parameters_and_workspace_for_derivs,
		       error_message),
	     error_message,
	     error_message);
  stepstat[2] += 1;

  /*I assume that a full jacobi matrix is always calculated in the beginning...*/
  for(i=0; i<neq; i++){
    dfdt[i]=0.0;
    for(j=0; j<neq; j++){
      dfdt[i] += ((jac.dfdy[i+1][j+1])*f0[j]+
		  (ftmp[i] - f0[i]) / tdel);
    }
  }
  rh = error_norm(y0, dfdt, threshold, neq);
  rh = 1.25*sqrt(0.5*rh/rtol);

  absh = fabs(t_final-t_ini);
  if (absh * rh > 1.0) 
    absh = 1.0 / rh;
  absh = max(absh, abshmin);
  h = tdir * absh;
   /* Done calculating initial step
     Get ready to do the loop:*/
  /** If we have a pattern, we apply it here and do one 
      call to numjac with the correct pattern: */
  if ((Ap!=NULL)&&(Ai!=NULL)){
    for (i=0; i<=neq; i++){
      jac.spJ->Ap[i] = Ap[i];
    }
    for (i=0; i<Ap[neq]; i++){
      jac.spJ->Ai[i] = Ai[i];
    }
    jac.pattern_supplied = _TRUE_;
    lasagna_call(numjac((*derivs),
			t,
			y0-1,
			f0-1,
			&jac,
			&nj_ws,
			abstol,
			neq,
			&nfenj,
			parameters_and_workspace_for_derivs,
			error_message),
		 error_message,
		 error_message);
    stepstat[3] += 1;
    stepstat[2] += nfenj;
    J_current = _TRUE_;
    calc_C(&jac);
    /* Calculate the optimal ordering: */
    sp_amd(jac.Cp, jac.Ci, neq, jac.cnzmax,
	   jac.Numerical->q,jac.Numerical->wamd);
  }
  new_linearisation_radau5(&jac, &jac_plus, h, neq, error_message);  
  stepstat[4] +=1;

  //Main loop:
  while ((t_final-t)*tdir>0.0){
    abshmin = 16.0*fabs(t)*DBL_EPSILON;
    lasagna_test(absh<abshmin, error_message,"Step size too small in evolver_radau5 at t=%g!. |h|=%g < |hmin|=%g", t, absh, abshmin);
    got_ynew = _FALSE_;
    //Loop for getting ynew:
    while (got_ynew == _FALSE_){
      //Initialise Newton method;
      Newton_converged = _FALSE_;
      gamma = gamma_hat/h;
      alpha = alpha_hat/h;
      beta = beta_hat/h;
      /** Set initial conditions for simplified Newton iteration.
	  We set the initial condition on Z, then transform to get
	  the initial condition on W, and then we add [y0;y0;y0] to Z so 
	  we don't need the transform W->Z in the first iteration.
      */
      if (first_step == _TRUE_){
	//Set trivial iniital condition:
	for (i=0; i<3*neq; i++)
	  Y0pZ[i] = 0.0;
      }
      else{
	//Extrapolate the colocation polynomial to get starting value:
	for (i=0; i<neq; i++){
	  z1 = Zlast[i]; z2 = Zlast[neq+i]; z3 = Zlast[2*neq+i];
	  for (j=0; j<3; j++){
	    tau = 1.0+ci[j]*absh/abshlast;
	    Y0pZ[j*neq+i] = ylast[i]-y0[i]+
	      25.0*(tau-1)*(sr6+4.0-10.0*tau)*tau*z1/(6.0-9.0*sr6)-
	      25.0*(tau-1)*(sr6-4.0+10.0*tau)*tau*z2/(6.0+9.0*sr6)+
	      tau*(1.0-8.0*tau+10.0*tau*tau)*z3/3.0;
	  }
	}
      }
      //We have a guess for Z stored in Y0pZ. Calculate W:
      transform_C_tensor_I(Tinv,3,Y0pZ,W,NULL,neq);
      //Add y0 x I to Y0pZ:
      for (i=0; i<3; i++){
	for (j=0; j<neq; j++){
	  Y0pZ[i*neq+j] += y0[j];
	}
      }
      debug_stop=1;
      theta = 1e-3;
      abshnew = 0.5*absh;
      for (newt_iter = 1; newt_iter<=newt_iter_max; newt_iter++){
	//printf("Newton iteration: %d/%d.\n",newt_iter,newt_iter_max);
	if (newt_iter != 1){
	  //We must update Y0pZ from W:
	  transform_C_tensor_I(T,3,W,Y0pZ,y0,neq);
	}
	//Call the function at the 3 nodes:
	for (i=0; i<3; i++){
	  lasagna_call((*derivs)(t+ci[i]*h,
				 Y0pZ+i*neq,
				 Fi+i*neq,
				 parameters_and_workspace_for_derivs,error_message),
		       error_message,
		       error_message);
	}
	stepstat[2] += 3;
	//Start forming the right hand side by doing transformation on Fi:
	transform_C_tensor_I(Tinv, 3, Fi, rhs, NULL, neq);
	for (i=0; i<neq; i++){
	  rhs[i] -= gamma*W[i];
	  rhs_cx[i] = rhs[neq+i]-alpha*W[neq+i]+beta*W[2*neq+i]+
	    I*(rhs[2*neq+i]-beta*W[neq+i]-alpha*W[2*neq+i]);
	}
	//Use backsubstitution to calculate delta W:
	if (jac.use_sparse){
	  sp_lusolve(jac.Numerical, rhs, dW);
	  sp_lusolve_cx(jac_plus.Numerical_cx, rhs_cx, delta_w_cx);
	}
	else{
	  for (i=0; i<neq; i++){
	    dW[i] = rhs[i];
	    delta_w_cx[i] = rhs_cx[i];
	  }
	  lubksb(jac.LU,neq,jac.luidx,dW-1);
	  lubksb_cx(jac_plus.LU_cx,neq,jac_plus.luidx,delta_w_cx-1);
	}
	stepstat[5]+=1;
	//Form dW:
	for (i=0; i<neq; i++){
	  dW[neq+i] = creal(delta_w_cx[i]);
	  dW[2*neq+i] = cimag(delta_w_cx[i]);
	}
	//Test for convergence rate:
	norm_dW = error_norm(W, dW, threshold, 3*neq);
	/**
	printf("norm_dW = %g\n",norm_dW);
	printf("dW = [%g %g %g %g %g %g]\n",
	       dW[0],dW[1],dW[2],dW[3],dW[4],dW[5]);
	printf("W = [%g %g %g %g %g %g]\n",
	       W[0],W[1],W[2],W[3],W[4],W[5]);
	*/
	if (newt_iter>1){
	  theta_k = norm_dW/norm_dW_old;
	  if (newt_iter==2)
	    theta = theta_k;
	  else
	    theta = sqrt(theta_k*theta_k_old);
	  theta_k_old = theta_k;
	}
	if (theta>0.99){
	  //Diverging!
	  break;
	}
	else{
	  eta = theta/(1.0-theta);
	  conv_est = eta*norm_dW*pow(theta,newt_iter_max-newt_iter-1)/tol_newton;
	  if (conv_est>1.0){
	    //Convergence is too slow..
	    conv_cor = max(1e-4,min(20.0,conv_est));
	    abshnew = absh*0.8*pow(conv_cor,-1.0/(4.0+newt_iter_max-newt_iter-1));
	    break;
	  }
	  norm_dW_old = max(norm_dW, 100*DBL_EPSILON);
	  //Update W:
	  for (i=0; i<3*neq; i++){
	    W[i] += dW[i];
	  }
	  if (eta*norm_dW<=tol_newton){
	    //Newton has converged.
	    Newton_converged = _TRUE_;
	    break;
	  }
	}
	lasagna_test(absh<=abshmin, error_message,"Step size too small in evolver_radau5 at t=%g!. |h|=%g <= |hmin|=%g", t, absh, abshmin);
      }
      if (Newton_converged == _FALSE_){
	//Newton is too slow
	//printf("Newton failed to converge.\n");
	if (J_current == _FALSE_){
	  //Recompute jacobian and try again
	  nfenj=0;
	  lasagna_call(numjac((*derivs),
			      t,
			      y0-1,
			      f0-1,
			      &jac,
			      &nj_ws,
			      abstol,
			      neq,
			      &nfenj,
			      parameters_and_workspace_for_derivs,
			      error_message),
		       error_message,
		       error_message);
	  stepstat[3] += 1;
	  stepstat[2] += nfenj;
	  J_current = _TRUE_;
	}
	//Reduce step size:
	absh = abshnew;
	h = tdir*absh;
	//We need a new linearisation:
	new_linearisation_radau5(&jac, &jac_plus, h, neq, error_message);
        stepstat[4] +=1;
      }
      else{
	//printf("Newton converged!\n");
	//Newton converged, so we can calculate ynew and err:
	transform_C_tensor_I(T, 3, W, Y0pZ, NULL, neq);
	for (i=0; i<neq; i++){
	  ynew[i] = y0[i] + Y0pZ[2*neq+i];
	  diff[i] = f0[i] + 
	    (ei[0]*Y0pZ[i]+ei[1]*Y0pZ[neq+i]+ei[2]*Y0pZ[2*neq+i])/h;
	}
	got_ynew = _TRUE_;
	// Solve for error err:
	if (jac.use_sparse){
	  sp_lusolve(jac.Numerical, diff, err);
	}
	else{
	  for (i=0; i<neq; i++)
	    err[i] = diff[i];
	  lubksb(jac.LU,neq,jac.luidx,err-1);
	}
	//stepstat[5]+=1;
	norm_err = error_norm(ynew, err, threshold, neq);
	if ((norm_err>=rtol)&&(last_failed == _TRUE_)){
	  //Improve error estimate:
	  for (i=0; i<neq; i++){
	    ytemp[i] = err[i]+y0[i];
	  }
	  lasagna_call((*derivs)(t,
				 ytemp,
				 ftmp,
				 parameters_and_workspace_for_derivs,
				 error_message),
		       error_message,
		       error_message);
	  stepstat[2]++;
	  for (i=0; i<neq; i++){
	    diff[i] += (ftmp[i]-f0[i]);
	  }
	  //Solve for err again:
	  if (jac.use_sparse){
	    sp_lusolve(jac.Numerical, diff, err);
	  }
	  else{
	    for (i=0; i<neq; i++)
	      err[i] = diff[i];
	    lubksb(jac.LU,neq,jac.luidx,err-1);
	  }
	  //stepstat[5]+=1;
	  norm_err = error_norm(ynew, err, threshold, neq);
	}
	abshnew = safety_base*(2.0*newt_iter_max+1.0)/
	  (2.0*newt_iter_max+newt_iter)*absh*pow(rtol/norm_err,0.25);
	abshnew = min(abshnew,5.0*absh);
	abshnew = max(abshnew,0.1*absh);
	got_ynew = _TRUE_;
      }
    }
    //We have found a new y-value. Take step?
    if (norm_err>rtol){
      //Step failed
      //printf("Step failed... \n");
      stepstat[1]++;
      last_failed = _TRUE_;
      absh = max(0.1*absh,abshnew);
      h = tdir*absh;      if (J_current == _FALSE_){
	//Recalculate jacobian:
	  nfenj=0;
	  lasagna_call(numjac((*derivs),
			      t,
			      y0-1,
			      f0-1,
			      &jac,
			      &nj_ws,
			      abstol,
			      neq,
			      &nfenj,
			      parameters_and_workspace_for_derivs,
			      error_message),
		       error_message,
		       error_message);
	  stepstat[3] += 1;
	  stepstat[2] += nfenj;
	  J_current = _TRUE_;
      }
      new_linearisation_radau5(&jac, &jac_plus, h, neq, error_message);
      stepstat[4] +=1;
    }
    else{  
      //Step accepted!
      //printf("Step accepted!\n");
      stepstat[0]++;
      first_step = _FALSE_;
      last_failed = _FALSE_;
      abshlast = absh;
      //Do output?
   
      
      while((next<tres)&&((t+h-t_vec[next])*tdir >= 0.0)){
	
	dense_output_radau5(t_vec[next],
			    ytemp,
			    t,
			    h,
			    y0,
			    Y0pZ,
			    interpidx,
			    neq);
	
	//We are not interpolating the derivative at the moment..
	lasagna_call((*output)(t_vec[next],ytemp,f0,next,
			       parameters_and_workspace_for_derivs,
			     error_message),error_message,error_message);
	
	//printf("%.16e %.16e %.16e\n",t+h,ynew[0],ynew[1]);
	next++;
      }
      //Update parameters:
      t += h;
      for (i=0; i<neq; i++){
	ylast[i] = y0[i];
	y0[i] = ynew[i];
      }
      for (i=0; i<3*neq; i++){
	Zlast[i] = Y0pZ[i];
      }
      lasagna_call((*derivs)(t,
			     y0,
			     f0,
			     parameters_and_workspace_for_derivs,error_message),
		   error_message,
		   error_message);
      stepstat[2]++;
      
      //Decide how to proceed in the next step:
      reuse_jacobian = _FALSE_;
      reuse_stepsize = _FALSE_;
      if ((newt_iter == 1)||(theta < theta_reuse_jacobian)){
	//printf("Re-using jacobian..\n");
	//Reuse jacobian.
	//printf("Reusing jacobian...\n");
	reuse_jacobian = _TRUE_;
	if ((abshnew/absh >= same_step_lower)&&
	    (abshnew/absh <= same_step_upper)){
	  reuse_stepsize = _TRUE_;
	}
      }
      if (reuse_jacobian == _TRUE_){
	J_current = _FALSE_;
      }
      else{
	//Recalculate jacobian:
	nfenj=0;
	lasagna_call(numjac((*derivs),
			    t,
			    y0-1,
			    f0-1,
			    &jac,
			    &nj_ws,
			    abstol,
			    neq,
			    &nfenj,
			    parameters_and_workspace_for_derivs,
			    error_message),
		     error_message,
		     error_message);
	stepstat[3] += 1;
	stepstat[2] += nfenj;
	J_current = _TRUE_;
      }
      if (abshnew>fabs(t_final-t)){
	abshnew = fabs(t_final-t);
	reuse_stepsize = _FALSE_;
      }
      if (reuse_stepsize == _FALSE_){
	//Change step size and do a new linearisation:
	absh = abshnew;
	h = tdir*absh;
	new_linearisation_radau5(&jac, &jac_plus, h, neq, error_message);
        stepstat[4] +=1;
      }
    }
    if (fabs(y0[0]*1e-15)>5e-6){
      //Stop condition
      lasagna_call((*output)(t,y0,f0,next,
			     parameters_and_workspace_for_derivs,
			     error_message),error_message,error_message);
      printf("Stop condition met...\n");
      break;
    }
  }
  printf("\n End of evolver. Next=%d, t=%e and tnew=%e.",next,t,t+h);
  printf("\n Statistics: [%d %d %d %d %d %d] \n",stepstat[0],stepstat[1],
	 stepstat[2],stepstat[3],stepstat[4],stepstat[5]);
	
  /** Deallocate memory */
  uninitialize_jacobian_plus(&jac_plus);
  uninitialize_jacobian(&jac);
  uninitialize_numjac_workspace(&nj_ws);

  free(W);
  free(dW_buf);
  free(Y0pZ);
  free(Fi);
  free(Zlast);
  free(rhs_buf);
  free(err_buf);
  free(diff);
  free(xtemp);
  free(ytemp);
  free(ynew);
  free(f0);
  free(ylast);
  free(ftmp);
  free(dfdt);
  free(delta_w);
  free(rhs_cx);
  free(delta_w_cx_buf);


  return _SUCCESS_;
}

int initialize_jacobian_plus(struct jacobian *jac, 
			     struct jacobian_plus *jac_plus,
			     int neq,
			     ErrorMsg error_message){
  int i,neqp1 = neq+1;

  lasagna_alloc(jac_plus->LU_cx,
		sizeof(double complex*)*neqp1,
		error_message); 
  /* Allocate vector of pointers to rows of matrix.*/
  lasagna_alloc(jac_plus->LU_cx[1],sizeof(double complex)*neqp1*neq,
		error_message);
  jac_plus->LU_cx[0] = NULL;
  for(i=2;i<=neq;i++) 
    jac_plus->LU_cx[i] = jac_plus->LU_cx[i-1]+neq; /* Set row pointers... */ 
	
  lasagna_alloc(jac_plus->luidx,sizeof(int)*neqp1,error_message);

  if (jac->use_sparse){
    lasagna_call(sp_num_alloc_cx(&jac_plus->Numerical_cx, neq,error_message),
		 error_message,error_message);
		
    lasagna_call(sp_mat_alloc_cx(&jac_plus->spJ_cx, neq, neq, jac->max_nonzero,
			      error_message),error_message,error_message);

    /** Immediately free spJ_cx->Ai, spJ_cx->Ap and Numerical_cx->q, 
	and set pointers to the corresponding locations in jac:
    */
    free(jac_plus->spJ_cx->Ai);
    free(jac_plus->spJ_cx->Ap);
    free(jac_plus->Numerical_cx->q);
    jac_plus->spJ_cx->Ai = jac->spJ->Ai;
    jac_plus->spJ_cx->Ap = jac->spJ->Ap;
    jac_plus->Numerical_cx->q = jac->Numerical->q;
    jac_plus->sparse_stuff_initialised = _TRUE_;
  }
  else{
    jac_plus->sparse_stuff_initialised = _FALSE_;
  }
  return _SUCCESS_;
}

int uninitialize_jacobian_plus(struct jacobian_plus *jac_plus){
  free(jac_plus->LU_cx[1]);
  free(jac_plus->LU_cx);
  free(jac_plus->luidx);
  
  if (jac_plus->sparse_stuff_initialised == _TRUE_){
    //Free rest of spJ_cx manually:
    free(jac_plus->spJ_cx->Ax);
    free(jac_plus->spJ_cx);
    //The column ordering is not allocated but is a copy of the pointer in
    //jac, so we are not allowed to free it here:
    jac_plus->Numerical_cx->q = NULL;
    sp_num_free_cx(jac_plus->Numerical_cx);
  }
  return _SUCCESS_;
}


int new_linearisation_radau5(struct jacobian *jac, 
			     struct jacobian_plus *jac_plus,
			     double hnew,
			     int neq,
			     ErrorMsg error_message){

  double gamma_hat = 3.0-pow(3.0,1.0/3.0)+pow(3.0,2.0/3.0);
  double alpha_hat = 0.5*(6.0+pow(3.0,1.0/3.0)-pow(3.0,2.0/3.0));
  double beta_hat = 0.5*pow(3.0,1.0/6.0)*(3.0+pow(3.0,2.0/3.0));
  double gamma = gamma_hat/hnew;
  double complex alpha_ibeta = alpha_hat/hnew+I*(beta_hat/hnew);

  double luparity, *Ax;
  double complex *Az;
  int i,j,*Ap,*Ai,funcreturn;
  if(jac->use_sparse==1){
    Ap = jac->spJ->Ap; Ai = jac->spJ->Ai; Ax = jac->spJ->Ax;
    Az = jac_plus->spJ_cx->Ax;
    /* Construct jac->spJ->Ax from jac->xjac, the jacobian:*/
    for(j=0;j<neq;j++){
      for(i=Ap[j];i<Ap[j+1];i++){
	if(Ai[i]==j){
	  /* I'm at the diagonal */
	  Ax[i] = gamma-jac->xjac[i];
	  Az[i] = alpha_ibeta-jac->xjac[i];
	}
	else{
	  Ax[i] = -jac->xjac[i];
	  Az[i] = -jac->xjac[i];
	}
      }
    }
    /* Matrix constructed... */
    if (jac->pattern_supplied == _TRUE_){
      //      if (jac->refactor_count%jac->refactor_max == 0){
      if (jac->new_jacobian == _TRUE_){
	//printf("Full calculation...\n");
	//We should do a full LU decomposition again:
	funcreturn = sp_ludcmp(jac->Numerical, jac->spJ, 1e-3);
	lasagna_test(funcreturn == _FAILURE_,error_message,
		     "Failure in sp_ludcmp. Possibly singular matrix!");
	/** Column ordering in spJ_cx->Numerical is (should be) pointing 
	    to column_ordering in spJ->Numerical.
	    Same with Ai and Ap.
	*/
	funcreturn = sp_ludcmp_cx(jac_plus->Numerical_cx, jac_plus->spJ_cx, 1e-3);
	lasagna_test(funcreturn == _FAILURE_,error_message,
		     "Failure in sp_ludcmp_cx. Possibly singular matrix!");
	jac->new_jacobian = _FALSE_;
      }
      else{
	//printf("Refactorisation..\n");
	sp_refactor(jac->Numerical, jac->spJ);
	sp_refactor_cx(jac_plus->Numerical_cx, jac_plus->spJ_cx);
      }
      jac->refactor_count++;
    }
    else{
      if((jac->new_jacobian==_TRUE_)&&(jac->repeated_pattern<1)){
        /*I have a new pattern, and I have not done a LU decomposition 
	  since the last jacobian calculation, so	I need to do a full 
	  sparse LU-decomposition: */
	/* Find the sparsity pattern C = J + J':*/
	calc_C(jac);
	/* Calculate the optimal ordering: */
	sp_amd(jac->Cp, jac->Ci, neq, jac->cnzmax,
	       jac->Numerical->q,jac->Numerical->wamd);
	/* if the next line is uncomented, the code uses natural ordering instead of AMD ordering */
	/*jac->Numerical->q = NULL;*/
	funcreturn = sp_ludcmp(jac->Numerical, jac->spJ, 1e-3);
	lasagna_test(funcreturn == _FAILURE_,error_message,
		     "Failure in sp_ludcmp. Possibly singular matrix!");
	/** Column ordering in spJ_cx->Numerical is (should be) pointing 
	    to column_ordering in spJ->Numerical.
	    Same with Ai and Ap.
	*/
	funcreturn = sp_ludcmp_cx(jac_plus->Numerical_cx, jac_plus->spJ_cx, 1e-3);
	lasagna_test(funcreturn == _FAILURE_,error_message,
		     "Failure in sp_ludcmp_cx. Possibly singular matrix!");

	/**printf("Non-zero elements in real LU-decomposition: %d+%d. \nIn complex LU: %d+%d.\n The base matrix is [%d x %d]=%d\n",jac->Numerical->L->Ap[neq],jac->Numerical->U->Ap[neq],jac_plus->Numerical_cx->L->Ap[neq],jac_plus->Numerical_cx->U->Ap[neq],neq,neq,neq*neq);
	 */
	jac->new_jacobian = _FALSE_;
      }
      else{
	/* I have a repeated pattern, so I can just refactor:*/
	sp_refactor(jac->Numerical, jac->spJ);
	sp_refactor_cx(jac_plus->Numerical_cx, jac_plus->spJ_cx);
      }
    }
  }
  else{
    /* Normal calculation: */
    for(i=1;i<=neq;i++){
      for(j=1;j<=neq;j++){
	jac->LU[i][j] = -jac->dfdy[i][j];
	jac_plus->LU_cx[i][j] = -jac->dfdy[i][j];
	if(i==j){
	  jac->LU[i][j] +=gamma;
	  jac_plus->LU_cx[i][j] +=alpha_ibeta;
	}
      }
    }
    /*Dense LU decomposition: */
    funcreturn = ludcmp(jac->LU,neq,jac->luidx,&luparity,jac->LUw);
    lasagna_test(funcreturn == _FAILURE_,error_message,
	       "Failure in ludcmp. Possibly singular matrix!");
    funcreturn = ludcmp_cx(jac_plus->LU_cx,neq,jac_plus->luidx,&luparity,jac->LUw);
    lasagna_test(funcreturn == _FAILURE_,error_message,
	       "Failure in ludcmp. Possibly singular matrix!");
  }
  return _SUCCESS_;
}

int transform_C_tensor_I(double *C, 
			 int s,
			 double *vec_in,
			 double *vec_out,
			 double *vec_init,
			 int n){

  /** Apply C\otimesI on a vector. If vec_init is non-null, we
      form vec_out = vec_init x I + CxI*vec_in.
  */
  int i, j, k;
  double w;
  
  //Set initial value for vector:
  if (vec_init==NULL){
    for (i=0; i<s*n; i++){
      vec_out[i] = 0.0;
    }
  }
  else{
    for (j=0; j<s; j++){
      for (i=0; i<n; i++){
	vec_out[j*n+i] = vec_init[i];
      }
    }
  }
  //Do the transformation:
  for (i=0; i<n; i++){
    for (j=0; j<s; j++){
      w = vec_in[j*n+i];
      for (k=0; k<s; k++){
	vec_out[k*n+i] += w*C[k*s+j];
      }
    }
  }
  return _SUCCESS_;
}
 



int dense_output_radau5(double tinterp, 
			double *yi, 
			double t0, 
			double h,
			double *y0, 
			double *Z,
			int *interpidx,
			int neq){
  double theta, theta2, theta3;
  int i;
  double sr6;

  theta = (tinterp-t0)/h;
  theta2 = theta*theta;
  theta3 = theta*theta2;
  sr6 = sqrt(6.0);

  if (y0 != NULL){
    for (i=0; i<neq; i++){
      if (interpidx != NULL){
	if (interpidx[i] == 0){
	  continue;
	}
      }
      yi[i] = y0[i] + 
	((13.0+7.0*sr6)*theta/3.0+
	 (-23.0-22.0*sr6)*theta2/3.0+
	 (10.0/3.0+5.0*sr6)*theta3)*Z[i] +
	((13.0-7.0*sr6)*theta/3.0+
	 (-23.0+22.0*sr6)*theta2/3.0+
	 (10.0/3.0-5.0*sr6)*theta3)*Z[neq+i] +
	(theta/3.0-8.0*theta2/3.0+10.0*theta3/3.0)*Z[2*neq+i]; 
    }
  }
  else{
    for (i=0; i<neq; i++){
      if (interpidx != NULL){
	if (interpidx[i] == 0){
	  continue;
	}
      }
      yi[i] =
	((13.0+7.0*sr6)*theta/3.0+
	 (-23.0-22.0*sr6)*theta2/3.0+
	 (10.0/3.0+5.0*sr6)*theta3)*Z[i] +
	((13.0-7.0*sr6)*theta/3.0+
	 (-23.0+22.0*sr6)*theta2/3.0+
	 (10.0/3.0-5.0*sr6)*theta3)*Z[neq+i] +
	(theta/3.0-8.0*theta2/3.0+10.0*theta3/3.0)*Z[2*neq+i]; 
    }
  }
  return _SUCCESS_;
}
       
double norm_inf(double *y, 
		double *err_y, 
		double threshold,
		int neq){
  int i;
  double wt,max_tmp=0.0;
  
  for (i=0; i<neq; i++){
    wt = max(threshold, fabs(y[i]));
    max_tmp = max(max_tmp,fabs(err_y[i]/wt));
    //    printf("wt=%g, max_tmp=%g\n",wt,max_tmp);
  }
  return max_tmp;
}

double norm_L2(double *y, 
	       double *err_y,
	       double threshold,
	       int neq){
  int i;
  double t, wt, sum = 0.0;
  
  for (i=0; i<neq; i++){
    wt = max(threshold, fabs(y[i]));
    t = fabs(err_y[i]/wt);
    sum += t*t;
  }
  return sqrt(sum/((double) neq));
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
