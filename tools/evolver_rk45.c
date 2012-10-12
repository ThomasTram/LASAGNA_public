#include "evolver_rk45.h"
int evolver_rk45(int (*derivs)(double x,double * y,double * dy,
			       void * parameters_and_workspace, ErrorMsg error_message),
		 void * ppaw,
		 double x_ini,
		 double x_final,
		 double * y_inout, 
		 int neq, 
		 EvolverOptions *options,
		 ErrorMsg error_message){
	
  /** Handle options: */
  int *interpidx, *stats, verbose, t_res; 
  double abstol, rtol, *t_vec;
  int (*output)(double t, double *y, double *dy, int i, void *p, ErrorMsg err);
  int (*print_variables)(double t, double *y, double *dy, void *p, ErrorMsg err);
  interpidx = options->used_in_output; stats = &(options->Stats[0]);
  verbose = options->EvolverVerbose; t_res = options->tres; abstol = options->AbsTol; 
  rtol = options->RelTol; t_vec = options->t_vec; 
  output = options->output; print_variables=options->print_variables; 
 
  double *dy,*err,*ynew,*ytemp, *ki;
  double h,errmax,errtemp,hmin,hnew;
  double t;
  int tdir, i, j, k, idx=0;
  double ci[6];
  double bi[6];
  double bi_diff[6];
  double ai[6][6];
  int s=0;
  double threshold = 1e-6/rtol;
  int nofailed;
  double pow_grow=0.5;
  double rh,maxtmp,absh,hmax;

  dy = malloc(sizeof(double)*neq);
  ytemp = malloc(sizeof(double)*neq);
  ynew = malloc(sizeof(double)*neq);
  err = malloc(sizeof(double)*neq);
  if (_FALSE_){
    //Runge-Kutta Fehlberg method:
    s = 6;
    pow_grow = 0.2;
    ci[0] = 0.0;
    ci[1] = 0.25;
    ci[2] = 3.0/8.0;
    ci[3] = 12.0/13.0;
    ci[4] = 1.0;
    ci[5] = 0.5;
    bi[0] = 16.0/135.0;
    bi[1] = 0.0;
    bi[2] = 6656.0/12825.0;
    bi[3] = 28561.0/56430.0;
    bi[4] = -9.0/50.0;
    bi[5] = 2.0/55.0;
    bi_diff[0] = 25.0/216.0-16.0/135.0;
    bi_diff[1] = 0.0;
    bi_diff[2] = 1408.0/2565.0-6656.0/12825.0;
    bi_diff[3] = 2197.0/4104.0-28561.0/56430.0;
    bi_diff[4] = -0.2+9.0/50;
    bi_diff[5] = -2.0/55.0;
    ai[1][0] = 0.25;
    ai[2][0] = 3.0/32.0;
    ai[2][1] = 9.0/32.0;
    ai[3][0] = 1932.0/2197.0;
    ai[3][1] = -7200.0/2197.0;
    ai[3][2] = 7296.0/2197.0;
    ai[4][0] = 439.0/216.0;
    ai[4][1] = -8.0;
    ai[4][2] = 3680.0/513.0;
    ai[4][3] = -845.0/4104.0;
    ai[5][0] = -8.0/27.0;
    ai[5][1] = 2.0;
    ai[5][2] = -3544.0/2565.0;
    ai[5][3] = 1859.0/4104.0;
    ai[5][4] = -11.0/40.0;
  }
  else{
    //Cash-Karp Runge-Kutta method
    s = 6;
    pow_grow = 0.2;
    ci[0] = 0.0;
    ci[1] = 0.2;
    ci[2] = 0.3;
    ci[3] = 0.6;
    ci[4] = 1.0;
    ci[5] = 7.0/8.0;
    bi[0] = 37.0/378.0;
    bi[1] = 0.0;
    bi[2] = 250.0/621.0;
    bi[3] = 125.0/594.0;
    bi[4] = 0.0;
    bi[5] = 512.0/1771.0;
    bi_diff[0] = 37.0/378.0-2825.0/27648.0;
    bi_diff[1] = 0.0;
    bi_diff[2] = 250.0/621.0-18575.0/48384.0;
    bi_diff[3] = 125.0/594.0-13525.0/55296.0;
    bi_diff[4] = -277.0/14336.0;
    bi_diff[5] = 512.0/1771.0-0.25;
    ai[1][0] = 0.2;
    ai[2][0] = 3.0/40.0;
    ai[2][1] = 9.0/40.0;
    ai[3][0] = 0.3;
    ai[3][1] = -0.9;
    ai[3][2] = 1.2;
    ai[4][0] = -11.0/54.0;
    ai[4][1] = 2.5;
    ai[4][2] = -70.0/27.0;
    ai[4][3] = 35.0/27.0;
    ai[5][0] = 1631.0/55296.0;
    ai[5][1] = 175.0/512.0;
    ai[5][2] = 575.0/13824.0;
    ai[5][3] = 44275.0/110592.0;
    ai[5][4] = 253.0/4096.0;
  }
  ki = malloc(sizeof(double)*s*neq);
  t = x_ini;
  derivs(t,y_inout,dy,ppaw, error_message);
  stats[2]++;
  output(t,y_inout,dy,0,ppaw, error_message);
  //
  hmin = 100.0*DBL_MIN*fabs(t);
  hmax = fabs(x_final-x_ini)/10.0;
  absh = hmax;
  //Compute h_initial:
  for (k=0,rh=0.0; k<neq; k++){
    maxtmp = max(fabs(y_inout[k]),threshold);
    rh = max(rh,fabs(dy[k])/maxtmp);
  }
  rh /= 0.8*pow(rtol,pow_grow);
  if (absh*rh>1.0)
    absh = 1.0/rh;
  if (x_final>x_ini)
    tdir = 1;
  else
    tdir = -1;
  hnew = absh*tdir;

  //
  h = hnew;
  nofailed = _TRUE_;

  for (idx=1; idx<t_res; idx++){
    while ((t-t_vec[idx])*tdir<0.0){
      hmin = 100.0*DBL_MIN*fabs(t);
      h = hnew;
      //printf("Estimated steps remaining: %g\n",fabs((x_final-x_ini)/h));
      if (fabs(h)>fabs(t_vec[idx]-t)){
	h = t_vec[idx]-t;
      }
      if (fabs(h)<hmin)
	return _FAILURE_;
      //Try to take a step h:
      for (k=0; k<neq; k++){
	ynew[k] = y_inout[k];
	err[k] = 0.0;
      }
      for (i=0; i<s; i++){
	//Reset ytemp.
	for (k=0; k<neq; k++)
	  ytemp[k] = y_inout[k];
	if (verbose>2)
	  printf("i: %d\n",i);
	for (j=0; j<i; j++){
	  if (verbose >2)
	    printf("j: %d, a[i][j] = %g\n",j,ai[i][j]);
	  for (k=0; k<neq; k++){
	    ytemp[k] += h*ai[i][j]*ki[j*neq+k];
	  }
	}
	// Calculate k_i:
	if (verbose>1){
	  printf("Evaluating ODE at t=%g, ci[i] = %g. h=%g\n",
		 t+ci[i]*h,ci[i],h);
	  printf("y_inout = [%g,%g]. ytemp = [%g,%g]\n",
		 y_inout[0],y_inout[1],ytemp[0],ytemp[1]);
	}
	derivs(t+ci[i]*h,ytemp,ki+i*neq,ppaw, error_message);
	stats[2]++;
	// Update ynew and err:
	for (k=0; k<neq; k++){  
	  ynew[k] += h*bi[i]*ki[i*neq+k];
	  err[k] += h*bi_diff[i]*ki[i*neq+k];
	}
      }
      if (verbose>1)
	printf("Finished loop over i, new y has been found.\n");
      // Got new y and error estimate.
      for (k=0,errmax = 0.0; k<neq; k++){
	errtemp = fabs(err[k]/max(threshold,fabs(ynew[k])));
	if (errtemp>errmax)
	  errmax = errtemp;
      }
      //printf("h: %g, errmax = %g\n",h,errmax);
      if (errmax>rtol){
	if (verbose>1) printf("Step rejected..\n");
	stats[1]++;
	//Step rejected.
	if (nofailed == _TRUE_){
	  hnew = tdir*max(hmin, fabs(h) * max(0.1, 0.8*pow(rtol/errmax,pow_grow)));
	}
	else{
	  hnew = tdir*max(hmin,0.5*fabs(h));
	}
	nofailed = _FALSE_;
      }
      else{
	//Step accepted.
	stats[0]++;
	if (verbose>1)
	  printf("Step accepted..\n");
	nofailed = _TRUE_;
	hnew = tdir*max(hmin, fabs(h) * max(0.1, 0.8*pow(rtol/errmax,pow_grow)));
	for (k=0; k<neq; k++){
	  y_inout[k] = ynew[k];
	}
	t += h;
      }
    }
    // Store values at this point:
    derivs(t,y_inout,dy,ppaw, error_message);
    output(t,y_inout,dy,idx,ppaw, error_message);
  }
  if (verbose>0)
    printf(" Successful steps: %d\n Failed steps: %d\n Function evaluations: %d\n",
	   stats[0],stats[1],stats[2]);
  free(dy);
  free(ynew);
  free(ytemp);
  free(ki);
  free(err);
  return _SUCCESS_;
}

int evolver_rkdp45(int (*derivs)(double x,
				 double * y,
				 double * dy,
				 void * ppaw, 
				 ErrorMsg error_message),
		   void *ppaw,
		   double t_ini,
		   double t_final,
		   double * y_inout, 
		   int neq, 
		   EvolverOptions *options,
		   ErrorMsg error_message){
	
  /** Handle options: */
  int *used_in_output, *stats, verbose, tres; 
  double abstol, rtol, *t_vec;
  int (*output)(double t, double *y, double *dy, int i, void *p, ErrorMsg err);
  int (*print_variables)(double t, double *y, double *dy, void *p, ErrorMsg err);
  int (*stop_function)(double t, double *y, double *dy, void *p, ErrorMsg err);
  used_in_output = options->used_in_output; stats = &(options->Stats[0]);
  verbose = options->EvolverVerbose; tres = options->tres; abstol = options->AbsTol; 
  rtol = options->RelTol; t_vec = options->t_vec; 
  output = options->output; print_variables=options->print_variables; 
  stop_function = options->stop_function;

  double *dy,*err,*ynew,*ytemp, *ki;
  double h,absh,hmax,errmax,errtemp,hmin,hnew;
  double t,tnew;
  int tdir, i, j, k, idx=0;
  int s=7;
  double ci[s];
  double bi[s];
  double bi_diff[s];
  double ai[s][s];
  double bi_vec_y[s],bi_vec_dy[s];
  double threshold = abstol/rtol;
  int nofailed;
  double pow_grow=0.2;
  double rh,maxtmp;
  //Interpolation variables:
  double i01,i02,i03;
  double ixx[5][3];
  double ti,ss1,ss2,ss3,ss4,*yinterp,*dyinterp;
  dy = malloc(sizeof(double)*neq);
  ytemp = malloc(sizeof(double)*neq);
  yinterp = malloc(sizeof(double)*neq);
  dyinterp = malloc(sizeof(double)*neq);
  ynew = malloc(sizeof(double)*neq);
  err = malloc(sizeof(double)*neq);
  pow_grow = 0.2;
  ki = malloc(sizeof(double)*s*neq);  
  /** Set method parameters for Runge-Kutta Dormand-Prince method:
      ------------------------------------------------------------------
  */
  ci[0] = 0.0; ci[1] = 0.2; ci[2] = 0.3; ci[3] = 0.8; 
  ci[4] = 8.0/9.0; ci[5] = 1.0; ci[6] = 1.0;
  bi[0] = 35.0/384.0; bi[1] = 0.0; bi[2] = 500.0/1113.0; 
  bi[3] = 125.0/192.0; bi[4] = -2187.0/6784.0; bi[5] = 11.0/84.0; bi[6] = 0.0; 
  bi_diff[0] = 71.0/57600.0; bi_diff[1] = 0.0; bi_diff[2] = -71.0/16695.0;
  bi_diff[3] = 71.0/1920.0; bi_diff[4] = -17253.0/339200.0; 
  bi_diff[5] = 22.0/525.0; bi_diff[6] = -1.0/40.0;
  ai[1][0] = 0.2;
  ai[2][0] = 3.0/40.0; ai[2][1] = 9.0/40.0;
  ai[3][0] = 44.0/45.0; ai[3][1] = -56.0/15.0; ai[3][2] = 32.0/9.0;
  ai[4][0] = 19372.0/6561.0; ai[4][1] = -25360.0/2187.0;
  ai[4][2] = 64448.0/6561.0; ai[4][3] = -212.0/729.0;
  ai[5][0] = 9017.0/3168.0; ai[5][1] = -355.0/33.0; ai[5][2] = 46732.0/5247.0;
  ai[5][3] = 49.0/176.0; ai[5][4] = -5103.0/18656.0;
  ai[6][0] = 35.0/384.0; ai[6][1] = 0.0; ai[6][2] = 500.0/1113.0;
  ai[6][3] = 125.0/192.0; ai[6][4] = -2187.0/6784.0; ai[6][5] = 11.0/84.0;
  ixx[0][0] = 1500.0/371.0; ixx[0][1] = -1000.0/159.0; ixx[0][2] = 1000.0/371.0;
  ixx[1][0] = -125.0/32.0; ixx[1][1] = 125.0/12.0; ixx[1][2] = -375.0/64.0;
  ixx[2][0] = 9477.0/3392.0; ixx[2][1] = -729.0/106.0; ixx[2][2] = 25515.0/6784.0;
  ixx[3][0] = -11.0/7.0; ixx[3][1] = 11.0/3.0; ixx[3][2] = -55.0/28.0;
  ixx[4][0] = 1.5; ixx[4][1] = -4.0; ixx[4][2] = 2.5;
  i01 = -183.0/64.0; i02=37.0/12.0; i03 = -145.0/128.0;
  /**
     Done setting method parameters.
  */

  t = t_ini;
  //initialise ki
  derivs(t,y_inout,ki,ppaw, error_message);
  stats[2]++;
  hmin = 100.0*DBL_MIN*fabs(t);
  hmax = fabs(t_final-t_ini)/10.0;
  if (t_vec!=NULL)
    absh = min(hmax,fabs(t_vec[1]-t_vec[0]));
  else
    absh = hmax;
  if (absh==0.0)
    absh = hmax;
  //Compute h_initial:
  for (k=0,rh=0.0; k<neq; k++){
    maxtmp = max(fabs(y_inout[k]),threshold);
    rh = max(rh,fabs(ki[k])/maxtmp);
  }
  rh /= 0.8*pow(rtol,pow_grow);
  if (absh*rh>1.0)
    absh = 1.0/rh;
  if (t_final>t_ini)
    tdir = 1;
  else
    tdir = -1;

  hnew = absh*tdir;
  h = hnew;
  //Find current index:
  if (t_vec != NULL){
    //Output at specified points.
    for(idx = 0; (t_vec[idx]-t)*tdir<0.0; idx++);
  }
  nofailed = _TRUE_;
  
  while ((t-t_final)*tdir<0.0){
    h = hnew;
    hmin = 100.0*DBL_MIN*fabs(t);
    if (fabs(h)<hmin){
      printf("h = %e",fabs(h));
      return _FAILURE_;
    }
    if (fabs(h)>0.9*fabs(t_final-t))
      h = t_final-t;
    //Try to take a step h.
    for (k=0; k<neq; k++){
      //Set ynew and err to the value after i=0 iteration:
      ynew[k] = y_inout[k]+h*bi[0]*ki[k];
      err[k] = h*bi_diff[0]*ki[k];
    }
    for (i=1; i<s; i++){
      //Note: starting at i=1 since method has FSAL property.
      //Reset ytemp. 
      for (k=0; k<neq; k++)
	ytemp[k] = y_inout[k];
      for (j=0; j<i; j++){
	for (k=0; k<neq; k++){
	  ytemp[k] += h*ai[i][j]*ki[j*neq+k];
	}
      }
      // Calculate k_i:
      if (verbose>2){
	printf("Evaluating ODE at t=%g, ci[i] = %g. h=%g\n",
	       t+ci[i]*h,ci[i],h);
	printf("y_inout = [%g,%g]. ytemp = [%g,%g]\n",
	       y_inout[0],y_inout[1],ytemp[0],ytemp[1]);
      }
      derivs(t+ci[i]*h,ytemp,ki+i*neq,ppaw, error_message);
      stats[2]++;
      // Update ynew and err:
      for (k=0; k<neq; k++){  
	ynew[k] += h*bi[i]*ki[i*neq+k];
	err[k] += h*bi_diff[i]*ki[i*neq+k];
      }
    }
    if (verbose>3)
      printf("Finished loop over i, new y has been found.\n");
    // Got new y and error estimate.
    for (k=0,errmax = 0.0; k<neq; k++){
      errtemp = fabs(err[k]/max(threshold,fabs(ynew[k])));
      if (errtemp>errmax){
	errmax = errtemp;
      }
    }
    if (verbose>3)
      printf("h: %g, errmax = %g\n",h,errmax);
    if (errmax>rtol){
      if (verbose>4) 
	printf("Step rejected..\n");
      stats[1]++;
      //Step rejected.
      if (nofailed == _TRUE_){
	hnew = tdir*max(hmin, fabs(h) * max(0.1, 0.8*pow(rtol/errmax,pow_grow)));
      }
      else{
	hnew = tdir*max(hmin,0.5*fabs(h));
	nofailed = _FALSE_;
      }
    }
    else{
      //Step accepted.
      stats[0]++;
      if (print_variables!=NULL){
	print_variables(t+h,ynew,ki+6*neq,ppaw,error_message); 
      }
      if (verbose>1)
	printf("Step accepted. t=%g, h=%g\n",t,h);
      nofailed = _TRUE_;
      
      hnew = tdir*max(hmin, fabs(h) * max(0.1, 0.8*pow(rtol/errmax,pow_grow)));
      tnew = t+h;
      //Do we need to write output?
      if (t_vec==NULL){
	//Refined output
	for (idx=1; idx<tres; idx++){
	  ti = t+idx/((double) tres)*h;
	  ss1 = (ti-t)/h; ss2=ss1*ss1; ss3=ss2*ss1; ss4=ss2*ss2;
	  bi_vec_y[0] = ss1+i01*ss2+i02*ss3+i03*ss4;
	  bi_vec_dy[0] = 1.0+i01*2.0*ss1+i02*3.0*ss2+i03*4.0*ss3;
	  //bi_vec_y[1] = 0.0; bi_vec_dy[1] = 0.0;
	  for (i=2; i<7; i++){
	    bi_vec_y[i] = ixx[i-2][0]*ss2+ixx[i-2][1]*ss3+ixx[i-2][2]*ss4;
	    bi_vec_dy[i] = ixx[i-2][0]*2*ss1+ixx[i-2][1]*3*ss2+ixx[i-2][2]*4*ss3;
	  }
	  for (k=0; k<neq; k++){
	    if (used_in_output[k] == _TRUE_){
	      //Interpolate
	      yinterp[k] = y_inout[k];
	      dyinterp[k] = 0.0;
	      for (i=0; i<7; i++){
		if (i!=1){
		  yinterp[k] +=h*bi_vec_y[i]*ki[i*neq+k];
		  dyinterp[k] +=bi_vec_dy[i]*ki[i*neq+k];
		}
	      }
	    }
	  }
	  output(ti,yinterp,dyinterp,idx,ppaw, error_message);
	}
	output(tnew,ynew,ki+6*neq,tres,ppaw, error_message);
      }
      else{
	for(; (idx<tres)&&((tnew-t_vec[idx])*tdir>=0.0); idx++){
	  if (tnew==t_vec[idx]){
	    //We have hit the point exactly. Use ynew and dy=ki+6*neq
	    output(tnew,ynew,ki+6*neq,idx,ppaw, error_message);
	  }
	  else{
	    //Interpolate to get output using the information in the ki-matrix:
	    ti = t_vec[idx];
	    ss1 = (ti-t)/h; ss2=ss1*ss1; ss3=ss2*ss1; ss4=ss2*ss2;
	    bi_vec_y[0] = ss1+i01*ss2+i02*ss3+i03*ss4;
	    bi_vec_dy[0] = 1.0+i01*2.0*ss1+i02*3.0*ss2+i03*4.0*ss3;
	    //bi_vec_y[1] = 0.0; bi_vec_dy[1] = 0.0;
	    for (i=2; i<7; i++){
	      bi_vec_y[i] = ixx[i-2][0]*ss2+ixx[i-2][1]*ss3+ixx[i-2][2]*ss4;
	      bi_vec_dy[i] = ixx[i-2][0]*2*ss1+ixx[i-2][1]*3*ss2+ixx[i-2][2]*4*ss3;
	    }
	    for (k=0; k<neq; k++){
	      if (used_in_output[k] == _TRUE_){
		//Interpolate
		yinterp[k] = y_inout[k];
		dyinterp[k] = 0.0;
		for (i=0; i<7; i++){
		  if (i!=1){
		    yinterp[k] +=h*bi_vec_y[i]*ki[i*neq+k];
		    dyinterp[k] +=bi_vec_dy[i]*ki[i*neq+k];
		  }
		}
	      }
	    }
	    output(ti,yinterp,dyinterp,idx,ppaw, error_message);
	  }
	}
      }
      /* Perhaps use stop function: */
      if (stop_function != NULL){
	if (stop_function(tnew,ynew,ki+6*neq, ppaw, error_message)==_TRUE_){
	  output(tnew,ynew,ki+6*neq,idx,ppaw, error_message);
	  printf("Stop condition met...\n");
	  break;
	}
      }
      for (k=0; k<neq; k++){
	//Update y:
	y_inout[k] = ynew[k];
	//Update k0:
	ki[k] = ki[6*neq+k];
      }
      t = tnew;
    }
  }
  if (verbose>0)
    printf(" Successful steps: %d\n Failed steps: %d\n Function evaluations: %d\n",
	   stats[0],stats[1],stats[2]);
  free(dy);
  free(ynew);
  free(ytemp);
  free(ki);
  free(err);
  return _SUCCESS_;
}
