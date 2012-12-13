#include "newton.h"

int Newton(int (*vecfun)(double * y, double * Fy, void *param),
	   int (*jacfun)(double * y, double ** jac, void *param),
	   double *y0,
	   void *param,
	   double *maxstep,
	   double rtol,
	   int *iter,
	   int maxiter,
	   size_t neq,
	   ErrorMsg error_message){
  double *Fval,**jac,lusign,*vv,*mdelta_y;
  int i,*indx,err_idx=-1;
  int converged=_TRUE_;
  double reldif;
  double abstol=1e-15;

  Fval = malloc(sizeof(double)*neq);
  vv = malloc(sizeof(double)*(neq+1));
  jac = malloc(sizeof(double*)*(neq+1));
  indx = malloc(sizeof(int)*(neq+1));
  mdelta_y = malloc(sizeof(double)*(neq+1));
  
  for (i=1; i<=neq; i++) 
    jac[i] = malloc(sizeof(double)*(neq+1));

  for (*iter=1; *iter<=maxiter; (*iter)++){
    //Do Newton iteration:
    converged = _TRUE_;
    vecfun(y0,Fval,param);
    if (jacfun == NULL){
      jacobian_for_Newton(vecfun,y0,Fval,neq,param,jac);
    }
    else{
      jacfun(y0,jac,param);
    }

    /**
       lasagna_call(ludcmp(jac, neq, indx, &lusign, vv),
		 "Problem in ludcmp",error_message);
    */
    if (ludcmp(jac, neq, indx, &lusign, vv) == _FAILURE_){
      jacobian_for_Newton(vecfun,y0,Fval,neq,param,jac);
      printf("Fval = [%g %g %g]\n",Fval[0],Fval[1],Fval[2]);
      printf("y0 = [%g %g %g]\n",y0[0],y0[1],y0[2]);
      for (i=1; i<=neq; i++){
	printf("|%g %g %g|\n",jac[i][1],jac[i][2],jac[i][3]);
      }
    }


    for(i=0; i<neq; i++) 
      mdelta_y[i+1] = Fval[i];
    
    lasagna_call(lubksb(jac, neq, indx, mdelta_y),
		 "Problem in lubksb",error_message);
    if (maxstep!=NULL){
      //impose maximum stepping for Newton:
      for(i=0; i<neq; i++){
	mdelta_y[i+1] = min(mdelta_y[i+1],maxstep[i]);
	mdelta_y[i+1] = max(mdelta_y[i+1],-maxstep[i]);
      }
    }
    reldif = 0.0;
    for (i=0; i<neq; i++){
      y0[i] -= mdelta_y[i+1];
      reldif = max(reldif,fabs(mdelta_y[i+1]/(y0[i]+DBL_MIN)));
    }
    if (reldif<rtol)
      break;
  }
  free(Fval);
  free(vv);
  for (i=1; i<=neq; i++) free(jac[i]);
  free(jac);
  free(indx);
  free(mdelta_y);
  if (converged==_TRUE_)
    return _SUCCESS_;
  else{
    sprintf(error_message,"Newton failed to converge, index of failure: %d, dy=%g, y = %g\n",err_idx,-mdelta_y[err_idx+1],y0[err_idx]); 
    return _FAILURE_;
  }
}

int jacobian_for_Newton(int (*vecfun)(double * y, double * Fy, void *param),
			double *y0,
			double *Fval,
			size_t neq,
			void *param,
			double **jac){
  double del;
  double *Fy;
  int i,col;
  Fy = malloc(sizeof(double)*neq);

  for(col=0; col < neq; col++){
    del = 1e-6*fabs(y0[col]);
    if (del==0.0) del=1e-6;
    //del = 1e-6;
    y0[col] += del;
    //if (col==1) printf("y0[1]=%g, del[1] = %g\n",y0[1],del);
    vecfun(y0,Fy,param);
    y0[col] -= del;
    for (i=0; i<neq; i++){
      jac[i+1][col+1] = (Fy[i]-Fval[i])/del;
    }
  }
  free(Fy);
  return _SUCCESS_;
}
