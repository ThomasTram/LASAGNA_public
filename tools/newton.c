#include "newton.h"

int Newton(int (*vecfun)(double * y, double * Fy, void *param),
	   double *y0,
	   void *param,
	   double *maxstep,
	   double rtol,
	   int *iter,
	   int maxiter,
	   int neq,
	   ErrorMsg error_message){
  double *Fval,**jac,lusign,*vv,*mdelta_y;
  int i,*indx;
  int converged=_TRUE_;
  Fval = malloc(sizeof(double)*neq);
  vv = malloc(sizeof(double)*(neq+1));
  jac = malloc(sizeof(double*)*(neq+1));
  indx = malloc(sizeof(int)*(neq+1));
  mdelta_y = malloc(sizeof(double)*(neq+1));
  for (i=1; i<=neq; i++) jac[i] = malloc(sizeof(double)*(neq+1));

  for (*iter=1; *iter<=maxiter; (*iter)++){
    //Do Newton iteration:
    converged = _TRUE_;
    vecfun(y0,Fval,param);
    jacobian_for_Newton(vecfun,y0,Fval,neq,param,jac);
    lasagna_call(ludcmp(jac, neq, indx, &lusign, vv),
		 "Problem in ludcmp",error_message);
    for(i=0; i<neq; i++) mdelta_y[i+1] = Fval[i];
    double tmp;
    int j;
    for (j=1; j<neq; j++){
      for (i=1; i<neq; i++){
	tmp = jac[j][i];
      }
    }
    lasagna_call(lubksb(jac, neq, indx, mdelta_y),
		 "Problem in lubksb",error_message);
    if (maxstep!=NULL){
      //impose maximum stepping for Newton:
      for(i=0; i<neq; i++){
	mdelta_y[i+1] = min(mdelta_y[i+1],maxstep[i]);
	mdelta_y[i+1] = max(mdelta_y[i+1],-maxstep[i]);
      }
    }
    for (i=0; i<neq; i++){
      y0[i] -= mdelta_y[i+1];
      if (fabs(mdelta_y[i+1]/y0[i])>rtol)
	converged = _FALSE_;
    }
    if (converged == _TRUE_)
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
  else
    return _FAILURE_;
}

int jacobian_for_Newton(int (*vecfun)(double * y, double * Fy, void *param),
			double *y0,
			double *Fval,
			int neq,
			void *param,
			double **jac){
  double del;
  double *Fy;
  int i,col;
  Fy = malloc(sizeof(double)*neq);

  for(col=0; col < neq; col++){
    del = 1e-6*fabs(y0[col]);
    if (del==0.0) del=1e-15;
    y0[col] += del;
    vecfun(y0,Fy,param);
    y0[col] -= del;
    for (i=0; i<neq; i++){
      jac[i+1][col+1] = (Fy[i]-Fval[i])/del;
    }
  }
  free(Fy);
  return _SUCCESS_;
}
