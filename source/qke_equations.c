/** @file qke_equations.c 
 * Thomas Tram Bülow, 25.08.2011    
 * Modified Rasmus Sloth Hansen, 2013
 */
#include "common.h"
#include "qke_equations.h"
int init_qke_param(qke_param *pqke){
  int i,j,k,idx,nz;
  size_t neq;
  double k1,k2;
  double Nres, vres, Tres;
  int **J;
  Nres = pqke->Nres;
  vres = pqke->vres;
  Tres = pqke->Tres;
  pqke->Tvec = malloc(sizeof(double)*Tres);
  pqke->xi = malloc(sizeof(double)*Nres);
  pqke->ui = malloc(sizeof(double)*Nres);
  pqke->vi = malloc(sizeof(double)*Nres);
  pqke->duidT = malloc(sizeof(double)*Nres);
  pqke->dvidT = malloc(sizeof(double)*Nres);
  pqke->duidx = malloc(sizeof(double)*Nres);
  pqke->dxidT = malloc(sizeof(double)*Nres);
  pqke->a = malloc(sizeof(double)*Nres);
  pqke->y_0 = malloc(sizeof(double)*(1+Nres));
  pqke->maxstep = malloc(sizeof(double)*(1+Nres));
  pqke->x_grid = malloc(sizeof(double)*vres);
  pqke->u_grid = malloc(sizeof(double)*vres);
  pqke->v_grid = malloc(sizeof(double)*vres);
  pqke->dvdu_grid = malloc(sizeof(double)*vres);
  pqke->dudT_grid = malloc(sizeof(double)*vres);
  pqke->mat = malloc(sizeof(double*)*(Nres+2));
  for (i=0; i<(Nres+2); i++) 
    pqke->mat[i] = malloc(sizeof(double)*(Nres + 2));
  pqke->vv = malloc(sizeof(double)*(Nres + 2));
  pqke->indx = malloc(sizeof(int)*(Nres + 2));

  //Do some calculations for the u(x) mapping:
  k1 = pqke->xmin/pqke->xext;
  k2 = pqke->xmax/pqke->xext;
  pqke->eps2 = (1+k1)/(k2-k1);
  pqke->eps1 = k1*(1+pqke->eps2);

  //Some secondary initialisations:
  for(i=0; i<vres; i++){
    pqke->v_grid[i] = pqke->v_left+
      i*(pqke->v_right-pqke->v_left)/(vres-1.0); 
  }
  for(i=0; i<Tres; i++){
    pqke->Tvec[i] = pqke->T_initial+
      i*(pqke->T_final-pqke->T_initial)/(Tres-1.0); 
  }
  if (pqke->is_electron == _TRUE_){
     pqke->C_alpha = 1.27;
  }
  else{
     pqke->C_alpha = 0.92;
  }
  pqke->guess_exists = _FALSE_;
  
  //Set up the indices:
  idx = 0;
  pqke->index_L = idx;
  idx++;
  pqke->index_Pa_plus = idx;
  idx+=pqke->vres;
  pqke->index_Pa_minus =idx;
  idx +=pqke->vres;
  pqke->index_Ps_plus = idx;
  idx +=pqke->vres;
  pqke->index_Ps_minus = idx;
  idx +=pqke->vres;
  pqke->index_Px_plus = idx;
  idx +=pqke->vres;
  pqke->index_Px_minus = idx;
  idx +=pqke->vres;
  pqke->index_Py_plus = idx;
  idx +=pqke->vres;
  pqke->index_Py_minus = idx;
  idx +=pqke->vres;

  //Last thing to do:
  neq = idx;
  pqke->neq = neq;

  //Pattern for Jacobi matrix:
  pqke->Ap = malloc(sizeof(int)*(neq+1));
  pqke->Ai = malloc(sizeof(int)*neq*neq);

  //Construct Jacobian pattern:
  J = malloc(sizeof(int *)*neq);
  J[0] = calloc(neq*neq,sizeof(int)); 
 for (i=1; i<neq; i++)
   J[i] = J[i-1]+neq;
  
  //Establish diagonal (just to be sure)
  for (i=0; i<neq; i++)
    J[i][i] = 1;
  //Everything couples to L through the parametrisation:
  for (i=0; i<neq; i++)
    J[i][pqke->index_L] = 1;


  //F-L dependence:
  for (i=0; i<vres; i++)
    J[pqke->index_L][pqke->index_Py_minus+i] = 1;
  //Loop over grid:
  for (j=0; j<vres; j++){
    //F-Pa_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Pa_plus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Pa_plus+j][pqke->index_Py_plus+j] = 1;
    //F-Pa_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Pa_minus+j][pqke->index_Pa_minus+k] = 1;
    J[pqke->index_Pa_minus+j][pqke->index_Py_minus+j] = 1;
    //F-Ps_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Ps_plus+j][pqke->index_Ps_plus+k] = 1;
    J[pqke->index_Ps_plus+j][pqke->index_Py_plus+j] = 1;
    //R_nus_plus dependence:
    for (k=0; k<vres; k++){
      J[pqke->index_Ps_plus+j][pqke->index_Ps_minus+k] = 1;
      J[pqke->index_Ps_plus+j][pqke->index_Ps_plus+k] = 1;
    }
    //F-Ps_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Ps_minus+j][pqke->index_Ps_minus+k] = 1;
    J[pqke->index_Ps_minus+j][pqke->index_Py_minus+j] = 1;
    //R_nus_minus dependence:
    J[pqke->index_Ps_minus+j][pqke->index_Pa_minus+j] = 1;
    //F-Px_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Px_plus+j][pqke->index_Px_plus+k] = 1;
    J[pqke->index_Px_plus+j][pqke->index_L] = 1;
    for (k=0; k<vres; k++)
      J[pqke->index_Px_plus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Px_plus+j][pqke->index_Py_plus+j] = 1;
    J[pqke->index_Px_plus+j][pqke->index_Py_minus+j] = 1;
    //F-Px_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Px_minus+j][pqke->index_Px_minus+k] = 1;
    J[pqke->index_Px_minus+j][pqke->index_L] = 1;
    for (k=0; k<vres; k++)
      J[pqke->index_Px_minus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Px_minus+j][pqke->index_Py_plus+j] = 1;
    J[pqke->index_Px_minus+j][pqke->index_Py_minus+j] = 1;
    //F-Py_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Py_plus+j][pqke->index_Py_plus+k] = 1;
    J[pqke->index_Py_plus+j][pqke->index_L] = 1;
    for (k=0; k<vres; k++)
      J[pqke->index_Py_plus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Ps_plus+j] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Px_plus+j] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Px_minus+j] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Py_plus+j] = 1;
    //F-Py_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Py_minus+j][pqke->index_Py_minus+k] = 1;
    J[pqke->index_Py_minus+j][pqke->index_L] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Pa_minus+j] = 1;
    for (k=0; k<vres; k++)
      J[pqke->index_Py_minus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Ps_minus+j] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Px_plus+j] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Px_minus+j] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Py_minus+j] = 1;
  }


  //Store pattern in sparse column compressed form:
  pqke->Ap[0] = 0;
  nz = 0;
  for (i=0; i<neq; i++){
    for (j=0; j<neq; j++){
      if (J[j][i] == 1){
	pqke->Ai[nz] = j;
	nz++;
      }
    }
    pqke->Ap[i+1] = nz;
  }
  pqke->Ai = realloc(pqke->Ai,sizeof(int)*nz);

  free(J[0]);
  free(J);
  return _SUCCESS_;
};

int free_qke_param(qke_param *pqke){
  int i;
  free(pqke->xi);
  free(pqke->ui);
  free(pqke->vi);
  free(pqke->Tvec);
  free(pqke->duidT);
  free(pqke->dvidT);
  free(pqke->duidx);
  free(pqke->dxidT);
  free(pqke->a);
  free(pqke->y_0);
  free(pqke->maxstep);
  free(pqke->x_grid);
  free(pqke->u_grid);
  free(pqke->v_grid);
  free(pqke->dvdu_grid);
  free(pqke->dudT_grid);
  for (i=0; i<(pqke->Nres+2); i++) 
    free(pqke->mat[i]);
  free(pqke->mat);
  free(pqke->vv);
  free(pqke->indx);
  free(pqke->Ap);
  free(pqke->Ai);
  background_free_dof(&(pqke->pbs));
  return _SUCCESS_;
};

int init_qke_param_fixed_grid(qke_param *pqke){
  int i,j,k,idx,nz;
  size_t neq;
  double k1,k2;
  double Nres, vres, Tres;
  int **J;
  Nres = pqke->Nres;
  vres = pqke->vres;
  Tres = pqke->Tres;
  pqke->Tvec = malloc(sizeof(double)*Tres);
  pqke->xi = malloc(sizeof(double)*Nres);
  pqke->ui = malloc(sizeof(double)*Nres);
  pqke->vi = malloc(sizeof(double)*Nres);
  pqke->duidT = malloc(sizeof(double)*Nres);
  pqke->dvidT = malloc(sizeof(double)*Nres);
  pqke->duidx = malloc(sizeof(double)*Nres);
  pqke->dxidT = malloc(sizeof(double)*Nres);
  pqke->a = malloc(sizeof(double)*Nres);
  pqke->y_0 = malloc(sizeof(double)*(1+Nres));
  pqke->maxstep = malloc(sizeof(double)*(1+Nres));
  pqke->x_grid = malloc(sizeof(double)*vres);
  pqke->u_grid = malloc(sizeof(double)*vres);
  pqke->v_grid = malloc(sizeof(double)*vres);
  pqke->dvdu_grid = malloc(sizeof(double)*vres);
  pqke->dudT_grid = malloc(sizeof(double)*vres);
  pqke->mat = malloc(sizeof(double*)*(Nres+2));
  for (i=0; i<(Nres+2); i++) 
    pqke->mat[i] = malloc(sizeof(double)*(Nres + 2));
  pqke->vv = malloc(sizeof(double)*(Nres + 2));
  pqke->indx = malloc(sizeof(int)*(Nres + 2));

  //Do some calculations for the u(x) mapping:
  k1 = pqke->xmin/pqke->xext;
  k2 = pqke->xmax/pqke->xext;
  pqke->eps2 = (1+k1)/(k2-k1);
  pqke->eps1 = k1*(1+pqke->eps2);

  //Some secondary initialisations:
  for(i=0; i<vres; i++){
    pqke->v_grid[i] = pqke->v_left+
      i*(pqke->v_right-pqke->v_left)/(vres-1.0); 
  }
  for(i=0; i<Tres; i++){
    pqke->Tvec[i] = pqke->T_initial+
      i*(pqke->T_final-pqke->T_initial)/(Tres-1.0); 
  }
  if (pqke->is_electron == _TRUE_){
     pqke->C_alpha = 1.27;
  }
  else{
     pqke->C_alpha = 0.92;
  }
  pqke->guess_exists = _FALSE_;
  
  //Set up the indices:
  idx = 0;
  pqke->index_L = idx;
  idx++;
  pqke->index_Pa_plus = idx;
  idx+=pqke->vres;
  pqke->index_Pa_minus =idx;
  idx +=pqke->vres;
  pqke->index_Ps_plus = idx;
  idx +=pqke->vres;
  pqke->index_Ps_minus = idx;
  idx +=pqke->vres;
  pqke->index_Px_plus = idx;
  idx +=pqke->vres;
  pqke->index_Px_minus = idx;
  idx +=pqke->vres;
  pqke->index_Py_plus = idx;
  idx +=pqke->vres;
  pqke->index_Py_minus = idx;
  idx +=pqke->vres;

  //Last thing to do:
  neq = idx;
  pqke->neq = neq;

  //Pattern for Jacobi matrix:
  pqke->Ap = malloc(sizeof(int)*(neq+1));
  pqke->Ai = malloc(sizeof(int)*neq*neq);

  //Construct Jacobian pattern:
  J = malloc(sizeof(int *)*neq);
  J[0] = calloc(neq*neq,sizeof(int)); 
 for (i=1; i<neq; i++)
   J[i] = J[i-1]+neq;
  
  //Establish diagonal (just to be sure)
  for (i=0; i<neq; i++)
    J[i][i] = 1;
  
  //F-L dependence:
  for (i=0; i<vres; i++)
    J[pqke->index_L][pqke->index_Py_minus+i] = 1;
  //Loop over grid:
  for (j=0; j<vres; j++){
    //F-Pa_plus[j] dependence:
    J[pqke->index_Pa_plus+j][pqke->index_L] = 1;
    J[pqke->index_Pa_plus+j][pqke->index_Pa_plus+j] = 1;
    J[pqke->index_Pa_plus+j][pqke->index_Py_plus+j] = 1;
    //F-Pa_minus[j] dependence:
    J[pqke->index_Pa_minus+j][pqke->index_L] = 1;
    J[pqke->index_Pa_minus+j][pqke->index_Pa_minus+j] = 1;
    J[pqke->index_Pa_minus+j][pqke->index_Py_minus+j] = 1;
    //F-Ps_plus[j] dependence:
    J[pqke->index_Ps_plus+j][pqke->index_Py_plus+j] = 1;
    //F-Ps_minus[j] dependence:
    J[pqke->index_Ps_minus+j][pqke->index_Py_minus+j] = 1;
    //F-Px_plus[j] dependence:
    J[pqke->index_Px_plus+j][pqke->index_Px_plus+j] = 1;
    J[pqke->index_Px_plus+j][pqke->index_L] = 1;
    //    for (k=0; k<vres; k++)
    //  J[pqke->index_Px_plus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Px_plus+j][pqke->index_Py_plus+j] = 1;
    J[pqke->index_Px_plus+j][pqke->index_Py_minus+j] = 1;
    //F-Px_minus[j] dependence:
    J[pqke->index_Px_minus+j][pqke->index_Px_minus+j] = 1;
    J[pqke->index_Px_minus+j][pqke->index_L] = 1;
    for (k=0; k<vres; k++)
      J[pqke->index_Px_minus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Px_minus+j][pqke->index_Py_plus+j] = 1;
    J[pqke->index_Px_minus+j][pqke->index_Py_minus+j] = 1;
    //F-Py_plus[j] dependence:
    J[pqke->index_Py_plus+j][pqke->index_Py_plus+j] = 1;
    J[pqke->index_Py_plus+j][pqke->index_L] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Pa_plus+j] = 1;
    for (k=0; k<vres; k++)
      J[pqke->index_Py_plus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Ps_plus+j] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Px_plus+j] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Px_minus+j] = 1;
    J[pqke->index_Py_plus+j][pqke->index_Py_plus+j] = 1;
    //F-Py_minus[j] dependence:
    J[pqke->index_Py_minus+j][pqke->index_Py_minus+j] = 1;
    J[pqke->index_Py_minus+j][pqke->index_L] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Pa_minus+j] = 1;
    for (k=0; k<vres; k++)
      J[pqke->index_Py_minus+j][pqke->index_Pa_plus+k] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Ps_minus+j] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Px_plus+j] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Px_minus+j] = 1;
    J[pqke->index_Py_minus+j][pqke->index_Py_minus+j] = 1;
  }


  //Store pattern in sparse column compressed form:
  pqke->Ap[0] = 0;
  nz = 0;
  for (i=0; i<neq; i++){
    for (j=0; j<neq; j++){
      if (J[j][i] == 1){
	pqke->Ai[nz] = j;
	nz++;
      }
    }
    pqke->Ap[i+1] = nz;
  }
  pqke->Ai = realloc(pqke->Ai,sizeof(int)*nz);

  free(J[0]);
  free(J);
  return _SUCCESS_;
};


int get_resonances_xi(double T, 
		      double L,
		      qke_param *pqke){
  double *xi=pqke->xi;
  double x0, x0b,A;
  int i;

  x0 = sqrt(fabs(pqke->V0/pqke->V1));
  
  xi[0] = x0;
  xi[1] = x0;
  A = fabs(0.5*pqke->VL/sqrt(fabs(pqke->V0*pqke->V1)));
  if (pqke->delta_m2<0.0){
    //Case of inverted Hierarchy
    xi[0] *= (-A+sqrt(1.0+A*A));
    xi[1] *=  (A+sqrt(1.0+A*A));
  }
  else{
    //Case of normal hierarchy:
    if (A>1.0){
      xi[0] *=(A-sqrt(A*A-1.0));
      xi[1] *=(A+sqrt(A*A-1.0));
    }
  }
  // Protect against possible seg fault:
  for (i=0; i<pqke->Nres; i++){
    if(xi[i]>pqke->xmax){
      xi[i] = pqke->xmax;
    }
  }
  return _SUCCESS_;
};

int get_resonances_dxidT(double T, 
			 double L,
			 double dLdT,
			 qke_param *pqke){
  double x0, A,A2;
  double *dxidT=pqke->dxidT;
  double F;
  double one_plus_dlogLdlogT;
  int i;

  x0 = sqrt(fabs(pqke->V0/pqke->V1));  
 
  dxidT[0] = -3.0*x0/T;
  dxidT[1] = -3.0*x0/T;

  if (fabs(L)>1e-100){
 
    A = fabs(0.5*pqke->VL/sqrt(fabs(pqke->V0*pqke->V1)));
    A2 = A*A;
    one_plus_dlogLdlogT = 1.0+T/L*dLdT;
 
    if (pqke->delta_m2<0.0){
      //Case of inverted Hierarchy
      F = -A+sqrt(1.0+A2);
      dxidT[0] *= (F-A/3.0*one_plus_dlogLdlogT*(-1.0+pow(1.0+1.0/A2,-0.5)));
      F = A+sqrt(1.0+A2);
      dxidT[1] *= (F-A/3.0*one_plus_dlogLdlogT*(1.0+pow(1.0+1.0/A2,-0.5)));
    }
    else{
      //Case of normal hierarchy:
      if (A>1.0){
	F = A-sqrt(A2-1.0);
	dxidT[0] *= (F-A/3.0*one_plus_dlogLdlogT*(1.0-pow(1.0-1.0/A2,-0.5)));
	F = A+sqrt(A2-1.0);
	dxidT[1] *= (F-A/3.0*one_plus_dlogLdlogT*(1.0+pow(1.0-1.0/A2,-0.5))); 
      }
    }
  }
  //Make it selfconsistent:
  for (i=0; i<pqke->Nres; i++){
    if(pqke->xi[i]>=(pqke->xmax-1e-12)){
      dxidT[i] = 0.0;
    }
  }
  return _SUCCESS_;
};


int u_of_x(double x, double *u, double *dudx, qke_param *pqke){
  double K,xmin,xmax,xext;
  xmin = pqke->xmin;
  xmax = pqke->xmax;
  xext = pqke->xext;
  K = (xext+xmax)/(xmax-xmin);
  *u = K*(x-xmin)/(x+xext);
  *dudx = K*(xext+xmin)/pow(x+xext,2);

  return _SUCCESS_;
};

int x_of_u(double u, double *x, qke_param *param){
  *x = param->xext*(u+param->eps1)/(1+param->eps2-u);
  return _SUCCESS_;
};

int nonlinear_rhs(double *y, double *Fy, void *param){
  qke_param *pqke=param;
  double b=y[0];
  double *vi=y+1;
  double alpha = pqke->alpha;
  double *ui=pqke->ui;
  double *a=pqke->a;
  int i,n=pqke->Nres;

  for (i=0; i<n; i++){
    a[i] = ui[i]-alpha*vi[i];
  }
  Fy[0] = a[0]-b*pow(vi[0],3);
  Fy[1] = alpha+a[n-1]+b*pow(1-vi[n-1],3)-1.0;
  for (i=0; i<n-1; i++){
    Fy[2+i] = a[i]-a[i+1]+0.25*b*pow(vi[i+1]-vi[i],3);
  }
  return _SUCCESS_;
}

int nonlinear_rhs_jac(double *y, double **jac, void *param){
  qke_param *pqke=param;
  double b=y[0];
  double *vi=y+1;
  double alpha = pqke->alpha;
  int i,j,n=pqke->Nres;
  
  for (i=1; i<=n+1; i++){
    for (j=1; j<=n+1; j++){
      jac[i][j] = 0.0;
    }
  }
  jac[1][1]=-pow(vi[0],3); jac[1][2] = -alpha-3.0*b*vi[0]*vi[0];
  jac[2][1]=pow(1.0-vi[n-1],3);jac[2][n+1]=-alpha-3.0*b*pow(1.0-vi[n-1],2);
  for (i=0; i<n-1; i++){
    jac[i+3][1] = 0.25*pow(vi[i+1]-vi[i],3);
    jac[i+3][i+2] = -alpha-0.75*b*pow(vi[i+1]-vi[i],2);
    jac[i+3][i+3] = -jac[i+3][i+2];
  }
  
  return _SUCCESS_;
}


int qke_initial_conditions(double Ti, double *y, qke_param *pqke){
  /** Set initial conditions at temperature Ti: */
  int i;
  ErrorMsg error_message;
  double x, L, n_plus=2.0;
  double Vx,D,Vz,Vz_bar,Px,Py,Px_bar,Py_bar;

  //Assuming y is calloc'ed -- dangerous, better to zero it.
  for (i=0; i<pqke->neq; i++) y[i] = 0.0;

  L = pqke->L_initial;
  //Set standard equilibrium initial conditions:
  y[pqke->index_L] = L/_L_SCALE_;
  for (i=0; i<pqke->vres; i++)
    y[pqke->index_Pa_plus+i] = 4.0;
  
/** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  if (pqke->is_electron==_TRUE_)
    pqke->g_alpha = 1.0+4.0/((1.0-_SIN2_THETA_W_)*n_plus);
  else
    pqke->g_alpha = 1.0;
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(Ti,3)/_PI_/_PI_*L;
  pqke->Vx = pqke->delta_m2/(2.0*Ti)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*Ti)*cos(2.0*pqke->theta_zero);
  pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
    _G_F_/_M_Z_/_M_Z_*pow(Ti,5)*n_plus*pqke->g_alpha;

  get_resonances_xi(Ti,L,pqke);  
  //Get new x_grid
  lasagna_call(get_parametrisation(Ti,
				   pqke, 
				   error_message),
	       error_message,error_message);
  

  for (i=0; i<pqke->vres; i++){
    x = pqke->x_grid[i];
    Vx = pqke->Vx/x;
    Vz = pqke->V0/x + pqke->V1*x + pqke->VL;
    Vz_bar = pqke->V0/x + pqke->V1*x - pqke->VL;
    D = 0.5*pqke->C_alpha*_G_F_*_G_F_*x*pow(Ti,5);

    Px = Vx*Vz/(D*D+Vz*Vz);
    Px_bar = Vx*Vz_bar/(D*D+Vz_bar*Vz_bar);
    Py = -Vx*D/(D*D+Vz*Vz);
    Py_bar = -Vx*D/(D*D+Vz_bar*Vz_bar);

    //continue;
    y[pqke->index_Px_plus+i] =  Px + Px_bar;
    y[pqke->index_Px_minus+i] = Px - Px_bar;
    y[pqke->index_Py_plus+i] = Py + Py_bar;
    y[pqke->index_Py_minus+i] = Py - Py_bar;
 
  }
  //initial conditions for qke_stop_at_divL
  pqke->max_old = L;
  pqke->max_cur = L;
  pqke->should_break = _FALSE_;
  pqke->breakpoint = 0;
  return _SUCCESS_;
};

int qke_init_output(qke_param *pqke){
  int Tres=pqke->Tres;
  int vres=pqke->vres;
  int Nres=pqke->Nres;
  int handle;
  double tmp_array[3];
  int32_t tmp_array_int[3];
  char *outf=pqke->output_filename;
  //Initialises the output file and stores parameters
  mat_create_file(outf);
  //Add matrices that are defined on momentum grid:
  mat_add_matrix(outf,"Pa_plus",miDOUBLE,Tres,vres,&(pqke->Pa_plus_handle));
  mat_add_matrix(outf,"Pa_minus",miDOUBLE,Tres,vres,&(pqke->Pa_minus_handle));
  mat_add_matrix(outf,"Ps_plus",miDOUBLE,Tres,vres,&(pqke->Ps_plus_handle));
  mat_add_matrix(outf,"Ps_minus",miDOUBLE,Tres,vres,&(pqke->Ps_minus_handle));
  mat_add_matrix(outf,"Px_plus",miDOUBLE,Tres,vres,&(pqke->Px_plus_handle));
  mat_add_matrix(outf,"Px_minus",miDOUBLE,Tres,vres,&(pqke->Px_minus_handle));
  mat_add_matrix(outf,"Py_plus",miDOUBLE,Tres,vres,&(pqke->Py_plus_handle));
  mat_add_matrix(outf,"Py_minus",miDOUBLE,Tres,vres,&(pqke->Py_minus_handle));
  mat_add_matrix(outf,"x_grid",miDOUBLE,Tres,vres,&(pqke->x_grid_handle));
  mat_add_matrix(outf,"u_grid",miDOUBLE,Tres,vres,&(pqke->u_grid_handle));
  mat_add_matrix(outf,"v_grid",miDOUBLE,Tres,vres,&(pqke->v_grid_handle));
  //Add resonance dependent matrices:
  mat_add_matrix(outf,"xi_vec",miDOUBLE,Tres,Nres,&(pqke->xi_handle));
  mat_add_matrix(outf,"ui_vec",miDOUBLE,Tres,Nres,&(pqke->ui_handle));
  mat_add_matrix(outf,"vi_vec",miDOUBLE,Tres,Nres,&(pqke->vi_handle));
  mat_add_matrix(outf,"b_a_vec",miDOUBLE,Tres,(1+Nres),&(pqke->b_a_vec_handle));
  //Add other matrices:
  mat_add_matrix(outf,"L_vec",miDOUBLE,Tres,1,&(pqke->L_handle));
  mat_add_matrix(outf,"T_vec",miDOUBLE,Tres,1,&(pqke->T_handle));
  mat_add_matrix(outf,"I_conserved",miDOUBLE,Tres,1,&(pqke->I_conserved_handle));
  mat_add_matrix(outf,"V0_vec",miDOUBLE,Tres,1,&(pqke->V0_handle));
  mat_add_matrix(outf,"V1_vec",miDOUBLE,Tres,1,&(pqke->V1_handle));
  mat_add_matrix(outf,"Vx_vec",miDOUBLE,Tres,1,&(pqke->Vx_handle));
  mat_add_matrix(outf,"VL_vec",miDOUBLE,Tres,1,&(pqke->VL_handle));
  //Add constant parameters:
  mat_add_matrix(outf,"L_initial",miDOUBLE,1,1,&handle);
  mat_add_matrix(outf,"delta_m2_theta_zero",miDOUBLE,1,2,&handle);
  mat_add_matrix(outf,"is_electron",miINT32,1,1,&handle);
  mat_add_matrix(outf,"Tres_vres",miINT32,1,2,&handle);
  mat_add_matrix(outf,"xmin_xext_xmax",miDOUBLE,1,3,&handle);
  mat_add_matrix(outf,"alpha_rs",miDOUBLE,1,2,&handle);
  //Write parameters and stuff we know in advance:
  //Temperature - we could write it here, but its nice to have non-computed
  //T values = 0:
  //mat_write_data(outf,"T",pqke->Tvec,0,Tres);
  tmp_array[0] = pqke->delta_m2; tmp_array[1] = pqke->theta_zero;
  mat_write_data(outf,"delta_m2_theta_zero",&(tmp_array),0,2);
  mat_write_data(outf,"is_electron",&(pqke->is_electron),0,1);
  tmp_array_int[0]=pqke->Tres; tmp_array_int[1]=pqke->vres;
  mat_write_data(outf,"Tres_vres",tmp_array_int,0,2);
  tmp_array[0] = pqke->xmin; tmp_array[1] = pqke->xext; 
  tmp_array[2] = pqke->xmax;
  mat_write_data(outf,"xmin_xext_xmax",tmp_array,0,3);
  tmp_array[0] = pqke->alpha; tmp_array[1] = pqke->rs;
  mat_write_data(outf,"alpha_rs",tmp_array,0,2);
  return _SUCCESS_;
}


int qke_store_output(double T,
			    double *y,
			    double *dy,
			    int index_t,
			    void *param,
			    ErrorMsg error_message){
  qke_param *pqke=param;
  int vres=pqke->vres;
  int Nres=pqke->Nres;
  int i;
  double x,xp1,f0,f0p1,Pa_minus,Ps_minus,Ps_minusp1,Pa_minusp1,I_PaPs,L; 
  FILE *mat_file;
  /** Calculate integrated quantities for convenience: */
 
  I_PaPs = 0.0;
  for (i=0; i<pqke->vres-1; i++){
    x = pqke->x_grid[i];
    xp1 = pqke->x_grid[i+1];
    f0 = 1.0/(1.0+exp(x));
    f0p1 = 1.0/(1.0+exp(xp1));
    Ps_minus = y[pqke->index_Py_minus+i];
    Ps_minusp1 = y[pqke->index_Py_minus+i+1];
    Pa_minus = y[pqke->index_Pa_plus+i];
    Pa_minusp1 = y[pqke->index_Pa_plus+i+1];
    I_PaPs += 0.5*(xp1-x)*(x*x*f0*(Ps_minus+Pa_minus)+
			   xp1*xp1*f0p1*(Ps_minusp1+Pa_minusp1));
  }

  //printf("Storing output at index: %d\n",index_t);
  mat_file = fopen(pqke->output_filename,"r+b");
  //Write stuff from y_vector:
  L = y[pqke->index_L]*_L_SCALE_;
  mat_write_fast(&L,pqke->L_handle,8,1);
  mat_write_fast(y+pqke->index_Pa_plus,pqke->Pa_plus_handle,8,vres);
  mat_write_fast(y+pqke->index_Pa_minus,pqke->Pa_minus_handle,8,vres);
  mat_write_fast(y+pqke->index_Ps_plus,pqke->Ps_plus_handle,8,vres);
  mat_write_fast(y+pqke->index_Ps_minus,pqke->Ps_minus_handle,8,vres);
  mat_write_fast(y+pqke->index_Px_plus,pqke->Px_plus_handle,8,vres);
  mat_write_fast(y+pqke->index_Px_minus,pqke->Px_minus_handle,8,vres);
  mat_write_fast(y+pqke->index_Py_plus,pqke->Py_plus_handle,8,vres);
  mat_write_fast(y+pqke->index_Py_minus,pqke->Py_minus_handle,8,vres);
  //Write stuff from structure:
  mat_write_fast(pqke->x_grid,pqke->x_grid_handle,8,vres);
  mat_write_fast(pqke->u_grid,pqke->u_grid_handle,8,vres);
  mat_write_fast(pqke->v_grid,pqke->v_grid_handle,8,vres);
  mat_write_fast(pqke->xi,pqke->xi_handle,8,Nres);
  mat_write_fast(pqke->ui,pqke->ui_handle,8,Nres);
  mat_write_fast(pqke->vi,pqke->vi_handle,8,Nres);
  mat_write_fast(&(I_PaPs),pqke->I_conserved_handle,8,1);
  mat_write_fast(&(pqke->V0),pqke->V0_handle,8,1);
  mat_write_fast(&(pqke->V1),pqke->V1_handle,8,1);
  mat_write_fast(&(pqke->Vx),pqke->Vx_handle,8,1);
  mat_write_fast(&(pqke->VL),pqke->VL_handle,8,1);
  mat_write_fast(&(pqke->b),pqke->b_a_vec_handle,8,1);
  mat_write_fast(pqke->vi,pqke->b_a_vec_handle,8,Nres);//Wrong
  

  //Write temperature at last so we know that everything has been written:
  mat_write_fast(&T,pqke->T_handle,8,1);
  fclose(mat_file);
  return _SUCCESS_;
};


int qke_stop_at_L(double t,
		  double *y,
		  double *dy,
		  void *param,
		  ErrorMsg error_message){
  qke_param *pqke=param;
  if (fabs(y[pqke->index_L]*_L_SCALE_) >= pqke->L_final){
    printf("Value stop..\n");
    return _TRUE_;
  }
  /**  else if((y[pqke->index_L] != 0.0)&&
	  (dy[pqke->index_L]/y[pqke->index_L]>pqke->trigger_dLdT_over_L)){
    printf("Trigger stop..\n");
    return _TRUE_;
  }
  */
  else{
    return _FALSE_;
  }
};


int qke_stop_at_divL(double t,
		     double *y,
		     double *dy,
		     void *param,
		     ErrorMsg error_message){
  qke_param *pqke=param;
  double L, T, stopfactor;
  stopfactor = 2;
  L = y[pqke->index_L]*_L_SCALE_;
  T = t;
  // Check for sign change:
  if(L*pqke->max_cur < 0){
    pqke->max_old = max(pqke->max_old,fabs(pqke->max_cur));
    pqke->max_cur = 0;
    pqke->should_break = _FALSE_;
  }
  // Is the value larger than the current max:
  if((fabs(L)>fabs(pqke->max_cur)) && (pqke->should_break==_FALSE_)){
    pqke->max_cur = L;
    if((fabs(pqke->max_cur) > pqke->max_old*stopfactor)){
      pqke->breakpoint = T - pqke->T_wait;
      pqke->should_break = _TRUE_;
    }
  }
  else if(pqke->should_break==_TRUE_){
    if(fabs(L) < pqke->max_old)
      pqke->should_break = _FALSE_;
    else if(T<pqke->breakpoint)
      return _TRUE_;
  }
  return _FALSE_;
};


int qke_print_variables(double T,
			double *y,
			double *dy,
			void *param,
			ErrorMsg error_message){
  
  qke_param *pqke=param;
  int idx=53;
  double x,Vx,V0,V1,VL;
  double D,Gamma;
  double Pa_plus, Pa_minus, Ps_plus, Ps_minus, Px_plus, Px_minus;
  double Py_plus, Py_minus;
  x = pqke->x_grid[idx];
  Vx = pqke->Vx/x;
  V0 = pqke->V0/x;
  V1 = pqke->V1*x;
  VL = pqke->VL;
  Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
  D = 0.5*Gamma;
  
  Pa_plus = y[pqke->index_Pa_plus+idx];
  Pa_minus = y[pqke->index_Pa_minus+idx];
  Ps_plus = y[pqke->index_Ps_plus+idx];
  Ps_minus = y[pqke->index_Ps_minus+idx];
  Px_plus = y[pqke->index_Px_plus+idx];
  Px_minus = y[pqke->index_Px_minus+idx];
  Py_plus = y[pqke->index_Py_plus+idx];
  Py_minus = y[pqke->index_Py_minus+idx];


  fprintf(stderr,
	  "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	  1e3*T,
	  y[pqke->index_L]*_L_SCALE_,
	  Pa_plus,
	  Pa_minus,
	  Ps_plus,
	  Ps_minus,
	  Px_plus,
	  Px_minus,
	  Py_plus,
	  Py_minus,
	  dy[pqke->index_L]/y[pqke->index_L],
	  VL,
	  -(V0+V1)*Px_plus-VL*Px_minus, //15
	  0.5*Vx*(Pa_plus-Ps_plus),
	  D*Py_plus,
	  -(V0+V1)*Px_minus-VL*Px_plus,
	  0.5*Vx*(Pa_minus-Ps_minus),
	  D*Py_minus);

//qke_param *pqke=param;
  
  return _SUCCESS_;
};

int qke_print_L(double T,
			double *y,
			double *dy,
			void *param,
			ErrorMsg error_message){

  qke_param *pqke=param;
  fprintf(stderr,"%.16e %.16e\n",1e3*T,y[pqke->index_L]*_L_SCALE_);
  return _SUCCESS_;
}
int get_integrated_quantities(double *y,
			      qke_param *pqke,
			      double *I_VxPy_minus,
			      double *I_f0Pa_plus,
			      double *I_rho_ss,
			      double *I_rho_ss_bar,
			      double *I_f0,
			      ErrorMsg error_message){
  int i;
  double w_trapz, x, x2, f0, Vx;
  double Py_minus, Pa_plus, Ps, Ps_bar;
			      
  /** Integrated quantities needed. We integrate in x space: */
  *I_VxPy_minus = 0.0;
  *I_f0Pa_plus = 0.0;
  *I_rho_ss = 0.0;
  *I_rho_ss_bar = 0.0;
  *I_f0 = 0.0;
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(pqke->x_grid[i+1]-pqke->x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(pqke->x_grid[i]-pqke->x_grid[i-1]);
    else
      w_trapz = 0.5*(pqke->x_grid[i+1]-pqke->x_grid[i-1]);
    x = pqke->x_grid[i];
    x2 = x*x;
    f0 = 1.0/(1.0+exp(x));
    Vx = pqke->Vx/x;
    Py_minus = y[pqke->index_Py_minus+i];
    Pa_plus = y[pqke->index_Pa_plus+i];
    Ps = (y[pqke->index_Ps_plus+i] + y[pqke->index_Ps_minus+i]);
    Ps_bar = (y[pqke->index_Ps_plus+i] - y[pqke->index_Ps_minus+i]);
      
    *I_VxPy_minus += w_trapz*(x2*f0*Vx*Py_minus);
    *I_f0Pa_plus += w_trapz*(x2*f0*Pa_plus);
    *I_rho_ss += w_trapz*(x2*f0*Ps); //From equation 2.18 in KS01.
    *I_rho_ss_bar += w_trapz*(x2*f0*Ps_bar);
    *I_f0 += w_trapz*(x2*f0);
  }
  return _SUCCESS_;
}

int get_parametrisation(double T,qke_param *pqke, ErrorMsg error_message){
  double alpha=pqke->alpha;
  double *maxstep=pqke->maxstep;;
  double *y_0=pqke->y_0;
  double *ui=pqke->ui;
  double *vi=pqke->vi;
  double *x_grid=pqke->x_grid;
  double *u_grid=pqke->u_grid;
  double *v_grid=pqke->v_grid;
  double tol_newton=1e-12;
  double wi;
  int i,j;
  int niter;

  for (i=0; i<pqke->Nres; i++){
    u_of_x(pqke->xi[i],ui+i,pqke->duidx+i,pqke);
  }
  
  //Establish guess and set maximum steps for Newton method:
  /** This has been a very hard 'bug' to trace: Some of the evolvers rely
      on the derivative function to give exactly the same output for the same
      inputs, so making a non-deterministic guess here spoils this!
  */
  if (1==1){//(pqke->guess_exists == _FALSE_){
    y_0[0] = 1.0 - alpha;
    maxstep[0] = 100.0;
    for (i=1; i<=pqke->Nres; i++){
      y_0[i] = ui[i-1];
      maxstep[i] = 0.1;
    }
  }
  else{
    //We can make a more realistic guess:
    y_0[0] = pqke->b;
    maxstep[0] = 100.0;
    for (i=1; i<=pqke->Nres; i++){
      y_0[i] = vi[i-1];
      maxstep[i] = 0.1;
    }
  }
  //Find parametrisation parameters vi, a and b using Newton:
  lasagna_call(Newton(nonlinear_rhs,
		      nonlinear_rhs_jac,
		      y_0,
		      pqke,
		      maxstep,
		      tol_newton,
		      &niter,
		      100,
		      pqke->Nres+1,
		      error_message),
	       error_message,error_message);
  pqke->guess_exists = _TRUE_;
  pqke->b = y_0[0];
  for (i=0; i<pqke->Nres; i++){
    vi[i] = y_0[i+1];
    pqke->a[i] = ui[i]-alpha*vi[i];
  }

  //Now update grids:
  /** Find the splitting of v and store it temporarily in pqke->indx.
      We use the fact that v is uniform.
  */
  pqke->indx[0] = 0;
  pqke->indx[pqke->Nres] = pqke->vres;
  for (i=1; i<pqke->Nres; i++){
    wi = 0.5*(vi[i-1]+vi[i]); //Weighted average
    pqke->indx[i] = (int)((wi-pqke->v_left)/(v_grid[1]-v_grid[0]));
    if (pqke->indx[i]<pqke->indx[i-1])
      pqke->indx[i] = pqke->indx[i-1];
  }
  // Now loop over resonances:
  for (i=0; i<pqke->Nres; i++){
    if (pqke->indx[i]==pqke->indx[i+1])
      continue;
    //Loop over each segment:
    for (j=pqke->indx[i]; j<pqke->indx[i+1]; j++){
      u_grid[j] = alpha*v_grid[j]+pqke->a[i]+pqke->b*pow(v_grid[j]-vi[i],3);
      x_of_u(u_grid[j],&(x_grid[j]),pqke);
    }
  }

  return _SUCCESS_;
}

int get_partial_derivatives(double T,
			    double L,
			    double dLdT,
			    qke_param *pqke, 
			    ErrorMsg error_message){
  int i,j;
  double u1, v1, vN;
  double gamma_j, beta_j, wi;
  double alpha=pqke->alpha;
  double *dvidT=pqke->dvidT;
  double *duidT=pqke->duidT;
  double *v_grid=pqke->v_grid;
  double *ui=pqke->ui;
  double *vi=pqke->vi;
  double *dvdu_grid=pqke->dvdu_grid;
  double *dudT_grid=pqke->dudT_grid;
  double dbdT,daidT;
  double lu_sgn;

/** Get partial derivatives of vi with respect to T by solving a 
      linear system A*dvidT = B(duidT):
  */
  get_resonances_dxidT(T,L,dLdT,pqke);

  //Set partial derivatives of ui with respect to T:
  for (i=0; i<pqke->Nres; i++){ 
    pqke->duidT[i]=pqke->duidx[i]*pqke->dxidT[i];
  }
  for(j=0; j<pqke->Nres; j++){
    for(i=0; i<pqke->Nres; i++){
        pqke->mat[j+1][i+1] = 0.0;
    }
  }
  u1 = ui[0];
  v1 = vi[0];
  vN = vi[pqke->Nres-1];
  
  for (j=0; j<pqke->Nres-1; j++){
    /** Enter non-zero entries in A: */
    gamma_j = 0.25*pow(vi[j+1]/v1-vi[j]/v1,3)*(2.0*alpha-3.0*u1/v1);
    beta_j = alpha + 0.75*pqke->b*pow(vi[j+1]-vi[j],2);
    pqke->mat[j+1][1] = gamma_j;
    pqke->mat[j+1][j+1] -= beta_j;
    pqke->mat[j+1][j+2] += beta_j;
    /** Setup RHS: */
    dvidT[j] = duidT[j+1]-duidT[j]-0.25*pow(vi[j+1]/v1-vi[j]/v1,3)*duidT[0];
  }
  pqke->mat[pqke->Nres][1] = -pow((1.0-vN)/v1,3)*(2.0*alpha-3.0*u1/v1);
  pqke->mat[pqke->Nres][pqke->Nres] = alpha+3.0*pqke->b*pow(1.0-vN,2);
  dvidT[pqke->Nres-1] = duidT[pqke->Nres-1]+pow((1.0-vN)/v1,3)*duidT[0];

  //LU decomposition of matrix:
  lasagna_call(ludcmp(pqke->mat,pqke->Nres,pqke->indx,&lu_sgn,pqke->vv),
	       error_message,error_message);
  lasagna_call(lubksb(pqke->mat,pqke->Nres,pqke->indx,dvidT-1),
	       error_message,error_message);

  /** Calculate x and u on the v-grid along with derivatives dvdu and dudT:*/
  dbdT = (duidT[0]+(2.0*alpha-3.0*u1/v1)*dvidT[0])/pow(v1,3);
  /** Find the splitting of v and store it temporarily in pqke->indx.
      We use the fact that v is uniform.
  */
  pqke->indx[0] = 0;
  pqke->indx[pqke->Nres] = pqke->vres;
  for (i=1; i<pqke->Nres; i++){
    wi = 0.5*(vi[i-1]+vi[i]); //Weighted average
    pqke->indx[i] = (int)((wi-pqke->v_left)/(v_grid[1]-v_grid[0]));
    if (pqke->indx[i]<pqke->indx[i-1])
      pqke->indx[i] = pqke->indx[i-1];
  }

  // Now loop over resonances:
  for (i=0; i<pqke->Nres; i++){
    if (pqke->indx[i]==pqke->indx[i+1])
      continue;
    //Loop over each segment:
    for (j=pqke->indx[i]; j<pqke->indx[i+1]; j++){
      daidT = duidT[i]-alpha*dvidT[i];
      dvdu_grid[j] = 1.0/(alpha+3.0*pqke->b*pow(v_grid[j]-vi[i],2));
      dudT_grid[j] = daidT+dbdT*pow(v_grid[j]-vi[i],3)-
	3.0*pqke->b*pow(v_grid[j]-vi[i],2)*dvidT[i];
    }
  }
  /** Debug area*/ 
  /**
  fprintf(pqke->tmp,"%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	  T,
	  pqke->a[0],
	  pqke->a[1],
	  duidT[0]-alpha*dvidT[0],
	  duidT[1]-alpha*dvidT[1], 
	  pqke->b,
	  dbdT,
	  pqke->xi[0],
	  pqke->xi[1],
	  pqke->dxidT[0], 
	  pqke->dxidT[1],
	  pqke->ui[0],
	  pqke->ui[1],
	  pqke->duidT[0],
	  pqke->duidT[1],
	  pqke->vi[0],
	  pqke->vi[1],
	  pqke->dvidT[0],
	  pqke->dvidT[1]);
  */
return _SUCCESS_;
}


int qke_derivs(double T, 
	       double *y, 
	       double *dy, 
	       void *param,
	       ErrorMsg error_message){
  qke_param *pqke=param;
  double L;
  int i;
  double *v_grid=pqke->v_grid;
  double *dvdu_grid=pqke->dvdu_grid;
  double *dudT_grid=pqke->dudT_grid;
  double gentr,H;
  double Vx, VL;
  double x;
  double Gamma, D, V0, V1, Pa_plus, Pa_minus, Ps_plus, Ps_minus;
  double Px_plus, Px_minus, Py_plus, Py_minus, f0,  mu_div_T;
  double feq_plus, feq_minus, feq, feq_bar;
  double I_VxPy_minus, I_f0Pa_plus, I_rho_ss, I_rho_ss_bar, I_f0;
  double dudTdvdu, delta_v;
  double rs;
  int idx, stencil_method;
  double n_plus = 2.0;
  double dLdT;

  if (pqke->is_electron==_TRUE_)
    pqke->g_alpha = 1.0+4.0/((1.0-_SIN2_THETA_W_)*n_plus);
  else
    pqke->g_alpha = 1.0;
  
  L = y[pqke->index_L]*_L_SCALE_;
  
  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);
  pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
    _G_F_/_M_Z_/_M_Z_*pow(T,5)*n_plus*pqke->g_alpha;

  get_resonances_xi(T,L,pqke);
  
  //Get new x_grid
  lasagna_call(get_parametrisation(T,
				   pqke, 
				   error_message),
	       error_message,error_message);
  
  //Get integrated quantities
  lasagna_call(get_integrated_quantities(y,
					 pqke,
					 &I_VxPy_minus,
					 &I_f0Pa_plus,
					 &I_rho_ss,
					 &I_rho_ss_bar,
					 &I_f0,
					 error_message),
	       error_message,error_message);
   
  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(pqke->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */  
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  dLdT = -1.0/(8.0*H*T*_ZETA3_)*I_VxPy_minus;
  
  //Get partial derivatives:
  lasagna_call(get_partial_derivatives(T,
				       L,
				       dLdT,
				       pqke, 
				       error_message),
	       error_message,error_message);

  /** Calculate RHS: */
  dy[pqke->index_L] = dLdT/_L_SCALE_;
  //Solving mu from L, using the chebyshev cubic root:
  mu_div_T = -2*_PI_/sqrt(3.0)*
    sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
  
  /** All quantities defined on the grid: */
  delta_v = v_grid[1]-v_grid[0];
  for (i=0; i<pqke->vres; i++){
    x = pqke->x_grid[i];
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
    D = 0.5*Gamma;
    
    Pa_plus = y[pqke->index_Pa_plus+i];
    Pa_minus = y[pqke->index_Pa_minus+i];
    Ps_plus = y[pqke->index_Ps_plus+i];
    Ps_minus = y[pqke->index_Ps_minus+i];
    Px_plus = y[pqke->index_Px_plus+i];
    Px_minus = y[pqke->index_Px_minus+i];
    Py_plus = y[pqke->index_Py_plus+i];
    Py_minus = y[pqke->index_Py_minus+i];
  
    //Regulator for sterile population:
    rs = pqke->rs;
    //Distributions:
    feq = 1.0/(1.0+exp(x-mu_div_T));
    feq_bar = 1.0/(1.0+exp(x+mu_div_T));
    feq_plus = feq + feq_bar;
    //Use an expansion for feq_minus since mu_div_T is very small.
    feq_minus = exp(x)*2*mu_div_T/(1+exp(x-mu_div_T))/(1+exp(x+mu_div_T));
    //Use the unexpanded expression if the error grows too large.
    if(exp(x+mu_div_T)/6*pow(mu_div_T,3)/(1+exp(x-mu_div_T))/
       (1+exp(x+mu_div_T))/feq_minus > 1e-10)
      feq_minus = feq - feq_bar;

    f0 = 1.0/(1.0+exp(x));

    dudTdvdu = dudT_grid[i]*dvdu_grid[i];


    //Define index steps for calculating derivatives:
    if (i==0)
      stencil_method = 12; //First order, forward   12
    else if (i==pqke->vres-1)
      stencil_method = 10; //First order, backwards 10
    else if ((i==pqke->vres-2)||(i==1))
      stencil_method = 21; //Second order, centered 21
    else
      stencil_method = 21; //Fifth order, centered  51
    

    idx = pqke->index_Pa_plus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_plus+Gamma*(2.0*feq_plus/f0-Pa_plus))+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = pqke->index_Pa_minus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_minus+Gamma*(2.0*feq_minus/f0-Pa_minus))+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = pqke->index_Ps_plus+i;
    dy[idx] = 1.0/(H*T)*(Vx*Py_plus - 
			 rs*Gamma*(1.0/(4.0*I_f0)*I_rho_ss*feq + 
				   1.0/(4.0*I_f0)*I_rho_ss_bar*feq_bar -
				   0.5*f0*Ps_plus)) +
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = pqke->index_Ps_minus+i;
    dy[idx] = 1.0/(H*T)*(Vx*Py_minus -
			 rs*Gamma*(1.0/(4.0*I_f0)*I_rho_ss*feq - 
				   1.0/(4.0*I_f0)*I_rho_ss_bar*feq_bar -
				   0.5*f0*Ps_minus)) +
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = pqke->index_Px_plus+i;
    dy[idx] = 1.0/(H*T)*((V0+V1)*Py_plus+VL*Py_minus+D*Px_plus)+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = pqke->index_Px_minus+i;
    dy[idx] = 1.0/(H*T)*((V0+V1)*Py_minus+VL*Py_plus+D*Px_minus)+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = pqke->index_Py_plus+i;
    dy[idx] = 1.0/(H*T)*
      (-(V0+V1)*Px_plus-VL*Px_minus+0.5*Vx*(Pa_plus-Ps_plus)+D*Py_plus)+ 
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);
   
    idx = pqke->index_Py_minus+i;
    dy[idx] = 1.0/(H*T)*
      (-(V0+V1)*Px_minus-VL*Px_plus+0.5*Vx*(Pa_minus-Ps_minus)+D*Py_minus)+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);
  }
  return _SUCCESS_;
}

int qke_derivs_fixed_grid(double T, 
			  double *y, 
			  double *dy, 
			  void *param,
			  ErrorMsg error_message){
  qke_param *pqke=param;
  double L;
  int i;
  double *x_grid=pqke->x_grid;
  double gentr,H;
  double Vx, VL;
  double x, w_trapz;
  double Gamma, D, V0, V1, n_plus;
  double f0,  mu_div_T;
  double Pa_plus,Pa_minus,Ps_plus,Ps_minus;
  double Px_plus,Px_minus,Py_plus,Py_minus;
  double feq_plus, feq_minus, feq, feq_bar;
  double I_VxPy_minus, I_f0Pa_plus, I_rho_ss, I_rho_ss_bar, I_f0;
  double rs;
  int idx;
  double dLdT;
 
  L = y[pqke->index_L]*_L_SCALE_;

  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(pqke->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);

  lasagna_call(get_integrated_quantities(y,
					 pqke,
					 &I_VxPy_minus,
					 &I_f0Pa_plus,
					 &I_rho_ss,
					 &I_rho_ss_bar,
					 &I_f0,
					 error_message),
	       error_message,error_message);

  n_plus = I_f0Pa_plus/(3.0*_ZETA3_); 
  //Set V1:
  if (pqke->is_electron==_TRUE_)
    pqke->g_alpha = 1.0+4.0/((1.0-_SIN2_THETA_W_)*n_plus);
  else
    pqke->g_alpha = 1.0;
  pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
    _G_F_/_M_Z_/_M_Z_*pow(T,5)*n_plus*pqke->g_alpha;

  /** Calculate RHS: */ 
  dLdT = -1.0/(8.0*H*T*_ZETA3_)*I_VxPy_minus;
  dy[pqke->index_L] = dLdT/_L_SCALE_;
  //Solving mu from L, using the chebyshev cubic root:
  mu_div_T = -2*_PI_/sqrt(3.0)*
    sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
  
  /** All quantities defined on the grid: */
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
  
    Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
    D = 0.5*Gamma;
    
    Pa_plus = y[pqke->index_Pa_plus+i];
    Pa_minus = y[pqke->index_Pa_minus+i];
    Ps_plus = y[pqke->index_Ps_plus+i];
    Ps_minus = y[pqke->index_Ps_minus+i];
    Px_plus = y[pqke->index_Px_plus+i];
    Px_minus = y[pqke->index_Px_minus+i];
    Py_plus = y[pqke->index_Py_plus+i];
    Py_minus = y[pqke->index_Py_minus+i];
   
    //Regulator for sterile population:
    rs = pqke->rs;
    //Distributions:
    feq_bar = 1.0/(1.0+exp(x+mu_div_T));
    feq_plus = feq + feq_bar;
    feq_plus = feq + feq_bar;
    //Use an expansion for feq_minus since mu_div_T is very small.
    feq_minus = exp(x)*2*mu_div_T/(1+exp(x-mu_div_T))/(1+exp(x+mu_div_T));
    //Use the unexpanded expression if the error grows too large.
    if(exp(x+mu_div_T)/6*pow(mu_div_T,3)/(1+exp(x-mu_div_T))/
       (1+exp(x+mu_div_T))/feq_minus > 1e-10)
      feq_minus = feq - feq_bar;

    f0 = 1.0/(1.0+exp(x));
    
    idx = pqke->index_Pa_plus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_plus+Gamma*(2.0*feq_plus/f0-Pa_plus));

    idx = pqke->index_Pa_minus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_minus+Gamma*(2.0*feq_minus/f0-Pa_minus));

    idx = pqke->index_Ps_plus+i;
    dy[idx] = 1.0/(H*T)*(Vx*Py_plus -
			 rs*Gamma*(1.0/(4.0*I_f0)*I_rho_ss*feq + 
				   1.0/(4.0*I_f0)*I_rho_ss_bar*feq_bar -
				   0.5*f0*Ps_plus));

    idx = pqke->index_Ps_minus+i;
    dy[idx] = 1.0/(H*T)*(Vx*Py_minus -
			 rs*Gamma*(1.0/(4.0*I_f0)*I_rho_ss*feq - 
				   1.0/(4.0*I_f0)*I_rho_ss_bar*feq_bar -
				   0.5*f0*Ps_minus));

    idx = pqke->index_Px_plus+i;
    dy[idx] = 1.0/(H*T)*((V0+V1)*Py_plus+VL*Py_minus+D*Px_plus);

    idx = pqke->index_Px_minus+i;
    dy[idx] = 1.0/(H*T)*((V0+V1)*Py_minus+VL*Py_plus+D*Px_minus);

    idx = pqke->index_Py_plus+i;
    dy[idx] = 1.0/(H*T)*
      (-(V0+V1)*Px_plus-VL*Px_minus+0.5*Vx*(Pa_plus-Ps_plus)+D*Py_plus);
    
    idx = pqke->index_Py_minus+i;
    dy[idx] = 1.0/(H*T)*
      (-(V0+V1)*Px_minus-VL*Px_plus+0.5*Vx*(Pa_minus-Ps_minus)+D*Py_minus);
    
  }
  return _SUCCESS_;
}

double drhodv(double *rho, double delta_v, int index, int stencil_method){
  double drho;
  if (stencil_method == 12)
    drho = (rho[index+1]-rho[index])/delta_v;
  else if (stencil_method == 10)
    drho = (rho[index]-rho[index-1])/delta_v;
  else if (stencil_method == 21)
    drho = (rho[index+1]-rho[index-1])/(2.0*delta_v);
  else
    drho = (-rho[index+2]+8.0*rho[index+1]
	    -8.0*rho[index-1]+rho[index-2])/(12*delta_v);
  return drho;
}


int qke_derivs_test_partial(double T, 
	       double *y, 
	       double *dy, 
	       void *param,
	       ErrorMsg error_message){
  double L,dLdT, delta_v,x,dudTdvdu,Ltau=1e3;
  int i, stencil_method;
  qke_param *pqke = param;
  L = 2e-6*pow(sin(T*Ltau),2);

  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);
  pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
    _G_F_/_M_Z_/_M_Z_*pow(T,5)*2.0;

  get_resonances_xi(T,L,pqke);

  //Get new x_grid
  lasagna_call(get_parametrisation(T,
				   pqke, 
				   error_message),
	       error_message,error_message);
  

  dLdT = 2e-6*2.0*sin(T*Ltau)*cos(T*Ltau)*Ltau;
  //Get partial derivatives:
  lasagna_call(get_partial_derivatives(T,
				       L,
				       dLdT,
				       pqke, 
				       error_message),
	       error_message,error_message);

  /** All quantities defined on the grid: */
  delta_v = pqke->v_grid[1]-pqke->v_grid[0];
  fprintf(stderr,"%.16e ",T);
  for (i=0; i<pqke->vres; i++){
    x = pqke->x_grid[i];
    dudTdvdu = pqke->dudT_grid[i]*pqke->dvdu_grid[i];

    //Define index steps for calculating derivatives:
    if (i==0)
      stencil_method = 12; //First order, forward   12
    else if (i==pqke->vres-1)
      stencil_method = 10; //First order, backwards 10
    else if ((i==pqke->vres-2)||(i==1))
      stencil_method = 21; //Second order, centered 21
    else
      stencil_method = 51; //Fifth order, centered  51

    dy[i] = dudTdvdu*drhodv(y, delta_v, i, stencil_method);
    fprintf(stderr,"%.16e %.16e %.16e ",-dudTdvdu,x,pqke->dudT_grid[i]);
    //fprintf(stderr,"%.16e ",drhodv(y, delta_v, i, stencil_method));
  }
  fprintf(stderr,"\n");
  //fprintf(stderr,"%g %g %g\n",T, dy[pqke->vres-2],dy[pqke->vres-1]);
  
  /**fprintf(stderr,"%g %g %g\n",T, 
	  drhodv(y, delta_v, pqke->vres-2, 10),
	  drhodv(y, delta_v, pqke->vres-1, 10));*/
  return _SUCCESS_;
}
  
int qke_test_partial_output(double T,
			    double *y,
			    double *dy,
			    int index_t,
			    void *param,
			    ErrorMsg error_message){

  FILE *output_file;
  qke_param *pqke=param;
  int i;
  double L,Ltau = 1e3;
  L = 1e-4*pow(sin(T*Ltau),2);

  output_file = fopen(pqke->output_filename,"a+");
  
  fprintf(output_file,"%.16e %.16e ",T*1e3, L);
  for (i=0; i<pqke->vres; i++){
    fprintf(output_file,"%.16e %.16e ",pqke->x_grid[i],y[i]);
  }
  fprintf(output_file,"\n");
  fclose(output_file);

  return _SUCCESS_;
}
