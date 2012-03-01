/** @file qke_equations.c 
 * Thomas Tram Bülow, 25.08.2011    
 */
#include "common.h"
#include "qke_equations.h"
int init_qke_param(qke_param *pqke){
  int i,j,k,idx,neq,nz;
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
  pqke->sigma0Pz = malloc(sizeof(double)*vres);
  pqke->sigma0Pz_bar = malloc(sizeof(double)*vres);
  pqke->P0_ini =  malloc(sizeof(double)*vres);
  pqke->P0_bar_ini =  malloc(sizeof(double)*vres);
  pqke->delta_Pz_ini =  malloc(sizeof(double)*vres);
  pqke->delta_Pz_bar_ini =  malloc(sizeof(double)*vres);
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
  if (pqke->is_electron == _TRUE_)
    pqke->C_alpha = 1.27;
  else
    pqke->C_alpha = 0.92;
  
  
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
  if (pqke->evolve_vi == _TRUE_){
    //    pqke->index_b = idx;
    //    idx++;
    pqke->index_vi = idx;
    idx += pqke->Nres;
  }

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
    //F-Ps_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[pqke->index_Ps_minus+j][pqke->index_Ps_minus+k] = 1;
    J[pqke->index_Ps_minus+j][pqke->index_Py_minus+j] = 1;
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

int init_qke_param_fixed_grid_low_temp(qke_param *pqke){
  int i,j,k,idx,neq,nz;
  double k1,k2;
  double Nres, vres, Tres;
  int **J;
  Nres = pqke->Nres;
  vres = pqke->vres;
  Tres = pqke->Tres;
  
  free(pqke->Ap);
  free(pqke->Ai);
  
  //Set up the indices:
  idx = 0;
  pqke->index_L = idx;
  idx++;
  pqke->index_n_plus = idx;
  idx++;
  pqke->index_Q0_plus = idx;
  idx+=pqke->vres;
  pqke->index_Q0_minus =idx;
  idx +=pqke->vres;
  pqke->index_Q1_plus = idx;
  idx +=pqke->vres;
  pqke->index_Q1_minus = idx;
  idx +=pqke->vres;
  pqke->index_Q2_plus = idx;
  idx +=pqke->vres;
  pqke->index_Q2_minus = idx;
  idx +=pqke->vres;
  pqke->index_Q3_plus = idx;
  idx +=pqke->vres;
  pqke->index_Q3_minus = idx;
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
  //Everything couples to integrated quantities n_plus and L through Vz:
  for (i=0; i<neq; i++)
    J[i][pqke->index_L] = 1;
  for (i=0; i<neq; i++)
    J[i][pqke->index_n_plus] = 1;

  //F-L dependence:
  for (i=0; i<vres; i++)
    J[pqke->index_L][pqke->index_Q2_minus+i] = 1;
  //F-n_plus dependence:
  for (i=0; i<vres; i++){
    J[pqke->index_n_plus][pqke->index_Q0_plus+i] = 1;
    J[pqke->index_n_plus][pqke->index_Q1_plus+i] = 1;
    J[pqke->index_n_plus][pqke->index_Q3_plus+i] = 1;
    J[pqke->index_n_plus][pqke->index_Q2_plus+i] = 1;
  }
  //Loop over grid:
  for (j=0; j<vres; j++){
    //F-Q0_plus[j] dependence:
    J[pqke->index_Q0_plus+j][pqke->index_Q1_plus+j] = 1;
    J[pqke->index_Q0_plus+j][pqke->index_Q3_plus+j] = 1;
    //F-Q0_minus[j] dependence:
    J[pqke->index_Q0_minus+j][pqke->index_Q1_minus+j] = 1;
    J[pqke->index_Q0_minus+j][pqke->index_Q3_minus+j] = 1;
    //F-Q1_plus[j] dependence:
    J[pqke->index_Q1_plus+j][pqke->index_Q3_plus+j] = 1;
    J[pqke->index_Q1_plus+j][pqke->index_Q2_plus+j] = 1;
    J[pqke->index_Q1_plus+j][pqke->index_Q2_minus+j] = 1;
    J[pqke->index_Q1_plus+j][pqke->index_Q1_plus+j] = 1;
    J[pqke->index_Q1_plus+j][pqke->index_Q0_plus+j] = 1;
    //F-Q1_minus[j] dependence:
    J[pqke->index_Q1_minus+j][pqke->index_Q3_minus+j] = 1;
    J[pqke->index_Q1_minus+j][pqke->index_Q2_minus+j] = 1;
    J[pqke->index_Q1_minus+j][pqke->index_Q2_plus+j] = 1;
    J[pqke->index_Q1_minus+j][pqke->index_Q1_minus+j] = 1;
    J[pqke->index_Q1_minus+j][pqke->index_Q0_minus+j] = 1;
    //F-Q2_plus[j] dependence:
    J[pqke->index_Q2_plus+j][pqke->index_Q1_plus+j] = 1;
    J[pqke->index_Q2_plus+j][pqke->index_Q1_minus+j] = 1;
    J[pqke->index_Q2_plus+j][pqke->index_Q3_minus+j] = 1;
    //F-Q2_minus[j] dependence:
    J[pqke->index_Q2_minus+j][pqke->index_Q1_minus+j] = 1;
    J[pqke->index_Q2_minus+j][pqke->index_Q1_plus+j] = 1;
    J[pqke->index_Q2_minus+j][pqke->index_Q3_plus+j] = 1;
    //F-Q3_plus[j] dependence:
    J[pqke->index_Q3_plus+j][pqke->index_Q1_plus+j] = 1;
    J[pqke->index_Q3_plus+j][pqke->index_Q2_minus+j] = 1;
    J[pqke->index_Q3_plus+j][pqke->index_Q0_plus+j] = 1;
    //F-Q3_minus[j] dependence:
    J[pqke->index_Q3_minus+j][pqke->index_Q1_minus+j] = 1;
    J[pqke->index_Q3_minus+j][pqke->index_Q2_plus+j] = 1;
    J[pqke->index_Q3_minus+j][pqke->index_Q0_minus+j] = 1;
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

int init_qke_param_fixed_grid_standard(qke_param *pqke){
  int i,j,k,idx,neq,nz;
  double k1,k2;
  double Nres, vres, Tres;
  int **J;
  Nres = pqke->Nres;
  vres = pqke->vres;
  Tres = pqke->Tres;
  
  free(pqke->Ap);
  free(pqke->Ai);
  
  //Set up the indices:
  idx = 0;
  pqke->index_L = idx;
  idx++;
  pqke->index_n_plus = idx;
  idx++;
  pqke->index_P0 = idx;
  idx+=pqke->vres;
  pqke->index_P0_bar =idx;
  idx +=pqke->vres;
  pqke->index_Px = idx;
  idx +=pqke->vres;
  pqke->index_Px_bar = idx;
  idx +=pqke->vres;
  pqke->index_Py = idx;
  idx +=pqke->vres;
  pqke->index_Py_bar = idx;
  idx +=pqke->vres;
  pqke->index_Pz = idx;
  idx +=pqke->vres;
  pqke->index_Pz_bar = idx;
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
  //Everything couples to integrated quantities n_plus and L through Vz:
  for (i=0; i<neq; i++)
    J[i][pqke->index_L] = 1;
  for (i=0; i<neq; i++)
    J[i][pqke->index_n_plus] = 1;

  //F-L dependence:
  for (i=0; i<vres; i++){
    J[pqke->index_L][pqke->index_Py+i] = 1;
    J[pqke->index_L][pqke->index_Py_bar+i] = 1;
  }
  //F-n_plus dependence:
  for (i=0; i<vres; i++){
    J[pqke->index_n_plus][pqke->index_P0+i] = 1;
    J[pqke->index_n_plus][pqke->index_P0_bar+i] = 1;
    J[pqke->index_n_plus][pqke->index_Pz+i] = 1;
    J[pqke->index_n_plus][pqke->index_Pz_bar+i] = 1;
    J[pqke->index_n_plus][pqke->index_Py+i] = 1;
    J[pqke->index_n_plus][pqke->index_Py_bar+i] = 1;
  } 
 //Loop over grid:
  for (j=0; j<vres; j++){
    //F-P0[j] dependence:
    J[pqke->index_P0+j][pqke->index_Pz+j] = 1;
    //F-P0_bar[j] dependence:
    J[pqke->index_P0_bar+j][pqke->index_Pz_bar+j] = 1;
    //F-Px[j] dependence:
    J[pqke->index_Px+j][pqke->index_Py+j] = 1;
    //F-Px_bar[j] dependence:
    J[pqke->index_Px_bar+j][pqke->index_Py_bar+j] = 1;
    //F-Py[j] dependence:
    J[pqke->index_Py+j][pqke->index_Px+j] = 1;
    J[pqke->index_Py+j][pqke->index_Pz+j] = 1;
    //F-Py_bar[j] dependence:
    J[pqke->index_Py_bar+j][pqke->index_Px_bar+j] = 1;
    J[pqke->index_Py_bar+j][pqke->index_Pz_bar+j] = 1;
    //F-Pz[j] dependence:
    J[pqke->index_Pz+j][pqke->index_Py+j] = 1;
    J[pqke->index_Pz+j][pqke->index_P0+j] = 1;
    //F-Pz_bar[j] dependence:
    J[pqke->index_Pz_bar+j][pqke->index_Py_bar+j] = 1;
    J[pqke->index_Pz_bar+j][pqke->index_P0_bar+j] = 1;
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

int init_qke_param_fixed_grid_standard_approx(qke_param *pqke){
  int i,j,k,idx,neq,nz;
  double k1,k2;
  double Nres, vres, Tres;
  int **J;
  Nres = pqke->Nres;
  vres = pqke->vres;
  Tres = pqke->Tres;
  
  free(pqke->Ap);
  free(pqke->Ai);
  
  //Set up the indices:
  idx = 0;
  //  pqke->index_L = idx;
  //idx++;
  pqke->index_P0 = idx;
  idx +=pqke->vres;
  pqke->index_P0_bar =idx;
  idx +=pqke->vres;
  pqke->index_Pz = idx;
  idx +=pqke->vres;
  pqke->index_Pz_bar = idx;
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
  //Everything couples to L..
  for (i=0; i<neq; i++){
    J[i][pqke->index_L] = 1;
  }
  
  //F-L dependence:
  //L is coupled to itself, and also n+nbar, but we neglect the latter.

 //Loop over grid:
  for (j=0; j<vres; j++){
    //F-P0[j] dependence:
    J[pqke->index_P0+j][pqke->index_Pz+j] = 1;
    //F-P0_bar[j] dependence:
    J[pqke->index_P0_bar+j][pqke->index_Pz_bar+j] = 1;
    //F-Pz[j] dependence:
    J[pqke->index_Pz+j][pqke->index_P0+j] = 1;
    for (i=0; i<vres; i++){
      J[pqke->index_Pz+j][pqke->index_Pz+i] = 1;
      J[pqke->index_Pz+j][pqke->index_Pz_bar+i] = 1;
    }
    //F-Pz_bar[j] dependence:
    J[pqke->index_Pz_bar+j][pqke->index_P0_bar+j] = 1;
    for (i=0; i<vres; i++){
      J[pqke->index_Pz_bar+j][pqke->index_Pz+i] = 1;
      J[pqke->index_Pz_bar+j][pqke->index_Pz_bar+i] = 1;
    }
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

int init_qke_param_fixed_grid_standard_approx2(qke_param *pqke){
  int i,j,k,idx,neq,nz;
  double k1,k2;
  double Nres, vres, Tres;
  int **J;
  Nres = pqke->Nres;
  vres = pqke->vres;
  Tres = pqke->Tres;
  
  //Set up the indices:
  idx = 0;
  pqke->index_P0 = idx;
  idx +=pqke->vres;
  pqke->index_P0_bar =idx;
  idx +=pqke->vres;

  //Last thing to do:
  neq = idx;
  pqke->neq = neq;

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
  free(pqke->sigma0Pz);
  free(pqke->sigma0Pz_bar);
  free(pqke->P0_ini);
  free(pqke->P0_bar_ini);
  free(pqke->delta_Pz_ini);
  free(pqke->delta_Pz_bar_ini);
  for (i=0; i<(pqke->Nres+2); i++) 
    free(pqke->mat[i]);
  free(pqke->mat);
  free(pqke->vv);
  free(pqke->indx);
  free(pqke->Ap);
  free(pqke->Ai);

  return _SUCCESS_;
};


int get_resonances_xi(double T, double L, qke_param *param){
	double phi,chi;
	double *xi=param->xi;
	double *dxidT = param->dxidT;
	double frac;
	int i;

	//Make sure the resonances are sorted:
	L = fabs(L);

	if (param->is_electron == _FALSE_){
	  phi = 45.0*_M_Z_*_M_Z_*_ZETA3_*L/(7.0*T*T*pow(_PI_,4));
	  chi = 45.0*_M_Z_*_M_Z_*fabs(param->delta_m2)*
	    cos(2.0*param->theta_zero)/
	    (14.0*sqrt(2.0)*_PI_*_PI_*_G_F_*pow(T,6));
	}
	else{
	  phi = 45.0*_M_W_*_M_W_*_ZETA3_*L/
	    (14.0*(1.5-0.5*_SIN2_THETA_W_)*T*T*pow(_PI_,4));
	  chi = 45.0*fabs(param->delta_m2)*cos(2.0*param->theta_zero)*_M_W_*_M_W_/
	    (28.0*sqrt(2.0)*_PI_*_PI_*_G_F_*pow(T,6)*(1.5-0.5*_SIN2_THETA_W_));
	}
	xi[0] = sqrt(phi*phi+chi)-phi;
	xi[1] = sqrt(phi*phi+chi)+phi;
	
	/**	if (xi[1]<xi[0])
	  printf("x1 = %g, x2 = %g, ???\n",xi[0],xi[1]);
	*/

	frac = (2.0*phi*phi+3*chi)/sqrt(phi*phi+chi);
	dxidT[0] = -(frac+2.0*phi)/T;
	dxidT[1] = -(frac-2.0*phi)/T;

	// Protect against possible seg fault:
	for (i=0; i<param->Nres; i++){
	  if ((1==0)&&(xi[i]<param->xmin)) {
	    /**printf("Note: Resonance at T=%g MeV is lower than xmin=%g. (It is %g.)\n",
	       T*1e3,param->xmin,xi[i]);*/
	    xi[i] = param->xmin;
	    dxidT[i] = 0.0;
	  }
	  else if(xi[i]>param->xmax){
	    /**	    printf("Note: Resonance at T=%g MeV is higher than xmax=%g. (It is %g.)\n",
		    T*1e3,param->xmax,xi[i]);*/
	    xi[i] = param->xmax;
	    dxidT[i] = 0.0;
	  }
	}

	return _SUCCESS_;
};


int u_of_x(double x, double *u, double *dudx, qke_param *param){
  double x_plus_xext;
  x_plus_xext = x+param->xext;
  *u = (x*(1.0+param->eps2)-param->xext*param->eps1)/x_plus_xext;
  *dudx = param->xext*(1+param->eps1+param->eps2)/x_plus_xext/x_plus_xext;
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

int nonlinear_rhs2(double *y, double *Fy, void *param){
  qke_param *pqke=param;
  double b=y[0];
  double *wi=y+1;
  double sum_wi=0.0;
  double alpha = pqke->alpha;
  double *ui=pqke->ui;
  int i,n=pqke->Nres;

  for (i=0; i<n; i++){
    sum_wi += wi[i];
  }
  
  Fy[0] = ui[0]-alpha*sum_wi-b*pow(sum_wi,3);
  for (i=1; i<n; i++){
    Fy[i] = ui[i-1]-ui[i]-0.25*b*pow(wi[i-1],3)-alpha*wi[i-1];
  }
  Fy[n] = ui[n-1]-1.0+alpha*(1.0-wi[n-1])+b*pow(1.0-wi[n-1],3);

  return _SUCCESS_;
}

int qke_initial_conditions(double Ti, double *y, qke_param *pqke){
  /** Set initial conditions at temperature Ti: */
  int i, evolve_vi;
  ErrorMsg error_message;
  double x, *dy;
  double Vx,D,Vz,Vz_bar,Px,Py,Px_bar,Py_bar;
  /** Calculate chemical potential from initial L:
  mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*pqke->L_initial/pow(_PI_,3)));
  */

  //Assuming y is calloc'ed -- dangerous, better to zero it.
  for (i=0; i<pqke->neq; i++) y[i] = 0.0;

  //Set standard equilibrium initial conditions:
  y[pqke->index_L] = pqke->L_initial/_L_SCALE_;
  for (i=0; i<pqke->vres; i++)
    y[pqke->index_Pa_plus+i] = 4.0;
  
  //Do one call to qke_derivs to get the coordinate transformation.
  evolve_vi = pqke->evolve_vi;
  //We must always solve the non-linear equation in the first step:
  pqke->evolve_vi = _FALSE_;
  dy = malloc(sizeof(double)*pqke->neq);
  qke_derivs(Ti,y,dy,pqke,error_message);
  pqke->evolve_vi = evolve_vi;
  if (pqke->evolve_vi == _TRUE_){
    //We must set initial conditions for these:
    //y[pqke->index_b] = pqke->b;
    for (i=0; i<pqke->Nres; i++)
      y[pqke->index_vi+i] = pqke->vi[i];
  }
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
    /**
       Old code for setting initial condition for non-zero L_nu:
       if(_FALSE_){ //(mu_div_T < 0.05){
       Use series expansion in mu_div_T to set initial values:
       y[pqke->index_Pa_plus+i] = 4.0+
       2.0*exp(x)*(exp(x)-1.0)/pow(exp(x)+1.0,2)*pow(mu_div_T,2);
       printf("Pa_plus initial: %e\n",y[pqke->index_Pa_plus+i]);
       y[pqke->index_Pa_minus+i] = 4.0*exp(x)/(exp(x)+1.0)*mu_div_T;
       }
       else{
       y[pqke->index_Pa_plus+i] = 2.0*(1.0+exp(x))*
       (1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T)));
       y[pqke->index_Pa_minus+i] = 2.0*(1.0+exp(x))*
       (1.0/(1.0+exp(x-mu_div_T))-1.0/(1.0+exp(x+mu_div_T)));
       }
    */
  free(dy);
  return _SUCCESS_;
};

int qke_initial_conditions_fixed_grid(double Ti, double *y, qke_param *pqke){
  /** Set initial conditions at temperature Ti: */
  int i, evolve_vi;
  ErrorMsg error_message;
  double x, *dy;
  double Vx,D,Vz,Vz_bar,Px,Py,Px_bar,Py_bar;
  /** Calculate chemical potential from initial L:
  mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*pqke->L_initial/pow(_PI_,3)));
  */

  //Assuming y is calloc'ed -- dangerous, better to zero it.
  for (i=0; i<pqke->neq; i++) y[i] = 0.0;

  //Set standard equilibrium initial conditions:
  y[pqke->index_L] = pqke->L_initial/_L_SCALE_;
  for (i=0; i<pqke->vres; i++)
    y[pqke->index_Pa_plus+i] = 4.0;
  
  //Do one call to qke_derivs to get the coordinate transformation.
  evolve_vi = pqke->evolve_vi;
  //We must always solve the non-linear equation in the first step:
  pqke->evolve_vi = _FALSE_;
  dy = malloc(sizeof(double)*pqke->neq);
  qke_derivs(Ti,y,dy,pqke,error_message);

  //Debug stuff:
  /**  for (i=0; i<40; i++){
    pqke->x_grid[i] = max(1e-10,pqke->xi[0]+(-20+i)*0.2*pqke->xi[0]*pqke->theta_zero);
  }
  */
  pqke->evolve_vi = evolve_vi;
  if (pqke->evolve_vi == _TRUE_){
    //We must set initial conditions for these:
    //y[pqke->index_b] = pqke->b;
    for (i=0; i<pqke->Nres; i++)
      y[pqke->index_vi+i] = pqke->vi[i];
  }
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
    /**
       Old code for setting initial condition for non-zero L_nu:
       if(_FALSE_){ //(mu_div_T < 0.05){
       Use series expansion in mu_div_T to set initial values:
       y[pqke->index_Pa_plus+i] = 4.0+
       2.0*exp(x)*(exp(x)-1.0)/pow(exp(x)+1.0,2)*pow(mu_div_T,2);
       printf("Pa_plus initial: %e\n",y[pqke->index_Pa_plus+i]);
       y[pqke->index_Pa_minus+i] = 4.0*exp(x)/(exp(x)+1.0)*mu_div_T;
       }
       else{
       y[pqke->index_Pa_plus+i] = 2.0*(1.0+exp(x))*
       (1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T)));
       y[pqke->index_Pa_minus+i] = 2.0*(1.0+exp(x))*
       (1.0/(1.0+exp(x-mu_div_T))-1.0/(1.0+exp(x+mu_div_T)));
       }
    */
  free(dy);
  return _SUCCESS_;
};

int qke_copy_IC_fixed_grid_to_fixed_grid_low_temp(double T,
						  double *y_in,
						  double *y_out,
						  qke_param *pqke){
  int i;
  double *x_grid = pqke->x_grid;
  double w_trapz,x,f0,Pa_plus,I_f0Pa_plus=0.0;
  double x2,Vx,V0,V1,w,s;
  double Px_plus,Px_minus,Py_plus,Py_minus,Pz_plus,Pz_minus;
  /** This is not so nice, since I should not require the index
      to be unchanged from one system to the other. However,
      this is only an issue for L. */
  y_out[pqke->index_L] = y_in[pqke->index_L];
  //Calculate n+nbar:
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    f0 = 1.0/(1.0+exp(x));
    Pa_plus = y_in[pqke->index_Pa_plus+i];
    I_f0Pa_plus += w_trapz*(x*x*f0*Pa_plus);
  }
  y_out[pqke->index_n_plus] = I_f0Pa_plus/(3.0*_ZETA3_);
  printf("n_plus=%g\n",y_out[pqke->index_n_plus]);

  //Loop over grid:
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];
    x2 = x*x;
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    w = (V0+V1)/Vx;
    s = sqrt(1.0+w*w);

    Px_plus = y_in[pqke->index_Px_plus+i];
    Px_minus = y_in[pqke->index_Px_minus+i];
    Py_plus = y_in[pqke->index_Py_plus+i];
    Py_minus = y_in[pqke->index_Py_minus+i];
    Pz_plus = 0.5*(y_in[pqke->index_Pa_plus+i]-y_in[pqke->index_Ps_plus+i]);
    Pz_minus = 0.5*(y_in[pqke->index_Pa_minus+i]-y_in[pqke->index_Ps_minus+i]);


    y_out[pqke->index_Q0_plus+i] = 
      0.5*(y_in[pqke->index_Pa_plus+i]+y_in[pqke->index_Ps_plus+i]);
    y_out[pqke->index_Q0_minus+i] = 
      0.5*(y_in[pqke->index_Pa_minus+i]+y_in[pqke->index_Ps_minus+i]);
    y_out[pqke->index_Q1_plus+i] = -w*Px_plus+Pz_plus;
    y_out[pqke->index_Q1_minus+i] = -w*Px_minus+Pz_minus;
    y_out[pqke->index_Q2_plus+i] = s*Py_plus;
    y_out[pqke->index_Q2_minus+i] = s*Py_minus;
    y_out[pqke->index_Q3_plus+i] = w*Px_plus+w*w*Pz_plus;
    y_out[pqke->index_Q3_minus+i] = w*Px_minus+w*w*Pz_minus;
    
  }
  return _SUCCESS_;
}

int qke_copy_IC_fixed_grid_to_fixed_grid_standard(double T,
						  double *y_in,
						  double *y_out,
						  qke_param *pqke){
  int i;
  double *x_grid = pqke->x_grid;
  double w_trapz,x,f0,Pa_plus,I_f0Pa_plus=0.0;
  double x2,Vx,V0,V1,w,s;
  double Px_plus,Px_minus,Py_plus,Py_minus,Pz_plus,Pz_minus,P0_plus,P0_minus;
  /** This is not so nice, since I should not require the index
      to be unchanged from one system to the other. However,
      this is only an issue for L. */
  y_out[pqke->index_L] = y_in[pqke->index_L];

  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    f0 = 1.0/(1.0+exp(x));
    Pa_plus = y_in[pqke->index_Pa_plus+i];
    I_f0Pa_plus += w_trapz*(x*x*f0*Pa_plus);
  }
  y_out[pqke->index_n_plus] = I_f0Pa_plus/(3.0*_ZETA3_);
  printf("n_plus=%g\n",y_out[pqke->index_n_plus]);

  //Loop over grid:
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];
    x2 = x*x;
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    w = (V0+V1)/Vx;
    s = sqrt(1.0+w*w);

    Px_plus = y_in[pqke->index_Px_plus+i];
    Px_minus = y_in[pqke->index_Px_minus+i];
    Py_plus = y_in[pqke->index_Py_plus+i];
    Py_minus = y_in[pqke->index_Py_minus+i];
    Pz_plus = 0.5*(y_in[pqke->index_Pa_plus+i]-y_in[pqke->index_Ps_plus+i]);
    Pz_minus = 0.5*(y_in[pqke->index_Pa_minus+i]-y_in[pqke->index_Ps_minus+i]);
    P0_plus = 0.5*(y_in[pqke->index_Pa_plus+i]+y_in[pqke->index_Ps_plus+i]);
    P0_minus = 0.5*(y_in[pqke->index_Pa_minus+i]+y_in[pqke->index_Ps_minus+i]);


    y_out[pqke->index_P0+i] = 0.5*(P0_plus+P0_minus);
    y_out[pqke->index_P0_bar+i] = 0.5*(P0_plus-P0_minus);
    y_out[pqke->index_Px+i] = 0.5*(Px_plus+Px_minus);
    y_out[pqke->index_Px_bar+i] = 0.5*(Px_plus-Px_minus);
    y_out[pqke->index_Py+i] = 0.5*(Py_plus+Py_minus);
    y_out[pqke->index_Py_bar+i] = 0.5*(Py_plus-Py_minus);
    y_out[pqke->index_Pz+i] = 0.5*(Pz_plus+Pz_minus);
    y_out[pqke->index_Pz_bar+i] = 0.5*(Pz_plus-Pz_minus);
    
  }
  return _SUCCESS_;
}

int qke_copy_IC_fixed_grid_standard_to_fixed_grid_standard_approx(double T,
								  double *y_in,
								  double *y_out,
								  qke_param *pqke){
  int i;
  double *x_grid = pqke->x_grid;
  double w_trapz,x,f0,Pa_plus,I_f0Pa_plus=0.0;
  double x2,Vx,V0,V1,w,s,Vz,VL,Vz_bar;
  double P0,P0_bar,Pz,Pz_bar,sin_alpha,sin_alpha_bar;
  double Px_plus,Px_minus,Py_plus,Py_minus,Pz_plus,Pz_minus,P0_plus,P0_minus;
  /** This is not so nice, since I should not require the index
      to be unchanged from one system to the other. However,
      this is only an issue for L. */
  double I_Paminus,L,L_int;

  L = y_in[pqke->index_L]*_L_SCALE_;
  VL = pqke->VL;
  I_Paminus = 0.0;
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];

    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    f0 = 1.0/(exp(x)+1.0);

    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Vz = V0+V1+VL;
    Vz_bar = V0+V1-VL;

    P0 = y_in[pqke->index_P0+i+2];
    P0_bar = y_in[pqke->index_P0_bar+i+2];
    if (i==0)
      printf("P0(0) = %g\n",P0);

    //Cheating, since I haven't figured out how to do this consistently.
    Pz = y_in[pqke->index_Py+2*pqke->vres+i];
    Pz_bar = y_in[pqke->index_Py_bar+2*pqke->vres+i];
    
    I_Paminus += w_trapz*x*x*f0*(P0-P0_bar+Pz-Pz_bar);
    
 
    sin_alpha = Vz/sqrt(Vz*Vz+Vx*Vx);
    sin_alpha_bar = Vz_bar/sqrt(Vz_bar*Vz_bar+Vx*Vx);
    
    pqke->sigma0Pz[i] = sin_alpha*Pz;
    pqke->sigma0Pz_bar[i] = sin_alpha_bar*Pz_bar;
    
    y_out[pqke->index_P0+i] = P0;
    y_out[pqke->index_P0_bar+i] = P0_bar;
    y_out[pqke->index_Pz+i] = Pz;
    y_out[pqke->index_Pz_bar+i] = Pz_bar;
  }
  L_int = I_Paminus/(8.0*_ZETA3_);
  pqke->L_fudge = L/L_int;
  printf("Fudge factor for L: %g\n",pqke->L_fudge);

  return _SUCCESS_;
}

int qke_copy_IC_fixed_grid_standard_to_fixed_grid_standard_approx2(
      double T,
      double *y_in,
      double *y_out,
      qke_param *pqke){
  int i;
  double *x_grid = pqke->x_grid;
  double w_trapz,x,f0,Pa_plus,I_f0Pa_plus=0.0;
  double x2,Vx,V0,V1,w,s,Vz,VL,Vz_bar;
  double P0,P0_bar,Pz,Pz_bar,sin_alpha,sin_alpha_bar;
  double Px_plus,Px_minus,Py_plus,Py_minus,Pz_plus,Pz_minus,P0_plus,P0_minus;
  /** This is not so nice, since I should not require the index
      to be unchanged from one system to the other. However,
      this is only an issue for L. */
  double I_Pzpure,L,L_int,I_Paplus;
  double P0_ini,P0_bar_ini,gamma,gamma_bar,delta_Pz,delta_Pz_bar,Pz_pure,Pz_bar_pure;

  L = y_in[pqke->index_L]*_L_SCALE_;
  VL = pqke->VL;
  
  I_Paplus = 0.0;
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];

    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    f0 = 1.0/(exp(x)+1.0);

    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Vz = V0+V1+VL;
    Vz_bar = V0+V1-VL;

    //Cheating, since I haven't figured out how to do this consistently.
    P0 = y_in[pqke->index_P0+i+2];
    P0_bar = y_in[pqke->index_P0_bar+i+2];
    Pz = y_in[pqke->index_Py+2*pqke->vres+i];
    Pz_bar = y_in[pqke->index_Py_bar+2*pqke->vres+i];
    

    sin_alpha = Vz/sqrt(Vz*Vz+Vx*Vx);
    sin_alpha_bar = Vz_bar/sqrt(Vz_bar*Vz_bar+Vx*Vx);
    
    Pz_pure = pow(sin_alpha,2)*Pz;
    Pz_bar_pure = pow(sin_alpha_bar,2)*Pz_bar;
    
    I_Pzpure +=w_trapz*x*x*f0*(Pz_pure-Pz_bar_pure);
    I_Paplus +=w_trapz*x*x*f0*(P0+P0_bar+Pz+Pz_bar);
    
    pqke->delta_Pz_ini[i] = Pz-Pz_pure;
    pqke->delta_Pz_bar_ini[i] = Pz_bar-Pz_bar_pure;
    pqke->P0_ini[i] = P0;
    pqke->P0_bar_ini[i] = P0_bar;
    pqke->sigma0Pz[i] = sin_alpha*Pz;
    pqke->sigma0Pz_bar[i] = sin_alpha_bar*Pz_bar;
    
    y_out[pqke->index_P0+i] = P0;
    y_out[pqke->index_P0_bar+i] = P0_bar;
  }
  pqke->L0 = L-I_Pzpure/(8.0*_ZETA3_);
  pqke->Lz = I_Pzpure/(8.0*_ZETA3_);
  pqke->n_plus = I_Paplus/(3.0*_ZETA3_);
  

  return _SUCCESS_;
}


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
  mat_add_matrix(outf,"xi",miDOUBLE,Tres,Nres,&(pqke->xi_handle));
  mat_add_matrix(outf,"ui",miDOUBLE,Tres,Nres,&(pqke->ui_handle));
  mat_add_matrix(outf,"vi",miDOUBLE,Tres,Nres,&(pqke->vi_handle));
  mat_add_matrix(outf,"b_a_vec",miDOUBLE,Tres,(1+Nres),&(pqke->b_a_vec_handle));
  //Add other matrices:
  mat_add_matrix(outf,"L",miDOUBLE,Tres,1,&(pqke->L_handle));
  mat_add_matrix(outf,"T",miDOUBLE,Tres,1,&(pqke->T_handle));
  mat_add_matrix(outf,"I_conserved",miDOUBLE,Tres,1,&(pqke->I_conserved_handle));
  mat_add_matrix(outf,"V0",miDOUBLE,Tres,1,&(pqke->V0_handle));
  mat_add_matrix(outf,"V1",miDOUBLE,Tres,1,&(pqke->V1_handle));
  mat_add_matrix(outf,"Vx",miDOUBLE,Tres,1,&(pqke->Vx_handle));
  mat_add_matrix(outf,"VL",miDOUBLE,Tres,1,&(pqke->VL_handle));
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

int qke_store_output_fixed_grid_low_temp(double T,
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
  double *out_vec, *Pa_plus_ptr, *Pa_minus_ptr, *Ps_plus_ptr, *Ps_minus_ptr;
  double *Px_plus_ptr, *Px_minus_ptr, *Py_plus_ptr, *Py_minus_ptr;
  double Vx,V0,V1,w,s;
  double Q0_plus,Q0_minus,Q1_plus,Q1_minus,Q2_plus,Q2_minus,Q3_plus,Q3_minus;

  Pa_plus_ptr = malloc(8*sizeof(double)*vres);
  Pa_minus_ptr = Pa_plus_ptr+vres;
  Ps_plus_ptr = Pa_plus_ptr+2*vres;
  Ps_minus_ptr = Pa_plus_ptr+3*vres;
  Px_plus_ptr = Pa_plus_ptr+4*vres;
  Px_minus_ptr = Pa_plus_ptr+5*vres;
  Py_plus_ptr = Pa_plus_ptr+6*vres;
  Py_minus_ptr = Pa_plus_ptr+7*vres;

  for (i=0; i<vres; i++){
    x = pqke->x_grid[i];
   
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    w = (V0+V1)/Vx;
    s = sqrt(1.0+w*w);
  
   
    Q0_plus = y[pqke->index_Q0_plus+i];
    Q1_plus = y[pqke->index_Q1_plus+i];
    Q2_plus = y[pqke->index_Q2_plus+i];
    Q3_plus = y[pqke->index_Q3_plus+i];
    Q0_minus = y[pqke->index_Q0_minus+i];
    Q1_minus = y[pqke->index_Q1_minus+i];
    Q2_minus = y[pqke->index_Q2_minus+i];   
    Q3_minus = y[pqke->index_Q3_minus+i];
   
    Pa_plus_ptr[i] = Q0_plus+(Q1_plus+Q3_plus)/s/s;
    Pa_minus_ptr[i] = Q0_minus+(Q1_minus+Q3_minus)/s/s;
    Ps_plus_ptr[i] = Q0_plus-(Q1_plus+Q3_plus)/s/s;
    Ps_minus_ptr[i] = Q0_minus-(Q1_minus+Q3_minus)/s/s;
    Px_plus_ptr[i] = -w/s/s*Q1_plus+1.0/w/s/s*Q3_plus;
    Px_minus_ptr[i] = -w/s/s*Q1_minus+1.0/w/s/s*Q3_minus;
    Py_plus_ptr[i] = 1.0/s*Q2_plus;
    Py_minus_ptr[i] = 1.0/s*Q2_minus;
  }
    
  //printf("Storing output at index: %d\n",index_t);
  mat_file = fopen(pqke->output_filename,"r+b");
  //Write stuff from y_vector:
  L = y[pqke->index_L]*_L_SCALE_;
  mat_write_fast(&L,pqke->L_handle,8,1);
  mat_write_fast(Pa_plus_ptr,pqke->Pa_plus_handle,8,vres);
  mat_write_fast(Pa_minus_ptr,pqke->Pa_minus_handle,8,vres);
  mat_write_fast(Ps_plus_ptr,pqke->Ps_plus_handle,8,vres);
  mat_write_fast(Ps_minus_ptr,pqke->Ps_minus_handle,8,vres);
  mat_write_fast(Px_plus_ptr,pqke->Px_plus_handle,8,vres);
  mat_write_fast(Px_minus_ptr,pqke->Px_minus_handle,8,vres);
  mat_write_fast(Py_plus_ptr,pqke->Py_plus_handle,8,vres);
  mat_write_fast(Py_minus_ptr,pqke->Py_minus_handle,8,vres);
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
  free(Pa_plus_ptr);
  return _SUCCESS_;
};

int qke_store_output_fixed_grid_standard(double T,
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
  double *out_vec, *Pa_plus_ptr, *Pa_minus_ptr, *Ps_plus_ptr, *Ps_minus_ptr;
  double *Px_plus_ptr, *Px_minus_ptr, *Py_plus_ptr, *Py_minus_ptr;
  double Vx,V0,V1,w,s;
  double P0, P0_bar, Px, Px_bar, Py, Py_bar, Pz, Pz_bar;

  Pa_plus_ptr = malloc(8*sizeof(double)*vres);
  Pa_minus_ptr = Pa_plus_ptr+vres;
  Ps_plus_ptr = Pa_plus_ptr+2*vres;
  Ps_minus_ptr = Pa_plus_ptr+3*vres;
  Px_plus_ptr = Pa_plus_ptr+4*vres;
  Px_minus_ptr = Pa_plus_ptr+5*vres;
  Py_plus_ptr = Pa_plus_ptr+6*vres;
  Py_minus_ptr = Pa_plus_ptr+7*vres;

  for (i=0; i<vres; i++){
    
    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    Px = y[pqke->index_Px+i];
    Px_bar = y[pqke->index_Px_bar+i];
    Py = y[pqke->index_Py+i];
    Py_bar = y[pqke->index_Py_bar+i];
    Pz = y[pqke->index_Pz+i];
    Pz_bar = y[pqke->index_Pz_bar+i];


    Pa_plus_ptr[i] = P0+P0_bar+Pz+Pz_bar;
    Pa_minus_ptr[i] = P0-P0_bar+Pz-Pz_bar;;
    Ps_plus_ptr[i] = P0+P0_bar-(Pz+Pz_bar);
    Ps_minus_ptr[i] = P0-P0_bar-(Pz-Pz_bar);
    Px_plus_ptr[i] = Px+Px_bar;
    Px_minus_ptr[i] = Px-Px_bar;
    Py_plus_ptr[i] = Py+Py_bar;
    Py_minus_ptr[i] = Py-Py_bar;
  }
    
  //printf("Storing output at index: %d\n",index_t);
  mat_file = fopen(pqke->output_filename,"r+b");
  //Write stuff from y_vector:
  L = y[pqke->index_L]*_L_SCALE_;
  mat_write_fast(&L,pqke->L_handle,8,1);
  mat_write_fast(Pa_plus_ptr,pqke->Pa_plus_handle,8,vres);
  mat_write_fast(Pa_minus_ptr,pqke->Pa_minus_handle,8,vres);
  mat_write_fast(Ps_plus_ptr,pqke->Ps_plus_handle,8,vres);
  mat_write_fast(Ps_minus_ptr,pqke->Ps_minus_handle,8,vres);
  mat_write_fast(Px_plus_ptr,pqke->Px_plus_handle,8,vres);
  mat_write_fast(Px_minus_ptr,pqke->Px_minus_handle,8,vres);
  mat_write_fast(Py_plus_ptr,pqke->Py_plus_handle,8,vres);
  mat_write_fast(Py_minus_ptr,pqke->Py_minus_handle,8,vres);
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
  free(Pa_plus_ptr);
  return _SUCCESS_;
};

int qke_store_output_fixed_grid_standard_approx(double T,
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
  double *out_vec, *Pa_plus_ptr, *Pa_minus_ptr, *Ps_plus_ptr, *Ps_minus_ptr;
  double *Px_plus_ptr, *Px_minus_ptr, *Py_plus_ptr, *Py_minus_ptr;
  double *x_grid = pqke->x_grid;
  double Vx,V0,V1,w,s;
  double P0, P0_bar, Px, Px_bar, Py, Py_bar, Pz, Pz_bar;
  double w_trapz,I_Paminus;

  Pa_plus_ptr = malloc(4*sizeof(double)*vres);
  Pa_minus_ptr = Pa_plus_ptr+vres;
  Ps_plus_ptr = Pa_plus_ptr+2*vres;
  Ps_minus_ptr = Pa_plus_ptr+3*vres;

  I_Paminus = 0.0;
  for (i=0; i<vres; i++){
    
    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    Pz = y[pqke->index_Pz+i];
    Pz_bar = y[pqke->index_Pz_bar+i];


    Pa_plus_ptr[i] = P0+P0_bar+Pz+Pz_bar;
    Pa_minus_ptr[i] = P0-P0_bar+Pz-Pz_bar;;
    Ps_plus_ptr[i] = P0+P0_bar-(Pz+Pz_bar);
    Ps_minus_ptr[i] = P0-P0_bar-(Pz-Pz_bar);

    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    f0 = 1.0/(exp(x)+1.0);
      
    I_Paminus += w_trapz*x*x*f0*Pa_minus_ptr[i];
  }

  L = pqke->L_fudge*I_Paminus/(8.0*_ZETA3_);    
  //printf("Storing output at index: %d\n",index_t);
  mat_file = fopen(pqke->output_filename,"r+b");
  //Write stuff from y_vector:
//  L = y[pqke->index_L]*_L_SCALE_;
  mat_write_fast(&L,pqke->L_handle,8,1);
  mat_write_fast(Pa_plus_ptr,pqke->Pa_plus_handle,8,vres);
  mat_write_fast(Pa_minus_ptr,pqke->Pa_minus_handle,8,vres);
  mat_write_fast(Ps_plus_ptr,pqke->Ps_plus_handle,8,vres);
  mat_write_fast(Ps_minus_ptr,pqke->Ps_minus_handle,8,vres);
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
  free(Pa_plus_ptr);
  return _SUCCESS_;
};

int qke_store_output_fixed_grid_standard_approx2(double T,
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
  double *out_vec, *Pa_plus_ptr, *Pa_minus_ptr, *Ps_plus_ptr, *Ps_minus_ptr;
  double *Px_plus_ptr, *Px_minus_ptr, *Py_plus_ptr, *Py_minus_ptr;
  double *x_grid = pqke->x_grid;
  double Vx,V0,V1,w,s;
  double P0, P0_bar, Px, Px_bar, Py, Py_bar, Pz, Pz_bar;
  double w_trapz,I_Paminus;
  double Vz,VL,Vz_bar,gamma,gamma_bar,delta_Pz_ini,delta_Pz_bar_ini,Pz_pure;
  double P0_bar_ini, delta_Pz, delta_Pz_bar,P0_ini,Pz_bar_pure;

  Pa_plus_ptr = malloc(4*sizeof(double)*vres);
  Pa_minus_ptr = Pa_plus_ptr+vres;
  Ps_plus_ptr = Pa_plus_ptr+2*vres;
  Ps_minus_ptr = Pa_plus_ptr+3*vres;

  I_Paminus = 0.0;
  for (i=0; i<vres; i++){
    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Vz = V0+V1+VL;
    Vz_bar = V0+V1-VL;
    gamma = Vz/sqrt(Vx*Vx+Vz*Vz);
    gamma_bar = Vz_bar/sqrt(Vx*Vx+Vz_bar*Vz_bar);

    delta_Pz_ini = pqke->delta_Pz_ini[i];
    delta_Pz_bar_ini = pqke->delta_Pz_bar_ini[i];
    P0_ini = pqke->P0_ini[i];
    P0_bar_ini = pqke->P0_bar_ini[i];

    delta_Pz = delta_Pz_ini+P0-P0_ini;
    delta_Pz_bar = delta_Pz_bar_ini+P0_bar-P0_bar_ini;

    Pz_pure = gamma*pqke->sigma0Pz[i];
    Pz_bar_pure = gamma_bar*pqke->sigma0Pz_bar[i];
      
    Pz = Pz_pure+delta_Pz;
    Pz_bar = Pz_bar_pure+delta_Pz_bar;


    Pa_plus_ptr[i] = P0+P0_bar+Pz+Pz_bar;
    Pa_minus_ptr[i] = P0-P0_bar+Pz-Pz_bar;;
    Ps_plus_ptr[i] = P0+P0_bar-(Pz+Pz_bar);
    Ps_minus_ptr[i] = P0-P0_bar-(Pz-Pz_bar);
  }

  L = pqke->L0+pqke->Lz;    
  //printf("Storing output at index: %d\n",index_t);
  mat_file = fopen(pqke->output_filename,"r+b");
  //Write stuff from y_vector:
//  L = y[pqke->index_L]*_L_SCALE_;
  mat_write_fast(&L,pqke->L_handle,8,1);
  mat_write_fast(Pa_plus_ptr,pqke->Pa_plus_handle,8,vres);
  mat_write_fast(Pa_minus_ptr,pqke->Pa_minus_handle,8,vres);
  mat_write_fast(Ps_plus_ptr,pqke->Ps_plus_handle,8,vres);
  mat_write_fast(Ps_minus_ptr,pqke->Ps_minus_handle,8,vres);
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
  free(Pa_plus_ptr);
  return _SUCCESS_;
};


int qke_stop_at_L(double t,
		  double *y,
		  double *dy,
		  void *param,
		  ErrorMsg error_message){
  qke_param *pqke=param;
  if (fabs(y[pqke->index_L]*_L_SCALE_) <= pqke->L_final){
    return _TRUE_;
  }
  else{
    return _FALSE_;
  }
};

int qke_print_variables(double T,
			double *y,
			double *dy,
			void *param,
			ErrorMsg error_message){
  
  qke_param *pqke=param;
  int i,idx=53;
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
	  V0+V1,
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

int qke_print_variables_fixed_grid_low_temp(double T,
					    double *y,
					    double *dy,
					    void *param,
					    ErrorMsg error_message){
  
  qke_param *pqke=param;
  int i,idx=53;
  double x,Vx,V0,V1,VL,w,s;
  double D,Gamma;
  double Pa_plus, Pa_minus, Ps_plus, Ps_minus, Px_plus, Px_minus;
  double Py_plus, Py_minus;
  double Q0_plus,Q0_minus,Q1_plus,Q1_minus,Q2_plus,Q2_minus,Q3_plus,Q3_minus;

  x = pqke->x_grid[idx];
    
  Vx = pqke->Vx/x;
  V0 = pqke->V0/x;
  V1 = pqke->V1*x;
  w = (V0+V1)/Vx;
  s = sqrt(1.0+w*w);
  VL = pqke->VL;
  Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
  D = 0.5*Gamma;
  
  
  Q0_plus = y[pqke->index_Q0_plus+idx];
  Q1_plus = y[pqke->index_Q1_plus+idx];
  Q2_plus = y[pqke->index_Q2_plus+idx];
  Q3_plus = y[pqke->index_Q3_plus+idx];
  Q0_minus = y[pqke->index_Q0_minus+idx];
  Q1_minus = y[pqke->index_Q1_minus+idx];
  Q2_minus = y[pqke->index_Q2_minus+idx];   
  Q3_minus = y[pqke->index_Q3_minus+idx];
  
  Pa_plus= Q0_plus+(Q1_plus+Q3_plus)/s/s;
  Pa_minus = Q0_minus+(Q1_minus+Q3_minus)/s/s;
  Ps_plus = Q0_plus-(Q1_plus+Q3_plus)/s/s;
  Ps_minus = Q0_minus-(Q1_minus+Q3_minus)/s/s;
  Px_plus = -w/s/s*Q1_plus+1.0/(w*s*s)*Q3_plus;
  Px_minus = -w/s/s*Q1_minus+1.0/(w*s*s)*Q3_minus;
  Py_plus = 1.0/s*Q2_plus;
  Py_minus = 1.0/s*Q2_minus;
  
  

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
	  V0+V1,
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

int qke_print_variables_fixed_grid_standard(double T,
					    double *y,
					    double *dy,
					    void *param,
					    ErrorMsg error_message){
  
  qke_param *pqke=param;
  int i,idx=53;
  double x,Vx,V0,V1,VL,w,s;
  double D,Gamma;
  double Pa_plus, Pa_minus, Ps_plus, Ps_minus, Px_plus, Px_minus;
  double Py_plus, Py_minus;
  double P0, P0_bar, Px, Px_bar, Py, Py_bar, Pz, Pz_bar;
  double dlogVzdlogT,dlogVz_bardlogT;
  double L,dLdT;
  
  L = y[pqke->index_L]*_L_SCALE_;
  dLdT = dy[pqke->index_L]*_L_SCALE_;;
  x = pqke->x_grid[idx];
    
  Vx = pqke->Vx/x;
  V0 = pqke->V0/x;
  V1 = pqke->V1*x;
  w = (V0+V1)/Vx;
  s = sqrt(1.0+w*w);
  VL = pqke->VL;
  Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
  D = 0.5*Gamma;

  dlogVzdlogT = (-V0+5.0*V1+VL*(3.0+T/L*dLdT))/(V0+V1+VL);
  dlogVz_bardlogT = (-V0+5.0*V1-VL*(3.0+T/L*dLdT))/(V0+V1-VL);

  P0 = y[pqke->index_P0+idx];
  P0_bar = y[pqke->index_P0_bar+idx];
  Px = y[pqke->index_Px+idx];
  Px_bar = y[pqke->index_Px_bar+idx];
  Py = y[pqke->index_Py+idx];
  Py_bar = y[pqke->index_Py_bar+idx];
  Pz = y[pqke->index_Pz+idx];
  Pz_bar = y[pqke->index_Pz_bar+idx];
    
  
  Pa_plus = P0+P0_bar+Pz+Pz_bar;
  Pa_minus = P0-P0_bar+Pz-Pz_bar;
  Ps_plus = P0+P0_bar-(Pz+Pz_bar);
  Ps_minus = P0-P0_bar-(Pz-Pz_bar);
  Px_plus = Px+Px_bar;
  Px_minus = Px-Px_bar;
  Py_plus = Py+Py_bar;
  Py_minus = Py-Py_bar;
  
  
  
  fprintf(stderr,
	  "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	  1e3*T,
	  L,
	  P0,
	  P0_bar,
	  Pz,
	  Pz_bar,
	  Px,
	  Px_bar,
	  Py,
	  Py_bar,
	  dlogVzdlogT,
	  dlogVz_bardlogT,
	  -V0/T,
	  5.0*V1/T,
	  VL*(3.0+T/L*dLdT)/T,
	  VL*3.0/T,
	  VL/L*dLdT,
	  D*Py_minus);
 
  /*
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
	  V0+V1,
	  VL,
	  -(V0+V1)*Px_plus-VL*Px_minus, //15
	  0.5*Vx*(Pa_plus-Ps_plus),
	  D*Py_plus,
	  -(V0+V1)*Px_minus-VL*Px_plus,
	  0.5*Vx*(Pa_minus-Ps_minus),
	  D*Py_minus);

  */
//qke_param *pqke=param;
 
  return _SUCCESS_;
};

int qke_print_variables_fixed_grid_standard_approx(double T,
						   double *y,
						   double *dy,
						   void *param,
						   ErrorMsg error_message){
  
  qke_param *pqke=param;
  int i,idx=53;
  double x,Vx,V0,V1,VL,w,s;
  double D,Gamma;
  double Pa_plus, Pa_minus, Ps_plus, Ps_minus, Px_plus, Px_minus;
  double Py_plus, Py_minus;
  double P0, P0_bar, Px, Px_bar, Py, Py_bar, Pz, Pz_bar;
  double w_trapz,I_Pyminus,L,f0;
  double dlogVzdlogT,dlogVz_bardlogT,dLdT;

  dLdT = pqke->dLdT;
  I_Pyminus = 0.0;
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(pqke->x_grid[i+1]-pqke->x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(pqke->x_grid[i]-pqke->x_grid[i-1]);
    else
      w_trapz = 0.5*(pqke->x_grid[i+1]-pqke->x_grid[i-1]);
    x = pqke->x_grid[i];
    f0 = 1.0/(1.0+exp(x));

    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    Pz = y[pqke->index_Pz+i];
    Pz_bar = y[pqke->index_Pz_bar+i];

    Pa_minus = P0-P0_bar+Pz-Pz_bar;


    I_Pyminus += w_trapz*x*x*f0*Pa_minus;
  }
  L = pqke->L_fudge*I_Pyminus/(8.0*_ZETA3_);

  x = pqke->x_grid[idx];
    
  Vx = pqke->Vx/x;
  V0 = pqke->V0/x;
  V1 = pqke->V1*x;
  w = (V0+V1)/Vx;
  s = sqrt(1.0+w*w);
  VL = pqke->VL;
  Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
  D = 0.5*Gamma;

  dlogVzdlogT = (-V0+5.0*V1+VL*(3.0+T/L*dLdT))/(V0+V1+VL);
  dlogVz_bardlogT = (-V0+5.0*V1-VL*(3.0+T/L*dLdT))/(V0+V1-VL);

  
  P0 = y[pqke->index_P0+idx];
  P0_bar = y[pqke->index_P0_bar+idx];
  Pz = y[pqke->index_Pz+idx];
  Pz_bar = y[pqke->index_Pz_bar+idx];
    
  
  Pa_plus = P0+P0_bar+Pz+Pz_bar;
  Pa_minus = P0-P0_bar+Pz-Pz_bar;
  Ps_plus = P0+P0_bar-(Pz+Pz_bar);
  Ps_minus = P0-P0_bar-(Pz-Pz_bar);
  Py_minus = 0.0;
  Py_plus=0.0;
  Px_minus = 0.0;
  Px_plus = 0.0;

  fprintf(stderr,
	  "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	  1e3*T,
	  L,
	  P0,
	  P0_bar,
	  Pz,
	  Pz_bar,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  dlogVzdlogT,
	  dlogVz_bardlogT,
	  -V0/T,
	  5.0*V1/T,
	  VL*(3.0+T/L*dLdT)/T,
	  VL*3.0/T,
	  VL/L*dLdT,
	  D*Py_minus);
    
  /*
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
	  V0+V1,
	  VL,
	  -(V0+V1)*Px_plus-VL*Px_minus, //15
	  0.5*Vx*(Pa_plus-Ps_plus),
	  D*Py_plus,
	  -(V0+V1)*Px_minus-VL*Px_plus,
	  0.5*Vx*(Pa_minus-Ps_minus),
	  D*Py_minus);

  */
//qke_param *pqke=param;
 
  return _SUCCESS_;
};

int qke_print_variables_fixed_grid_standard_approx2(double T,
						    double *y,
						    double *dy,
						    void *param,
						    ErrorMsg error_message){
  
  qke_param *pqke=param;
  int i,idx=53;
  double x,Vx,V0,V1,VL,w,s;
  double D,Gamma;
  double Pa_plus, Pa_minus, Ps_plus, Ps_minus, Px_plus, Px_minus;
  double Py_plus, Py_minus;
  double P0, P0_bar, Px, Px_bar, Py, Py_bar, Pz, Pz_bar;
  double w_trapz,I_Pyminus,L,f0,Vz,Vz_bar,delta_Pz_ini,delta_Pz_bar_ini;
  double dlogVzdlogT,dlogVz_bardlogT,dLdT;
  double P0_ini,P0_bar_ini,gamma,gamma_bar,delta_Pz,delta_Pz_bar,Pz_pure,Pz_bar_pure;

  L = pqke->L0+pqke->Lz;

  x = pqke->x_grid[idx];
    
  Vx = pqke->Vx/x;
  V0 = pqke->V0/x;
  V1 = pqke->V1*x;
  w = (V0+V1)/Vx;
  s = sqrt(1.0+w*w);
  VL = pqke->VL;
  Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
  D = 0.5*Gamma;

  gamma = Vz/sqrt(Vx*Vx+Vz*Vz);
  gamma_bar = Vz_bar/sqrt(Vx*Vx+Vz_bar*Vz_bar);

  delta_Pz_ini = pqke->delta_Pz_ini[idx];
  delta_Pz_bar_ini = pqke->delta_Pz_bar_ini[idx];
  
  P0 = y[pqke->index_P0+idx];
  P0_bar = y[pqke->index_P0_bar+idx];


  P0_ini = pqke->P0_ini[idx];
  P0_bar_ini = pqke->P0_bar_ini[idx];

  delta_Pz = delta_Pz_ini+P0-P0_ini;
  delta_Pz_bar = delta_Pz_bar_ini+P0_bar-P0_bar_ini;

  Pz_pure = gamma*pqke->sigma0Pz[idx];
  Pz_bar_pure = gamma_bar*pqke->sigma0Pz_bar[idx];

  Pz = Pz_pure+delta_Pz;
  Pz_bar = Pz_bar_pure+delta_Pz_bar;


  
  Pa_plus = P0+P0_bar+Pz+Pz_bar;
  Pa_minus = P0-P0_bar+Pz-Pz_bar;
  Ps_plus = P0+P0_bar-(Pz+Pz_bar);
  Ps_minus = P0-P0_bar-(Pz-Pz_bar);
  Py_minus = 0.0;
  Py_plus=0.0;
  Px_minus = 0.0;
  Px_plus = 0.0;

  fprintf(stderr,
	  "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	  1e3*T,
	  L,
	  P0,
	  P0_bar,
	  Pz,
	  Pz_bar,
	  0.0,
	  0.0,
	  0.0,
	  0.0,
	  dlogVzdlogT,
	  dlogVz_bardlogT,
	  -V0/T,
	  5.0*V1/T,
	  VL*(3.0+T/L*dLdT)/T,
	  VL*3.0/T,
	  VL/L*dLdT,
	  D*Py_minus);
    
  /*
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
	  V0+V1,
	  VL,
	  -(V0+V1)*Px_plus-VL*Px_minus, //15
	  0.5*Vx*(Pa_plus-Ps_plus),
	  D*Py_plus,
	  -(V0+V1)*Px_minus-VL*Px_plus,
	  0.5*Vx*(Pa_minus-Ps_minus),
	  D*Py_minus);

  */
//qke_param *pqke=param;
 
  return _SUCCESS_;
};


int qke_derivs(double T, 
	       double *y, 
	       double *dy, 
	       void *param,
	       ErrorMsg error_message){
  qke_param *pqke=param;
  double L;
  int i,j,niter;
  double *y_0 = pqke->y_0;
  double *maxstep=pqke->maxstep;
  double *duidT=pqke->duidT;
  double *dvidT=pqke->dvidT;
  double *duidx=pqke->duidx;
  double *ui=pqke->ui;
  double *vi=pqke->vi;
  double *u_grid=pqke->u_grid;
  double *v_grid=pqke->v_grid;
  double *x_grid=pqke->x_grid;
  double *dvdu_grid=pqke->dvdu_grid;
  double *dudT_grid=pqke->dudT_grid;
  double tol_newton = 1e-12;
  double vN,u1,v1,alpha,wi;
  double gentr,H;
  double lu_sgn;
  double daidT,dbdT;
  double Vx, VL;
  double x, w_trapz;
  double Gamma, D, V0, V1, Pa_plus, Pa_minus, Ps_plus, Ps_minus;
  double Px_plus, Px_minus, Py_plus, Py_minus, f0,  mu_div_T;
  double PsPs,I_rho_ss, feq_plus, feq_minus, I_VxPy_minus, I_f0Pa_plus;
  double dudTdvdu, delta_v;
  double rs;
  double gamma_j, beta_j;
  int idx, stencil_method;

  alpha = pqke->alpha;
  L = y[pqke->index_L]*_L_SCALE_;
  get_resonances_xi(T,L,pqke);
  for (i=0; i<pqke->Nres; i++){
    u_of_x(pqke->xi[i],pqke->ui+i,duidx+i,pqke);
  }
  //Set partial derivatives of ui with respect to T:
  for (i=0; i<pqke->Nres; i++) 
    pqke->duidT[i]=pqke->duidx[i]*pqke->dxidT[i];


  /**
     Do we evolve the parametrisation parameters b and vi or
     do we use Newtons method?
  */
  if (pqke->evolve_vi == _FALSE_){
    //Establish guess and set maximum steps for Newton method:
    if (T==pqke->T_initial){
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
    
    pqke->b = y_0[0];
    for (i=0; i<pqke->Nres; i++)
      vi[i] = y_0[i+1];
  }
  else{
    //Enter vi and b from the solution vector:
    //pqke->b = y[pqke->index_b];
    for (i=0; i<pqke->Nres; i++)
      vi[i] = y[pqke->index_vi+i];
    pqke->b = ui[0]/pow(vi[0],3)-alpha/pow(vi[0],2);
  }

  //Set remaining values of the parametrisation:
  for (i=0; i<pqke->Nres; i++)
    pqke->a[i] = ui[i]-alpha*vi[i];

  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(pqke->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  /** Get partial derivatives of vi with respect to T by solving a 
      linear system A*dvidT = B(duidT):
  */
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
      u_grid[j] = alpha*v_grid[j]+pqke->a[i]+pqke->b*pow(v_grid[j]-vi[i],3);
      x_of_u(u_grid[j],&(x_grid[j]),pqke);
      dvdu_grid[j] = 1.0/(alpha+3.0*pqke->b*pow(v_grid[j]-vi[i],2));
      dudT_grid[j] = daidT+dbdT*pow(v_grid[j]-vi[i],3)-
	3.0*pqke->b*pow(v_grid[j]-vi[i],2);
    }
  }
  
  /** Done calculating coordinate transformations, now we can 
      calculate the RHS of the ODE: */

  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);
  
  /** Integrated quantities needed. We integrate in x space: */
  I_VxPy_minus = 0.0;
  I_f0Pa_plus = 0.0;
  I_rho_ss = 0.0;
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(pqke->x_grid[i+1]-pqke->x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(pqke->x_grid[i]-pqke->x_grid[i-1]);
    else
      w_trapz = 0.5*(pqke->x_grid[i+1]-pqke->x_grid[i-1]);
    x = pqke->x_grid[i];
    f0 = 1.0/(1.0+exp(x));
    Vx = pqke->Vx/x;
    Py_minus = y[pqke->index_Py_minus+i];
    Pa_plus = y[pqke->index_Pa_plus+i];
    PsPs = y[pqke->index_Ps_plus+i]+y[pqke->index_Ps_minus+i];
    
    I_VxPy_minus += w_trapz*(x*x*Vx*Py_minus*f0);
    I_f0Pa_plus += w_trapz*(x*x*f0*Pa_plus);
    I_rho_ss += w_trapz*(x*x*f0*PsPs);
  }
  //Set V1:
  if (pqke->is_electron == _TRUE_){
    pqke->V1 = -14.0*sqrt(2.0)*_PI_*_PI_/45.0*_G_F_/_M_W_/_M_W_*pow(T,5)*
      (1.0+1.0/(12.0*_ZETA3_)*(1.0-_SIN2_THETA_W_)*I_f0Pa_plus);
  }
  else{
    pqke->V1 = -7.0*_PI_*_PI_/(135.0*sqrt(2.0)*_ZETA3_)*
      _G_F_/_M_Z_/_M_Z_*pow(T,5)*I_f0Pa_plus;
  }

  /** Calculate RHS: */
  dy[pqke->index_L] = 0.0;
  //dy[pqke->index_L] = -1.0/(8.0*H*T*_ZETA3_)*I_VxPy_minus/_L_SCALE_;
  /** All quantities defined on the grid: */

  // Set perhaps flow of grid:
  if (pqke->evolve_vi == _TRUE_){
    //dy[pqke->index_b] = dbdT;
    for (i=0; i<pqke->Nres; i++){
      dy[pqke->index_vi+i] = dvidT[i];
    }
  }

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
  
    //Solving mu from L, using the chebyshev cubic root:
    mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
    //Regulator for sterile population:
    rs = pqke->rs;
    //Distributions:
    feq_plus = 1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T));
    feq_minus = 1.0/(1.0+exp(x-mu_div_T))-1.0/(1.0+exp(x+mu_div_T));
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
    dy[idx] = 1.0/(H*T)*
      (Vx*Py_plus-rs*Gamma*(1.0/(6.0*_ZETA3_)*I_rho_ss*feq_plus-0.5*f0*Ps_plus))+     dudTdvdu*drhodv(y, delta_v, idx, stencil_method);
    
    idx = pqke->index_Ps_minus+i;
    dy[idx] = 1.0/(H*T)*
      (Vx*Py_minus-rs*Gamma*f0*(1.0-0.5*(Pa_minus+Ps_minus)))+
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
  double Gamma, D, V0, V1, Pa_plus, Pa_minus, Ps_plus, Ps_minus;
  double Px_plus, Px_minus, Py_plus, Py_minus, f0,  mu_div_T;
  double PsPs,I_rho_ss, feq_plus, feq_minus, I_VxPy_minus, I_f0Pa_plus;
  double rs;
  int idx;

  double V_scale=1.0,Vz,erfres,T_fiducial=1e-4,k3, mu_fiducial=0.0023;
  double L_from_int, I_Paminus, dLdT;
  

  L = y[pqke->index_L]*_L_SCALE_;


 //L = 1e-2*(1.0-1.0/(exp((T*1e3-2.486022964891155)/0.09)+1.0));
  /**
  if (T<0.0025){
    if (pqke->L1<-1e99){
      pqke->L1 = L*(1.0+exp((-T+mu_fiducial)/T_fiducial));
      printf("L1: %.16e\n",pqke->L1);
      //return _FAILURE_;
    }
    L = pqke->L1/(exp((-T+mu_fiducial)/T_fiducial)+1.0);
    printf("L = %g\n",L);
  }
  */

  get_resonances_xi(T,L,pqke);

  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(pqke->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);
  
  /** Integrated quantities needed. We integrate in x space: */
  I_VxPy_minus = 0.0;
  I_f0Pa_plus = 0.0;
  I_rho_ss = 0.0;
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    f0 = 1.0/(1.0+exp(x));
    Vx = pqke->Vx/x;
    Py_minus = y[pqke->index_Py_minus+i];
    Pa_plus = y[pqke->index_Pa_plus+i];
    PsPs = y[pqke->index_Ps_plus+i]+y[pqke->index_Ps_minus+i];
    
    I_VxPy_minus += w_trapz*(x*x*Vx*Py_minus*f0);
    I_f0Pa_plus += w_trapz*(x*x*f0*Pa_plus);
    I_rho_ss += w_trapz*(x*x*f0*PsPs);
  }
  //Set V1:
  if (pqke->is_electron == _TRUE_){
    pqke->V1 = -14.0*sqrt(2.0)*_PI_*_PI_/45.0*_G_F_/_M_W_/_M_W_*pow(T,5)*
      (1.0+1.0/(12.0*_ZETA3_)*(1.0-_SIN2_THETA_W_)*I_f0Pa_plus);
  }
  else{
    pqke->V1 = -7.0*_PI_*_PI_/(135.0*sqrt(2.0)*_ZETA3_)*
      _G_F_/_M_Z_/_M_Z_*pow(T,5)*I_f0Pa_plus;
  }

  /** Calculate RHS: */ 
  dLdT = -1.0/(8.0*H*T*_ZETA3_)*I_VxPy_minus;
  dy[pqke->index_L] = dLdT/_L_SCALE_;
  
  if ((pqke->L_decay_trigger == _FALSE_)&&(dLdT/L>1e3)){
    pqke->L_decay_trigger = _TRUE_;
    printf("L decay trigger active at T=%g MeV.\n",T*1e3);
    printf("L=%g, dLdT=%g 1/MeV\n",L,dLdT*1e3);
  }

  //  if (fabs(L)<0.9*fabs(pqke->L_initial)){
  if (pqke->L_decay_trigger == _TRUE_){
    if (pqke->dLdT_approx == _FALSE_){
      pqke->mu = dy[pqke->index_L]/L;
      pqke->T1 = T;
      pqke->dLdT_approx = _TRUE_;
    }
    dy[pqke->index_L] = L*(pqke->mu+6.9e10*pow(pqke->T1-T,2)/_L_SCALE_);
  }
    

  /** All quantities defined on the grid: */

  
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Vz = V0+V1+VL;
    V_scale = max(fabs(V0),fabs(V1));
    V_scale = max(V_scale,fabs(VL));
    V_scale *= 0.5;
    erfres = erf(pow(fabs(Vz/V_scale),3.0));

    //VL = 0.5*VL*(erfres+1.0)-0.5*(V0+V1)*(1.0-erfres);

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
  
    //Solving mu from L, using the chebyshev cubic root:
    mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
    //Regulator for sterile population:
    rs = pqke->rs;
    //Distributions:
    feq_plus = 1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T));
    feq_minus = 1.0/(1.0+exp(x-mu_div_T))-1.0/(1.0+exp(x+mu_div_T));
    f0 = 1.0/(1.0+exp(x));
    
    idx = pqke->index_Pa_plus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_plus+Gamma*(2.0*feq_plus/f0-Pa_plus));

    idx = pqke->index_Pa_minus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_minus+Gamma*(2.0*feq_minus/f0-Pa_minus));

    idx = pqke->index_Ps_plus+i;
    dy[idx] = 1.0/(H*T)*
      (Vx*Py_plus-rs*Gamma*(1.0/(6.0*_ZETA3_)*I_rho_ss*feq_plus-0.5*f0*Ps_plus));
    
    idx = pqke->index_Ps_minus+i;
    dy[idx] = 1.0/(H*T)*
      (Vx*Py_minus-rs*Gamma*f0*(1.0-0.5*(Pa_minus+Ps_minus)));
    
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

int qke_derivs_fixed_grid_low_temp(double T, 
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
  double Gamma, D, V0, V1, Pa_plus, Pa_minus, Ps_plus, Ps_minus;
  double Px_plus, Px_minus, Py_plus, Py_minus, f0,  mu_div_T;
  double Q0_plus,Q0_minus,Q1_plus,Q1_minus,Q2_plus,Q2_minus,Q3_plus,Q3_minus;
  double I_df0Pa_plus, feq_plus, feq_minus, I_VxPy_minus, I_f0Pa_plus;
  double dPa_plusdt,dP0_plusdT,dP0_minusdT,dLdT,dn_plusdT,n_plus;
  double rs;
  int idx;
  double A_plus[3],A_minus[3],B_plus[3],B_minus[3],C_plus[3],C_minus[3];
  double c_A, c_B, c_C;
  double w,s,dwdT,x2;


  L = y[pqke->index_L]*_L_SCALE_;
  n_plus = y[pqke->index_n_plus];

  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(pqke->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);
  //Set V1:
  if (pqke->is_electron == _TRUE_){
    pqke->V1 = -14.0*sqrt(2.0)*_PI_*_PI_/45.0*_G_F_/_M_W_/_M_W_*pow(T,5)*
      (1.0+0.25*(1.0-_SIN2_THETA_W_)*n_plus);
  }
  else{
    pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
      _G_F_/_M_Z_/_M_Z_*pow(T,5)*n_plus;
  }
  
  /** Integrated quantities needed. We integrate in x space: */
  I_VxPy_minus = 0.0;
  I_df0Pa_plus = 0.0;
  
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    x2 = x*x;

    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    w = (V0+V1)/Vx;
    s = sqrt(1.0+w*w);
    //Solving mu from L, using the chebyshev cubic root:
    mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
    //Distributions:
    feq_plus = 1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T));
    feq_minus = 1.0/(1.0+exp(x-mu_div_T))-1.0/(1.0+exp(x+mu_div_T));
    f0 = 1.0/(1.0+exp(x));

    Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
    D = 0.5*Gamma;

    Q0_plus = y[pqke->index_Q0_plus+i];
    Q1_plus = y[pqke->index_Q1_plus+i];
    Q2_plus = y[pqke->index_Q2_plus+i];
    Q3_plus = y[pqke->index_Q3_plus+i];
    Q2_minus = y[pqke->index_Q2_minus+i];

    Pa_plus = Q0_plus+(Q1_plus+Q3_plus)/s/s;
    dPa_plusdt = Vx/s*Q2_plus+Gamma*(2.0*feq_plus/f0-Pa_plus);
    
    I_df0Pa_plus += w_trapz*(x2*f0*dPa_plusdt);
    I_VxPy_minus += w_trapz*(x2*f0*Vx/s*Q2_minus);
  }
  dLdT = -1.0/(8.0*H*T*_ZETA3_)*I_VxPy_minus;
  dn_plusdT = -1.0/(3.0*H*T*_ZETA3_)*I_df0Pa_plus;

  /** Calculate RHS: */
  dy[pqke->index_L] = dLdT/_L_SCALE_;
  //dy[pqke->index_L] = 0.0;
  dy[pqke->index_n_plus] = dn_plusdT;

  /** All quantities defined on the grid: */

  
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];
    x2 = x*x;
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    w = (V0+V1)/Vx;
    s = sqrt(1.0+w*w);
   
    //Calculate dw/dT:
    if (pqke->is_electron == _TRUE_){
      dwdT = V1/(T*Vx)*(6.0+T*
			(1.0+0.25*(1.0-_SIN2_THETA_W_)*dn_plusdT)/
			(1.0+0.25*(1.0-_SIN2_THETA_W_)*n_plus));
    }
    else{
      dwdT = V1/(T*Vx)*(6.0+T*dn_plusdT/n_plus);
    }

    Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
    D = 0.5*Gamma;
    
    Q0_plus = y[pqke->index_Q0_plus+i];
    Q0_minus = y[pqke->index_Q0_minus+i];
    Q1_plus = y[pqke->index_Q1_plus+i];
    Q1_minus = y[pqke->index_Q1_minus+i];
    Q2_plus = y[pqke->index_Q2_plus+i];
    Q2_minus = y[pqke->index_Q2_minus+i];
    Q3_plus = y[pqke->index_Q3_plus+i];
    Q3_minus = y[pqke->index_Q3_minus+i];

    dP0_plusdT = -Gamma/(H*T*f0)*
      (feq_plus-0.5*f0*(Q0_plus+(Q1_plus+Q3_plus)/s/s));
    dP0_minusdT = -Gamma/(H*T*f0)*
      (feq_minus-0.5*f0*(Q0_minus+(Q1_minus+Q3_minus)/s/s));
    
    //Calculate matrix operations:
    matrix_x_on_Q(A_plus, w, D/Vx, y,
		  pqke->index_Q1_plus+i,
		  pqke->index_Q2_plus-pqke->index_Q1_plus,
		  'A');
    matrix_x_on_Q(A_minus, w, D/Vx, y,
		  pqke->index_Q1_minus+i,
		  pqke->index_Q2_minus-pqke->index_Q1_minus,
		  'A');
    matrix_x_on_Q(B_plus, w, D/Vx, y,
		  pqke->index_Q1_plus+i,
		  pqke->index_Q2_plus-pqke->index_Q1_plus,
		  'B');
    matrix_x_on_Q(B_minus, w, D/Vx, y,
		  pqke->index_Q1_minus+i,
		  pqke->index_Q2_minus-pqke->index_Q1_minus,
		  'B');
    matrix_x_on_Q(C_plus, w, D/Vx, y,
		  pqke->index_Q1_minus+i,
		  pqke->index_Q2_minus-pqke->index_Q1_minus,
		  'C');
    matrix_x_on_Q(C_minus, w, D/Vx, y,
		  pqke->index_Q1_plus+i,
		  pqke->index_Q2_plus-pqke->index_Q1_plus,
		  'C');
  
    //Calculate coefficients for matrix equation:
    c_A = w*dwdT/s/s;
    c_B = -Vx/(H*T);
    c_C = -VL*w/(H*T*s);

    //Solving mu from L, using the chebyshev cubic root:
    mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
    //Distributions:
    feq_plus = 1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T));
    feq_minus = 1.0/(1.0+exp(x-mu_div_T))-1.0/(1.0+exp(x+mu_div_T));
    f0 = 1.0/(1.0+exp(x));
    
    idx = pqke->index_Q0_plus+i;
    dy[idx] = dP0_plusdT;

    idx = pqke->index_Q0_minus+i;
    dy[idx] = dP0_minusdT;

    idx = pqke->index_Q1_plus+i;
    dy[idx] = c_A*A_plus[0]+c_B*B_plus[0]+c_C*C_plus[0]+dP0_plusdT;
    
    idx = pqke->index_Q1_minus+i;
    dy[idx] = c_A*A_minus[0]+c_B*B_minus[0]+c_C*C_minus[0]+dP0_minusdT;
    
    idx = pqke->index_Q2_plus+i;
    dy[idx] = c_A*A_plus[1]+c_B*B_plus[1]+c_C*C_plus[1]+dP0_plusdT;

    idx = pqke->index_Q2_minus+i;
    dy[idx] = c_A*A_minus[1]+c_B*B_minus[1]+c_C*C_minus[1]+dP0_minusdT;

    idx = pqke->index_Q3_plus+i;
    dy[idx] = c_A*A_plus[2]+c_B*B_plus[2]+c_C*C_plus[2]+dP0_plusdT/w/w;
    
    idx = pqke->index_Q3_minus+i;
    dy[idx] = c_A*A_minus[2]+c_B*B_minus[2]+c_C*C_minus[2]+dP0_minusdT/w/w;
  }
  return _SUCCESS_;
}

void matrix_x_on_Q(double *output, 
		   double w,
		   double D,
		   double *input, 
		   int start_idx, 
		   int skip_idx,
		   char x){
  double w2=w*w;
  double s2=1+w2;
  double s=sqrt(s2);
  if (x=='A'){
    //Matrix A
    output[0] = input[start_idx]-input[start_idx+2*skip_idx]/w2;
    output[1] = input[start_idx+skip_idx];
    output[2] = input[start_idx]+input[start_idx+2*skip_idx]*(1.0+s2/w2);
    return;
  }
  else if (x=='B'){
    //Matrix B
    output[0] = -D*w2/s2*input[start_idx]+
      s*input[start_idx+skip_idx]+D/s2*input[start_idx+2*skip_idx];
    output[1] = -s*input[start_idx]-D*input[start_idx+skip_idx];
    output[2] = D*w2/s2*input[start_idx]-D/s2*input[start_idx+2*skip_idx];
    return;
  }
  else if (x=='C'){
    //Matrix C
    output[0] = input[start_idx+skip_idx];
    output[1] = -input[start_idx]+input[start_idx+2*skip_idx]/w2;
    output[2] = -output[0];
    return;
  }
  return;
}

int qke_derivs_fixed_grid_standard(double T, 
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
  double x, w_trapz, Py_plus, Py_minus, Pa_plus;
  double Gamma, D, V0, V1, Vz, Vz_bar,rho_aa,rho_aa_bar;
  double P0,P0_bar,Px,Px_bar,Py,Py_bar,Pz,Pz_bar;
  double f0,  mu_div_T, feq, feq_bar;
  double I_df0Pa_plus, I_VxPy_minus, I_f0Pa_plus,dPa_plusdt, n_plus;
  int idx;
  double dn_plusdT, dLdT, feq_plus;

  if (pqke->L_initial != 0.0)
    L = y[pqke->index_L]*_L_SCALE_;
  else
    L = 0.0;
  n_plus = y[pqke->index_n_plus];
  //L = 1e-2*(1.0-1.0/(exp((T*1e3-2.486022964891155)/0.09)+1.0))+1e-5*sin((T-0.050)*1e5);
  


  get_resonances_xi(T,L,pqke);

  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(pqke->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);
  //Set V1:
  if (pqke->is_electron == _TRUE_){
    pqke->V1 = -14.0*sqrt(2.0)*_PI_*_PI_/45.0*_G_F_/_M_W_/_M_W_*pow(T,5)*
      (1.0+0.25*(1.0-_SIN2_THETA_W_)*n_plus);
  }
  else{
    pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
      _G_F_/_M_Z_/_M_Z_*pow(T,5)*n_plus;
  }
  
  /** Integrated quantities needed. We integrate in x space: */
  I_VxPy_minus = 0.0;
  I_df0Pa_plus = 0.0;
   for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
    //Solving mu from L, using the chebyshev cubic root:
    mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
    //Distributions:
    feq_plus = 1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T));
    f0 = 1.0/(1.0+exp(x));

    Vx = pqke->Vx/x;

    Py_plus = y[pqke->index_Py+i]+y[pqke->index_Py_bar+i];
    Py_minus = y[pqke->index_Py+i]-y[pqke->index_Py_bar+i];
    Pa_plus = y[pqke->index_P0+i]+y[pqke->index_P0_bar+i]+
      y[pqke->index_Pz+i]+y[pqke->index_Pz_bar+i];
    dPa_plusdt = Vx*Py_plus+Gamma*(2.0*feq_plus/f0-Pa_plus);
    
    
    I_VxPy_minus += w_trapz*(x*x*Vx*Py_minus*f0);
    I_df0Pa_plus += w_trapz*(x*x*f0*dPa_plusdt);
  }
  dLdT = -1.0/(8.0*H*T*_ZETA3_)*I_VxPy_minus;
  pqke->dLdT = dLdT;
  dn_plusdT = -1.0/(3.0*H*T*_ZETA3_)*I_df0Pa_plus;


  /** Calculate RHS: */
  dy[pqke->index_L] = dLdT/_L_SCALE_;
  dy[pqke->index_n_plus] = dn_plusdT;

  if (fabs(L)<0.9*fabs(pqke->L_initial)){
    if (pqke->dLdT_approx == _FALSE_){
      pqke->mu = dy[pqke->index_L]/L;
      pqke->T1 = T;
      pqke->dLdT_approx = _TRUE_;
    }
    dy[pqke->index_L] = L*(pqke->mu+6.9e9*pow(pqke->T1-T,2)/_L_SCALE_);
  }


  /** All quantities defined on the grid: */

  
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Vz = V0+V1+VL;
    Vz_bar = V0+V1-VL;
  
    Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
    D = 0.5*Gamma;
    
    

    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    Px = y[pqke->index_Px+i];
    Px_bar = y[pqke->index_Px_bar+i];
    Py = y[pqke->index_Py+i];
    Py_bar = y[pqke->index_Py_bar+i];
    Pz = y[pqke->index_Pz+i];
    Pz_bar = y[pqke->index_Pz_bar+i];
  
    //Solving mu from L, using the chebyshev cubic root:
    mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
    //Distributions:
    feq = 1.0/(1.0+exp(x-mu_div_T));
    feq_bar = 1.0/(1.0+exp(x+mu_div_T));
    f0 = 1.0/(1.0+exp(x));
    
    rho_aa = 0.5*f0*(P0+Pz);
    rho_aa_bar = 0.5*f0*(P0_bar+Pz_bar);

    idx = pqke->index_P0+i;
    dy[idx] = -1.0/(H*T)*(Gamma/f0*(feq-rho_aa));

    idx = pqke->index_P0_bar+i;
    dy[idx] = -1.0/(H*T)*(Gamma/f0*(feq_bar-rho_aa_bar));

    idx = pqke->index_Px+i;
    dy[idx] = -1.0/(H*T)*(-Vz*Py-D*Px);

    idx = pqke->index_Px_bar+i;
    dy[idx] = -1.0/(H*T)*(-Vz_bar*Py_bar-D*Px_bar);
    
    idx = pqke->index_Py+i;
    dy[idx] = -1.0/(H*T)*(Vz*Px-Vx*Pz-D*Py);

    idx = pqke->index_Py_bar+i;
    dy[idx] = -1.0/(H*T)*(Vz_bar*Px_bar-Vx*Pz_bar-D*Py_bar);

    idx = pqke->index_Pz+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py+Gamma/f0*(feq-rho_aa));
    
    idx = pqke->index_Pz_bar+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_bar+Gamma/f0*(feq_bar-rho_aa_bar));
    if (i==-1){
      printf("dPzdT[40] = %g, dPz_bardT[4]=%g\n",
	     dy[pqke->index_Pz+i],dy[pqke->index_Pz_bar+i]);
    }

  }
  return _SUCCESS_;
}

int qke_derivs_fixed_grid_standard_approx(double T, 
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
  double x, w_trapz, Py_plus, Py_minus, Pa_plus;
  double Gamma, D, V0, V1, Vz, Vz_bar,rho_aa,rho_aa_bar;
  double P0,P0_bar,Px,Px_bar,Py,Py_bar,Pz,Pz_bar;
  double f0,  mu_div_T, feq, feq_bar;
  double I_df0Pa_plus, I_VxPy_minus, I_f0Pa_plus,dPa_plusdt, n_plus;
  int idx;
  double dn_plusdT, feq_plus,x2;

  double tol_dLdT=1e-12, dLdT, dLdT_old=0.0, dLdT_from_dy=0.0,dLdT_from_dy_old=0.0;
  int iter, max_iter=100, search_success;
  double I_dPzdT, dVzdT, dVz_bardT, s, s_bar,dsinalphadT,dsinalpha_bardT;
  double sinalphaPz_ini,sinalphaPz_bar_ini,I_Paplus,dPzdT,dPz_bardT;
  double gamma, gamma_bar, I_Paminus;
  double dlogVzdlogT, dlogVz_bardlogT, threshold = 100.0;
  int approx = _FALSE_, approx_bar = _FALSE_;
  //  L = y[pqke->index_L]*_L_SCALE_;
  //L = 1e-2*(1.0-1.0/(exp((T*1e3-2.486022964891155)/0.09)+1.0))+1e-5*sin((T-0.050)*1e5);

  I_Paplus = 0.0;
  I_Paminus = 0.0;
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    f0 = 1.0/(exp(x)+1.0);
    
    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    Pz = y[pqke->index_Pz+i];
    Pz_bar = y[pqke->index_Pz_bar+i];
    
    I_Paminus += w_trapz*x*x*f0*(P0-P0_bar+Pz-Pz_bar);
    I_Paplus += w_trapz*x2*f0*(P0+P0_bar+Pz+Pz_bar);
  }

  n_plus = I_Paplus/(3.0*_ZETA3_);
  L = pqke->L_fudge*I_Paminus/(8.0*_ZETA3_);
  //printf("L_fudge = %g, L=%g\n",pqke->L_fudge,L);
  
  /**
     get_resonances_xi(T,L,pqke);
  for(i=0; i<pqke->vres-1; i++){
    if (pqke->xi[0]>x_grid[i])
      break;
  }
  */

  //x = x_grid[i];
  //  dLdT = -1.0/(4.0*_ZETA3_)*x*x/(exp(x)+1.0)*y[pqke->index_Pz_bar+i+1]*pqke->dxidT[0];
  dLdT = pqke->dLdT;
  //dLdT = 0.0;

  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(pqke->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);
  
  search_success = _FALSE_;
  for(iter=1; iter<=max_iter;iter++){
    approx = _FALSE_; 
    approx_bar = _FALSE_;
    //printf("dLdT_guess = %g\n",dLdT);
     //Use dLdT to set dP0dT and find the dLdT resulting from this:
    //(We are searching for a self-consistent value of dLdT)

    /** Integrated quantities needed. We integrate in x space: */
  
    for (i=0; i<pqke->vres; i++){
      if (i==0)
	w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
      else if (i==pqke->vres-1)
	w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
      else
	w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
      x = x_grid[i];
      x2 = x*x;
      f0 = 1.0/(exp(x)+1.0);
      
      I_Paplus += w_trapz*x2*f0*
	(y[pqke->index_P0+i]+y[pqke->index_P0_bar+i]+
	 y[pqke->index_Pz+i]+y[pqke->index_Pz_bar+i]);
    }
    n_plus = I_Paplus/(3.0*_ZETA3_);
    //printf("n_plus = %g\n",n_plus);

    //Set V1:
    if (pqke->is_electron == _TRUE_){
      pqke->V1 = -14.0*sqrt(2.0)*_PI_*_PI_/45.0*_G_F_/_M_W_/_M_W_*pow(T,5)*
	(1.0+0.25*(1.0-_SIN2_THETA_W_)*n_plus);
    }
    else{
      pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
	_G_F_/_M_Z_/_M_Z_*pow(T,5)*n_plus;
    }

    /** All quantities defined on the grid: */
    I_dPzdT = 0.0;  
    for (i=0; i<pqke->vres; i++){
      x = x_grid[i];
      Vx = pqke->Vx/x;
      V0 = pqke->V0/x;
      V1 = pqke->V1*x;
      Vz = V0+V1+VL;
      Vz_bar = V0+V1-VL;

      Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
      D = 0.5*Gamma;
      
      //Perform integral for dLdT:
      if (i==0)
	w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
      else if (i==pqke->vres-1)
	w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
      else
	w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
      
      P0 = y[pqke->index_P0+i];
      P0_bar = y[pqke->index_P0_bar+i];
      Pz = y[pqke->index_Pz+i];
      Pz_bar = y[pqke->index_Pz_bar+i];
  
      //Solving mu from L, using the chebyshev cubic root:
      mu_div_T = -2*_PI_/sqrt(3.0)*
	sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
      //Distributions:
      feq = 1.0/(1.0+exp(x-mu_div_T));
      feq_bar = 1.0/(1.0+exp(x+mu_div_T));
      f0 = 1.0/(1.0+exp(x));
    
      rho_aa = 0.5*f0*(P0+Pz);
      rho_aa_bar = 0.5*f0*(P0_bar+Pz_bar);


      //Preparations for the dPzdT calculation:
      dVzdT = (-V0+5.0*V1+VL*(3.0+T/L*dLdT))/T;
      dVz_bardT = (-V0+5.0*V1-VL*(3.0+T/L*dLdT))/T;

      dlogVzdlogT = (-V0+5.0*V1+VL*(3.0+T/L*dLdT))/Vz;
      dlogVz_bardlogT = (-V0+5.0*V1-VL*(3.0+T/L*dLdT))/Vz_bar;
      /**
      if (dlogVzdlogT>threshold){
	approx = _TRUE_;
	dlogVzdlogT = threshold;
      }
      if (dlogVzdlogT<threshold){
	dlogVzdlogT = threshold;
	approx = _TRUE_;
      }
      
      if (dlogVz_bardlogT>threshold){
	dlogVz_bardlogT = threshold;
	approx_bar = _TRUE_;
      }
      if (dlogVz_bardlogT<-threshold){
	approx_bar = _TRUE_;
	dlogVz_bardlogT = -threshold;
      }
      */

      s = sqrt(Vx*Vx+Vz*Vz);
      s_bar = sqrt(Vx*Vx+Vz_bar*Vz_bar);
      dsinalphadT = Vx*Vx*Vz/pow(s,3)/T*(1.0+dlogVzdlogT);
      dsinalpha_bardT = Vx*Vx*Vz_bar/pow(s_bar,3)/T*(1.0+dlogVz_bardlogT);

      sinalphaPz_ini = pqke->sigma0Pz[i];
      sinalphaPz_bar_ini = pqke->sigma0Pz_bar[i];

      gamma = pow(1.0+pow(Vx/Vz,2),-0.5);
      gamma_bar = pow(1.0+pow(Vx/Vz_bar,2),-0.5);
      
      /**      
      dPzdT = Pz/T*(1.0-pow(Pz/sinalphaPz_ini,2))*
	(1.0+T/Vz*dVzdT);
      dPz_bardT = Pz_bar/T*(1.0-pow(Pz_bar/sinalphaPz_bar_ini,2))*
	(1.0+T/Vz_bar*dVz_bardT);
      */
      /**
      dPzdT = sinalphaPz_ini*gamma/T*
	(1.0-gamma*gamma)*(1+T/Vz*dVzdT);
      dPz_bardT = sinalphaPz_bar_ini*gamma_bar/T*
	(1.0-gamma_bar*gamma_bar)*(1+T/Vz_bar*dVz_bardT);
      */
      dPzdT = sinalphaPz_ini*dsinalphadT;
      dPz_bardT = sinalphaPz_bar_ini*dsinalpha_bardT;

      //Set integral for dLdT:
      I_dPzdT +=w_trapz*f0*x*x*(dPzdT-dPz_bardT);
      
      idx = pqke->index_P0+i;
      dy[idx] = -1.0/(H*T)*(Gamma/f0*(feq-rho_aa));

      idx = pqke->index_P0_bar+i;
      dy[idx] = -1.0/(H*T)*(Gamma/f0*(feq_bar-rho_aa_bar));
      
      idx = pqke->index_Pz+i;
      if (approx==_TRUE_)
	dy[idx] = dVzdT*sinalphaPz_ini*pow(Vx,2)/pow(s,3)+dy[pqke->index_P0+i];
      else
	dy[idx] = dPzdT+dy[pqke->index_P0+i];

      idx = pqke->index_Pz_bar+i;
      if (approx_bar == _TRUE_)
	dy[idx] =dVz_bardT*sinalphaPz_bar_ini*pow(Vx,2)/pow(s_bar,3)+dy[pqke->index_P0_bar+i];
      else
	dy[idx] =dPz_bardT+dy[pqke->index_P0_bar+i];
      
      /**
      idx = pqke->index_Pz+i;
      dy[idx] = Pz/T*(1.0-pow(Pz/sinalphaPz_ini,2))*(1+dlogVzdlogT)+
	dy[pqke->index_P0+i];

      idx = pqke->index_Pz_bar+i;
      dy[idx] =Pz_bar/T*(1.0-pow(Pz_bar/sinalphaPz_bar_ini,2))*(1+dlogVz_bardlogT)+
	dy[pqke->index_P0_bar+i];
      */
      if (i==-1){
      printf("dPzdT[40] = %g, dPz_bardT[4]=%g, dPzdT[40] = %g, dPz_bardT[4]=%g\n",
	     dPzdT,dPz_bardT,Pz/T,Pz_bar/T);
      printf("Pz/sinalphaPz: %g, %g\n",
	     Pz*Vz/s/sinalphaPz_ini,
	     Pz_bar*Vz_bar/s_bar/sinalphaPz_bar_ini);
      }
    }
    dLdT_from_dy = I_dPzdT/(8.0*_ZETA3_);
  
    if (abs(1.0-dLdT/dLdT_from_dy)<tol_dLdT){
      dLdT = 0.5*(dLdT+dLdT_from_dy);
      pqke->dLdT = dLdT;
      search_success = _TRUE_;
      if (iter>3)
	printf("Had to work a bit, iter = %d\n",iter);
      break;
    }
    //Set new guess for dLdT:
    if (iter==1){
      dLdT = dLdT_from_dy;
    }
    else{
      dLdT += (dLdT-dLdT_from_dy)*(dLdT_old-dLdT)/
	(dLdT-dLdT_from_dy-dLdT_old+dLdT_from_dy_old);
    }
    dLdT_old = dLdT;
    dLdT_from_dy_old = dLdT_from_dy;
  }
  if (search_success == _TRUE_){
    /** Calculate RHS: */
    //dy[pqke->index_L] = dLdT/_L_SCALE_;
    return _SUCCESS_;
  }
  else{
    return _FAILURE_;
  }
}

int qke_derivs_fixed_grid_standard_approx2(double T, 
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
  double x, w_trapz, Py_plus, Py_minus, Pa_plus;
  double Gamma, D, V0, V1, Vz, Vz_bar,rho_aa,rho_aa_bar;
  double P0,P0_bar,Px,Px_bar,Py,Py_bar,Pz,Pz_bar;
  double f0,  mu_div_T, feq, feq_bar;
  double I_df0Pa_plus, I_VxPy_minus, I_f0Pa_plus,dPa_plusdt, n_plus;
  int idx;
  double dn_plusdT, feq_plus,x2;

  double tol_dLdT=1e-12, tol_L = 1e-15, tol_dif_L=1e-12, dLdT, dLdT_old=0.0, dLdT_from_dy=0.0,dLdT_from_dy_old=0.0;
  int iter, max_iter=500, search_success=_FALSE_;
  double I_dPzdT, dVzdT, dVz_bardT, s, s_bar,dsinalphadT,dsinalpha_bardT;
  double sinalphaPz_ini,sinalphaPz_bar_ini,I_Paplus,dPzdT,dPz_bardT;
  double gamma, gamma_bar, I_Paminus, I_Pzpure;
  double dlogVzdlogT, dlogVz_bardlogT, threshold = 100.0;
  int approx = _FALSE_, approx_bar = _FALSE_;
  double L_left, L_right,Lnew,Lz;
  double P0_ini,P0_bar_ini,delta_Pz,delta_Pz_bar,Pz_pure,Pz_bar_pure;
  double delta_Pz_ini,delta_Pz_bar_ini;
  //  L = y[pqke->index_L]*_L_SCALE_;
  //L = 1e-2*(1.0-1.0/(exp((T*1e3-2.486022964891155)/0.09)+1.0))+1e-5*sin((T-0.050)*1e5);
  int sign_f_left, sign_f_right, factor;
  double dif_left,dif_right, dif_L;
  double L_last, steps = 100, sign_L_last, step_L = 1e-8;
  int bracket_success;
  n_plus = pqke->n_plus;

  pqke->Vx = pqke->delta_m2/(2.0*T)*sin(2.0*pqke->theta_zero);
  pqke->V0 = -pqke->delta_m2/(2.0*T)*cos(2.0*pqke->theta_zero);
  //Set V1:
  if (pqke->is_electron == _TRUE_){
    pqke->V1 = -14.0*sqrt(2.0)*_PI_*_PI_/45.0*_G_F_/_M_W_/_M_W_*pow(T,5)*
      (1.0+0.25*(1.0-_SIN2_THETA_W_)*n_plus);
  }
  else{
    pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
      _G_F_/_M_Z_/_M_Z_*pow(T,5)*n_plus;
  }

  //Find L:  
  L_last  = pqke->Lz+pqke->L0;
  if (L_last<0)
    sign_L_last = -1;
  else
    sign_L_last = 1;

  bracket_success = _FALSE_;
  L_left = L_last;
  L_right = L_last;
  //Bracket root:
  for (factor=1; factor<100000; factor++){
    L_left = 0.99*L_left-sign_L_last*1e-12;
    L_right *=(1.0+5e-5);
  
    dif_left = dif_L_from_exact(T, y, L_left, pqke);
    dif_right = dif_L_from_exact(T, y, L_right, pqke);
    if (dif_left*dif_right<0.0){
      bracket_success = _TRUE_;
      break;
    }
  }
  if (bracket_success == _FALSE_){
    printf("Bracketing failed! dif_L = %g, dif_R = %g\n",dif_left,dif_right);
    printf("L_left = %g, L_right = %g, L_last=%g, steps=%d\n",
	   L_left,L_right,L_last,factor);
    printf("Temeprature=%g MeV \n",T*1e3);
    for (i=0,L=L_last/pow(0.99,50); i<(factor+50); i++){
      L = 0.99*L-sign_L_last*1e-12;
      printf("%.16e %.16e\n",L,dif_L_from_exact(T,y,L,pqke));
    }
    printf("Rows: %d\n",factor+50);
    return _FAILURE_;
  }
  else{
    printf("Bracketing success! After %d steps..\n",factor);
  }
  printf("L_last=%g, bracket:[%g; %g], difL=%g, difR=%g\n",
	 L_last,L_left,L_right,dif_left, dif_right);
  do {
    L = 0.5*(L_left+L_right);
    dif_L = dif_L_from_exact(T, y, L, pqke);
    printf("%g ",fabs(L_right-L_left));
    if (dif_L*dif_right>0.0)
      L_right = L;
    else
      L_left = L;
  } while(fabs(L_right-L_left)>tol_dif_L);
  pqke->Lz = L-pqke->L0;

  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  I_Paplus = 0.0; 
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    f0 = 1.0/(exp(x)+1.0);

    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Vz = V0+V1+VL;
    Vz_bar = V0+V1-VL;
    gamma = Vz/sqrt(Vx*Vx+Vz*Vz);
    gamma_bar = Vz_bar/sqrt(Vx*Vx+Vz_bar*Vz_bar);
    
    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    
    delta_Pz_ini = pqke->delta_Pz_ini[i];
    delta_Pz_bar_ini = pqke->delta_Pz_bar_ini[i];
    P0_ini = pqke->P0_ini[i];
    P0_bar_ini = pqke->P0_bar_ini[i];
    
    delta_Pz = delta_Pz_ini+P0-P0_ini;
    delta_Pz_bar = delta_Pz_bar_ini+P0_bar-P0_bar_ini;
    
    Pz_pure = gamma*pqke->sigma0Pz[i];
    Pz_bar_pure = gamma_bar*pqke->sigma0Pz_bar[i];
      
    Pz = Pz_pure+delta_Pz;
    Pz_bar = Pz_bar_pure+delta_Pz_bar;
      
    I_Pzpure += w_trapz*x*x*f0*(Pz_pure-Pz_bar_pure);
    I_Paplus += w_trapz*x*x*f0*(P0+P0_bar+Pz+Pz_bar);
  }

  pqke->n_plus = I_Paplus/(3.0*_ZETA3_);
  
  //  printf("L = %g\n",L);
  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(pqke->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  pqke->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = pqke->VL;
  //Set V1:
  if (pqke->is_electron == _TRUE_){
    pqke->V1 = -14.0*sqrt(2.0)*_PI_*_PI_/45.0*_G_F_/_M_W_/_M_W_*pow(T,5)*
      (1.0+0.25*(1.0-_SIN2_THETA_W_)*pqke->n_plus);
  }
  else{
    pqke->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
      _G_F_/_M_Z_/_M_Z_*pow(T,5)*pqke->n_plus;
  }
  
  /** All quantities defined on the grid: */
  for (i=0; i<pqke->vres; i++){
    x = x_grid[i];
    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Vz = V0+V1+VL;
    Vz_bar = V0+V1-VL;

    Gamma = pqke->C_alpha*_G_F_*_G_F_*x*pow(T,5);
    D = 0.5*Gamma;
    
    gamma = Vz/sqrt(Vx*Vx+Vz*Vz);
    gamma_bar = Vz_bar/sqrt(Vx*Vx+Vz_bar*Vz_bar);

    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    
    delta_Pz_ini = pqke->delta_Pz_ini[i];
    delta_Pz_bar_ini = pqke->delta_Pz_bar_ini[i];
    P0_ini = pqke->P0_ini[i];
    P0_bar_ini = pqke->P0_bar_ini[i];
    
    delta_Pz = delta_Pz_ini+P0-P0_ini;
    delta_Pz_bar = delta_Pz_bar_ini+P0_bar-P0_bar_ini;

    Pz_pure = gamma*pqke->sigma0Pz[i];
    Pz_bar_pure = gamma_bar*pqke->sigma0Pz_bar[i];
      
    Pz = Pz_pure+delta_Pz;
    Pz_bar = Pz_bar_pure+delta_Pz_bar;
    
    //Solving mu from L, using the chebyshev cubic root:
    mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
    //Distributions:
    feq = 1.0/(1.0+exp(x-mu_div_T));
    feq_bar = 1.0/(1.0+exp(x+mu_div_T));
    f0 = 1.0/(1.0+exp(x));
    
    rho_aa = 0.5*f0*(P0+Pz);
    rho_aa_bar = 0.5*f0*(P0_bar+Pz_bar);

    idx = pqke->index_P0+i;
    dy[idx] = -1.0/(H*T)*(Gamma/f0*(feq-rho_aa));

    idx = pqke->index_P0_bar+i;
    dy[idx] = -1.0/(H*T)*(Gamma/f0*(feq_bar-rho_aa_bar));
  }
  return _SUCCESS_;
}


double dif_L_from_exact(double T, 
			double *y,
			double L_in,
			qke_param *pqke){
  double VL, I_Pzpure, w_trapz, x, f0, Vx, V0, Vz, Vz_bar, V1;
  double gamma, gamma_bar, P0, P0_bar, delta_Pz_ini, delta_Pz_bar_ini;
  double P0_ini, P0_bar_ini, delta_Pz, delta_Pz_bar, Pz_pure, Pz_bar_pure;
  double Pz, Pz_bar;
  double *x_grid = pqke->x_grid;
  int i;

  VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L_in;
  I_Pzpure = 0.0;
  for (i=0; i<pqke->vres; i++){
    if (i==0)
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i]);
    else if (i==pqke->vres-1)
      w_trapz = 0.5*(x_grid[i]-x_grid[i-1]);
    else
      w_trapz = 0.5*(x_grid[i+1]-x_grid[i-1]);
    x = x_grid[i];
    f0 = 1.0/(exp(x)+1.0);

    Vx = pqke->Vx/x;
    V0 = pqke->V0/x;
    V1 = pqke->V1*x;
    Vz = V0+V1+VL;
    Vz_bar = V0+V1-VL;
    gamma = Vz/sqrt(Vx*Vx+Vz*Vz);
    gamma_bar = Vz_bar/sqrt(Vx*Vx+Vz_bar*Vz_bar);

    P0 = y[pqke->index_P0+i];
    P0_bar = y[pqke->index_P0_bar+i];
    
    delta_Pz_ini = pqke->delta_Pz_ini[i];
    delta_Pz_bar_ini = pqke->delta_Pz_bar_ini[i];
    P0_ini = pqke->P0_ini[i];
    P0_bar_ini = pqke->P0_bar_ini[i];

    delta_Pz = delta_Pz_ini+P0-P0_ini;
    delta_Pz_bar = delta_Pz_bar_ini+P0_bar-P0_bar_ini;

    Pz_pure = gamma*pqke->sigma0Pz[i];
    Pz_bar_pure = gamma_bar*pqke->sigma0Pz_bar[i];
    
    Pz = Pz_pure+delta_Pz;
    Pz_bar = Pz_bar_pure+delta_Pz_bar;
    
    I_Pzpure += w_trapz*x*x*f0*(Pz_pure-Pz_bar_pure);
  }
  return L_in-(pqke->L0+I_Pzpure/(8.0*_ZETA3_));
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
