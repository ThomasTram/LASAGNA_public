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
  pqke->Tvec = malloc(sizeof(double)*pqke->Tres);
  pqke->xi = malloc(sizeof(double)*pqke->Nres);
  pqke->ui = malloc(sizeof(double)*pqke->Nres);
  pqke->vi = malloc(sizeof(double)*pqke->Nres);
  pqke->duidT = malloc(sizeof(double)*pqke->Nres);
  pqke->dvidT = malloc(sizeof(double)*(pqke->Nres));
  pqke->duidx = malloc(sizeof(double)*pqke->Nres);
  pqke->dxidT = malloc(sizeof(double)*pqke->Nres);
  pqke->a = malloc(sizeof(double)*pqke->Nres);
  pqke->y_0 = malloc(sizeof(double)*(1+pqke->Nres));
  pqke->maxstep = malloc(sizeof(double)*(1+pqke->Nres));
  pqke->x_grid = malloc(sizeof(double)*pqke->vres);
  pqke->u_grid = malloc(sizeof(double)*pqke->vres);
  pqke->v_grid = malloc(sizeof(double)*pqke->vres);
  pqke->dvdu_grid = malloc(sizeof(double)*pqke->vres);
  pqke->dudT_grid = malloc(sizeof(double)*pqke->vres);
  pqke->mat = malloc(sizeof(double*)*(pqke->Nres+2));
  for (i=0; i<(pqke->Nres+2); i++) 
    pqke->mat[i] = malloc(sizeof(double)*(pqke->Nres + 2));
  pqke->vv = malloc(sizeof(double)*(pqke->Nres + 2));
  pqke->indx = malloc(sizeof(int)*(pqke->Nres + 2));

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

int free_qke_param(qke_param *pqke){
  int i;
  background_free_dof(&(pqke->pbs));
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
	
	if (xi[1]<xi[0])
	  printf("x1 = %g, x2 = %g, ???\n",xi[0],xi[1]);

	frac = (2.0*phi*phi+3*chi)/sqrt(phi*phi+chi);
	dxidT[0] = -(frac+2.0*phi)/T;
	dxidT[1] = -(frac-2.0*phi)/T;

	// Protect against possible seg fault:
	for (i=0; i<param->Nres; i++){
	  if (xi[i]<param->xmin) {
	    printf("Note: Resonance at T=%g MeV is lower than xmin=%g. (It is %g.)\n",
		   T*1e3,param->xmin,xi[i]);
	    xi[i] = param->xmin;
	    dxidT[i] = 0.0;
	  }
	  else if(xi[i]>param->xmax){
	    printf("Note: Resonance at T=%g MeV is higher than xmax=%g. (It is %g.)\n",
		   T*1e3,param->xmax,xi[i]);
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

int qke_initial_conditions(double Ti, double *y, qke_param *pqke){
  /** Set initial conditions at temperature Ti: */
  int i, evolve_vi;
  double mu_div_T,x;
  ErrorMsg error_message;
  double *dy;
  double Vx,D,Vz,Vz_bar,Px,Py,Px_bar,Py_bar;
  //double sterile_ratio = 0.0;
  //Calculate chemical potential from initial L:
  mu_div_T = -2*_PI_/sqrt(3.0)*
      sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*pqke->L_initial/pow(_PI_,3)));
  pqke->mu_div_T_initial = mu_div_T;

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
  mat_add_matrix(outf,"L_initial",miDOUBLE,1,2,&handle);
  mat_add_matrix(outf,"delta_m2",miDOUBLE,1,2,&handle);
  mat_add_matrix(outf,"is_electron",miINT32,1,1,&handle);
  mat_add_matrix(outf,"Tres_vres",miINT32,1,2,&handle);
  mat_add_matrix(outf,"xmin_xext_xmax",miDOUBLE,1,3,&handle);
  mat_add_matrix(outf,"alpha_rs",miDOUBLE,1,2,&handle);
  //Write parameters and stuff we know in advance:
  //Temperature - we could write it here, but its nice to have non-computed
  //T values = 0:
  //mat_write_data(outf,"T",pqke->Tvec,0,Tres);
  mat_write_data(outf,"delta_m2",&(pqke->delta_m2),0,1);
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


int qke_print_variables(double t,
			double *y,
			double *dy,
			void *param,
			ErrorMsg error_message){
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
  double vN,uN,u1,v1,alpha,wi;
  double gentr,H;
  double lu_sgn;
  double daidT,dbdT;
  double Vx, Vxp1, VL;
  double x,xp1;
  double Gamma, D, V0, V1, Pa_plus, Pa_minus, Ps_plus, Ps_minus;
  double Px_plus, Px_minus, Py_plus, Py_minus, f0, f0p1,  mu_div_T;
  double Py_minusp1,Pa_plusp1,PsPs,PsPsp1,I_rho_ss;
  double feq_plus, feq_minus, feq_ini, I_VxPy_minus, I_f0Pa_plus;
  double dudTdvdu, delta_v;
  double rs, mu_div_T_ini;
  double gamma_j, beta_j;
  int idx,idx_step_l,idx_step_r;

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
    y_0[0] = 1.0 - alpha;
    maxstep[0] = 100.0;
    for (i=1; i<=pqke->Nres; i++){
      y_0[i] = ui[i-1];
      maxstep[i] = 0.1;
    }
    //Find parametrisation parameters vi, a and b using Newton:
    lasagna_call(Newton(nonlinear_rhs,y_0,pqke,maxstep,tol_newton,&niter,100,pqke->Nres+1,error_message),error_message,error_message);
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
  for (i=1; i<pqke->Nres; i++){
    /** Find the splitting of v and store it temporarily in pqke->indx.
	We use the fact that v is uniform.
    */
    wi = 0.5*(vi[i-1]+vi[i]); //Weighted average
    pqke->indx[i] = (int)(wi/(v_grid[1]-v_grid[0]));
  }
  pqke->indx[0] = 0;
  pqke->indx[pqke->Nres] = pqke->vres;
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
  for (i=0; i<pqke->vres-1; i++){
    x = pqke->x_grid[i];
    xp1 = pqke->x_grid[i+1];
    f0 = 1.0/(1.0+exp(x));
    f0p1 = 1.0/(1.0+exp(xp1));
    Vx = pqke->Vx/x;
    Vxp1 = pqke->Vx/xp1;
    Py_minus = y[pqke->index_Py_minus+i];
    Py_minusp1 = y[pqke->index_Py_minus+i+1];
    Pa_plus = y[pqke->index_Pa_plus+i];
    Pa_plusp1 = y[pqke->index_Pa_plus+i+1];
    PsPs = y[pqke->index_Ps_plus+i]+y[pqke->index_Ps_minus+i];
    PsPsp1 = y[pqke->index_Ps_plus+i+1]+y[pqke->index_Ps_minus+i+1];
    
    I_VxPy_minus += 0.5*(xp1-x)*(x*x*Vx*Py_minus*f0+xp1*xp1*Vxp1*Py_minusp1*f0p1);
    I_f0Pa_plus += 0.5*(xp1-x)*(x*x*f0*Pa_plus+xp1*xp1*f0p1*Pa_plusp1);
    I_rho_ss += 0.5*(xp1-x)*(x*x*f0*PsPs+xp1*xp1*f0p1*PsPsp1);
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
  dy[pqke->index_L] = -1.0/(8.0*H*T*_ZETA3_)*I_VxPy_minus/_L_SCALE_;
  /** All quantities defined on the grid: */

  // Set perhaps flow of grid:
  if (pqke->evolve_vi == _TRUE_){
    //dy[pqke->index_b] = dbdT;
    for (i=0; i<pqke->Nres; i++){
      dy[pqke->index_vi+i] = dvidT[i];
    }
  }

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
    mu_div_T_ini = pqke->mu_div_T_initial;
    //Distributions:
    feq_plus = 1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T));
    feq_minus = 1.0/(1.0+exp(x-mu_div_T))-1.0/(1.0+exp(x+mu_div_T));
    f0 = 1.0/(1.0+exp(x));
    feq_ini = 1.0/(1.0+exp(x-mu_div_T_ini));

    dudTdvdu = dudT_grid[i]*dvdu_grid[i];
    //Define index steps for calculating derivatives:
    if (i==0)
      idx_step_l = 0;
    else
      idx_step_l = -1;
    if (i==pqke->vres-1)
      idx_step_r = 0;
    else
      idx_step_r = 1;
    delta_v = v_grid[i+idx_step_r]-v_grid[i+idx_step_l];

    idx = pqke->index_Pa_plus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_plus+Gamma*(2.0*feq_plus/f0-Pa_plus))+
      dudTdvdu*(y[idx+idx_step_r]-y[idx+idx_step_l])/delta_v;

    idx = pqke->index_Pa_minus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_minus+Gamma*(2.0*feq_minus/f0-Pa_minus))+
      dudTdvdu*(y[idx+idx_step_r]-y[idx+idx_step_l])/delta_v;

    idx = pqke->index_Ps_plus+i;
    dy[idx] = 1.0/(H*T)*
      (Vx*Py_plus-rs*Gamma*(1.0/(6.0*_ZETA3_)*I_rho_ss*feq_plus-0.5*f0*Ps_plus))+
      dudTdvdu*(y[idx+idx_step_r]-y[idx+idx_step_l])/delta_v;
    
    idx = pqke->index_Ps_minus+i;
    dy[idx] = 1.0/(H*T)*
      (Vx*Py_minus-rs*Gamma*f0*(1.0-0.5*(Pa_minus+Ps_minus)))+
       dudTdvdu*(y[idx+idx_step_r]-y[idx+idx_step_l])/delta_v;
    
    idx = pqke->index_Px_plus+i;
    dy[idx] = 1.0/(H*T)*((V0+V1)*Py_plus+VL*Py_minus+D*Px_plus)+
       dudTdvdu*(y[idx+idx_step_r]-y[idx+idx_step_l])/delta_v;

    idx = pqke->index_Px_minus+i;
    dy[idx] = 1.0/(H*T)*((V0+V1)*Py_minus+VL*Py_plus+D*Px_minus)+
       dudTdvdu*(y[idx+idx_step_r]-y[idx+idx_step_l])/delta_v;

    idx = pqke->index_Py_plus+i;
    dy[idx] = 1.0/(H*T)*
      (-(V0+V1)*Px_plus-VL*Px_minus+0.5*Vx*(Pa_plus-Ps_plus)+D*Py_plus)+ 
      dudTdvdu*(y[idx+idx_step_r]-y[idx+idx_step_l])/delta_v;
   
    idx = pqke->index_Py_minus+i;
    dy[idx] = 1.0/(H*T)*
      (-(V0+V1)*Px_minus-VL*Px_plus+0.5*Vx*(Pa_minus-Ps_minus)+D*Py_minus)+
       dudTdvdu*(y[idx+idx_step_r]-y[idx+idx_step_l])/delta_v;
  }
  return _SUCCESS_;
}
