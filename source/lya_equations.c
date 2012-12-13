/** @file lya_equations.c 
 * Thomas Tram Bülow, 25.08.2011    
 */
#include "common.h"
#include "lya_equations.h"
int init_lya_param(lya_param *plya){
  int i,j,k,idx,nz;
  size_t neq;
  double k1,k2;
  double Nres, vres, Tres;
  int **J;
  ErrorMsg errmsg;
  lya_param *paw_derivs;
  double *Ax;
  Nres = plya->Nres;
  vres = plya->vres;
  Tres = plya->Tres;
  plya->Tvec = malloc(sizeof(double)*Tres);
  plya->xi = malloc(sizeof(double)*Nres);
  plya->ui = malloc(sizeof(double)*Nres);
  plya->vi = malloc(sizeof(double)*Nres);
  plya->duidT = malloc(sizeof(double)*Nres);
  plya->dvidT = malloc(sizeof(double)*Nres);
  plya->duidx = malloc(sizeof(double)*Nres);
  plya->dxidT = malloc(sizeof(double)*Nres);
  plya->a = malloc(sizeof(double)*Nres);
  plya->y_0 = malloc(sizeof(double)*(1+Nres));
  plya->maxstep = malloc(sizeof(double)*(1+Nres));
  plya->x_grid = malloc(sizeof(double)*vres);
  plya->u_grid = malloc(sizeof(double)*vres);
  plya->v_grid = malloc(sizeof(double)*vres);
  plya->dvdu_grid = malloc(sizeof(double)*vres);
  plya->dudT_grid = malloc(sizeof(double)*vres);
  plya->mat = malloc(sizeof(double*)*(Nres+2));
  for (i=0; i<(Nres+2); i++) 
    plya->mat[i] = malloc(sizeof(double)*(Nres + 2));
  plya->vv = malloc(sizeof(double)*(Nres + 2));
  plya->indx = malloc(sizeof(int)*(Nres + 2));

  //Do some calculations for the u(x) mapping:
  k1 = plya->xmin/plya->xext;
  k2 = plya->xmax/plya->xext;
  plya->eps2 = (1+k1)/(k2-k1);
  plya->eps1 = k1*(1+plya->eps2);

  //Some secondary initialisations:
  for(i=0; i<vres; i++){
    plya->v_grid[i] = plya->v_left+
      i*(plya->v_right-plya->v_left)/(vres-1.0); 
  }
  for(i=0; i<Tres; i++){
    plya->Tvec[i] = plya->T_initial+
      i*(plya->T_final-plya->T_initial)/(Tres-1.0); 
  }
  if (plya->is_electron == _TRUE_){
     plya->C_alpha = 1.27;
  }
  else{
     plya->C_alpha = 0.92;
  }
  plya->guess_exists = _FALSE_;
  
  //Set up the indices:
  idx = 0;
  plya->index_L = idx;
  idx++;
  plya->index_Pa_plus = idx;
  idx+=plya->vres;
  plya->index_Pa_minus =idx;
  idx +=plya->vres;
  plya->index_Ps_plus = idx;
  idx +=plya->vres;
  plya->index_Ps_minus = idx;
  idx +=plya->vres;
  plya->index_Px_plus = idx;
  idx +=plya->vres;
  plya->index_Px_minus = idx;
  idx +=plya->vres;
  plya->index_Py_plus = idx;
  idx +=plya->vres;
  plya->index_Py_minus = idx;
  idx +=plya->vres;

  //Last thing to do:
  neq = idx;
  plya->neq = neq;
  plya->index_v = neq;

  //Pattern for Jacobi matrix:
  plya->Ap = malloc(sizeof(int)*(2*neq+1));
  plya->Ai = malloc(sizeof(int)*2*neq*2*neq);

  //Construct Jacobian pattern:
  J = malloc(sizeof(int *)*2*neq);
  J[0] = calloc(2*neq*2*neq,sizeof(int)); 
  for (i=1; i<2*neq; i++)
    J[i] = J[i-1]+2*neq;
  

  /*
    J = [ J_QKE     0 ]
        [ J_X   J_QKE ]
    as an conservative choice set J_X = J_QKE.
  */

  /*
    Setup the upper left part that is used to solve the QKE.
  */

  //Establish diagonal (just to be sure)
  for (i=0; i<neq; i++)
    J[i][i] = 1;
  //Everything couples to L through the parametrisation:
  for (i=0; i<neq; i++)
    J[i][plya->index_L] = 1;


  //F-L dependence:
  for (i=0; i<vres; i++)
    J[plya->index_L][plya->index_Py_minus+i] = 1;
  //Loop over grid:
  for (j=0; j<vres; j++){
    //F-Pa_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Pa_plus+j][plya->index_Pa_plus+k] = 1;
    J[plya->index_Pa_plus+j][plya->index_Py_plus+j] = 1;
    //F-Pa_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Pa_minus+j][plya->index_Pa_minus+k] = 1;
    J[plya->index_Pa_minus+j][plya->index_Py_minus+j] = 1;
    //F-Ps_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Ps_plus+j][plya->index_Ps_plus+k] = 1;
    J[plya->index_Ps_plus+j][plya->index_Py_plus+j] = 1;
    //R_nus_plus dependence:
    for (k=0; k<vres; k++){
      J[plya->index_Ps_plus+j][plya->index_Ps_minus+k] = 1;
      J[plya->index_Ps_plus+j][plya->index_Ps_plus+k] = 1;
    }
    //F-Ps_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Ps_minus+j][plya->index_Ps_minus+k] = 1;
    J[plya->index_Ps_minus+j][plya->index_Py_minus+j] = 1;
    //R_nus_minus dependence:
    J[plya->index_Ps_minus+j][plya->index_Pa_minus+j] = 1;
    //F-Px_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Px_plus+j][plya->index_Px_plus+k] = 1;
    J[plya->index_Px_plus+j][plya->index_L] = 1;
    for (k=0; k<vres; k++)
      J[plya->index_Px_plus+j][plya->index_Pa_plus+k] = 1;
    J[plya->index_Px_plus+j][plya->index_Py_plus+j] = 1;
    J[plya->index_Px_plus+j][plya->index_Py_minus+j] = 1;
    //F-Px_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Px_minus+j][plya->index_Px_minus+k] = 1;
    J[plya->index_Px_minus+j][plya->index_L] = 1;
    for (k=0; k<vres; k++)
      J[plya->index_Px_minus+j][plya->index_Pa_plus+k] = 1;
    J[plya->index_Px_minus+j][plya->index_Py_plus+j] = 1;
    J[plya->index_Px_minus+j][plya->index_Py_minus+j] = 1;
    //F-Py_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Py_plus+j][plya->index_Py_plus+k] = 1;
    J[plya->index_Py_plus+j][plya->index_L] = 1;
    for (k=0; k<vres; k++)
      J[plya->index_Py_plus+j][plya->index_Pa_plus+k] = 1;
    J[plya->index_Py_plus+j][plya->index_Ps_plus+j] = 1;
    J[plya->index_Py_plus+j][plya->index_Px_plus+j] = 1;
    J[plya->index_Py_plus+j][plya->index_Px_minus+j] = 1;
    J[plya->index_Py_plus+j][plya->index_Py_plus+j] = 1;
    //F-Py_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Py_minus+j][plya->index_Py_minus+k] = 1;
    J[plya->index_Py_minus+j][plya->index_L] = 1;
    J[plya->index_Py_minus+j][plya->index_Pa_minus+j] = 1;
    for (k=0; k<vres; k++)
      J[plya->index_Py_minus+j][plya->index_Pa_plus+k] = 1;
    J[plya->index_Py_minus+j][plya->index_Ps_minus+j] = 1;
    J[plya->index_Py_minus+j][plya->index_Px_plus+j] = 1;
    J[plya->index_Py_minus+j][plya->index_Px_minus+j] = 1;
    J[plya->index_Py_minus+j][plya->index_Py_minus+j] = 1;
  }


  /*
    Setup the lower right corner that relates v to itself.
  */

  //Establish diagonal (just to be sure)
  for (i=0; i<neq; i++)
    J[i+neq][i+neq] = 1;
  //Everything couples to L through the parametrisation:
  for (i=0; i<neq; i++)
    J[i+neq][plya->index_L+neq] = 1;


  //F-L dependence:
  for (i=0; i<vres; i++)
    J[plya->index_L+neq][plya->index_Py_minus+i+neq] = 1;
  //Loop over grid:
  for (j=0; j<vres; j++){
    //F-Pa_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Pa_plus+j+neq][plya->index_Pa_plus+k+neq] = 1;
    J[plya->index_Pa_plus+j+neq][plya->index_Py_plus+j+neq] = 1;
    //F-Pa_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Pa_minus+j+neq][plya->index_Pa_minus+k+neq] = 1;
    J[plya->index_Pa_minus+j+neq][plya->index_Py_minus+j+neq] = 1;
    //F-Ps_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Ps_plus+j+neq][plya->index_Ps_plus+k+neq] = 1;
    J[plya->index_Ps_plus+j+neq][plya->index_Py_plus+j+neq] = 1;
    //R_nus_plus dependence:
    for (k=0; k<vres; k++){
      J[plya->index_Ps_plus+j+neq][plya->index_Ps_minus+k+neq] = 1;
      J[plya->index_Ps_plus+j+neq][plya->index_Ps_plus+k+neq] = 1;
    }
    //F-Ps_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Ps_minus+j+neq][plya->index_Ps_minus+k+neq] = 1;
    J[plya->index_Ps_minus+j+neq][plya->index_Py_minus+j+neq] = 1;
    //R_nus_minus dependence:
    J[plya->index_Ps_minus+j+neq][plya->index_Pa_minus+j+neq] = 1;
    //F-Px_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Px_plus+j+neq][plya->index_Px_plus+k+neq] = 1;
    J[plya->index_Px_plus+j+neq][plya->index_L+neq] = 1;
    for (k=0; k<vres; k++)
      J[plya->index_Px_plus+j+neq][plya->index_Pa_plus+k+neq] = 1;
    J[plya->index_Px_plus+j+neq][plya->index_Py_plus+j+neq] = 1;
    J[plya->index_Px_plus+j+neq][plya->index_Py_minus+j+neq] = 1;
    //F-Px_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Px_minus+j+neq][plya->index_Px_minus+k+neq] = 1;
    J[plya->index_Px_minus+j+neq][plya->index_L+neq] = 1;
    for (k=0; k<vres; k++)
      J[plya->index_Px_minus+j+neq][plya->index_Pa_plus+k+neq] = 1;
    J[plya->index_Px_minus+j+neq][plya->index_Py_plus+j+neq] = 1;
    J[plya->index_Px_minus+j+neq][plya->index_Py_minus+j+neq] = 1;
    //F-Py_plus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Py_plus+j+neq][plya->index_Py_plus+k+neq] = 1;
    J[plya->index_Py_plus+j+neq][plya->index_L+neq] = 1;
    for (k=0; k<vres; k++)
      J[plya->index_Py_plus+j+neq][plya->index_Pa_plus+k+neq] = 1;
    J[plya->index_Py_plus+j+neq][plya->index_Ps_plus+j+neq] = 1;
    J[plya->index_Py_plus+j+neq][plya->index_Px_plus+j+neq] = 1;
    J[plya->index_Py_plus+j+neq][plya->index_Px_minus+j+neq] = 1;
    J[plya->index_Py_plus+j+neq][plya->index_Py_plus+j+neq] = 1;
    //F-Py_minus[j] dependence:
    for (k=max(0,j-1); k<min(vres,j+2); k++)
      J[plya->index_Py_minus+j+neq][plya->index_Py_minus+k+neq] = 1;
    J[plya->index_Py_minus+j+neq][plya->index_L+neq] = 1;
    J[plya->index_Py_minus+j+neq][plya->index_Pa_minus+j+neq] = 1;
    for (k=0; k<vres; k++)
      J[plya->index_Py_minus+j+neq][plya->index_Pa_plus+k+neq] = 1;
    J[plya->index_Py_minus+j+neq][plya->index_Ps_minus+j+neq] = 1;
    J[plya->index_Py_minus+j+neq][plya->index_Px_plus+j+neq] = 1;
    J[plya->index_Py_minus+j+neq][plya->index_Px_minus+j+neq] = 1;
    J[plya->index_Py_minus+j+neq][plya->index_Py_minus+j+neq] = 1;
  }

  //Store pattern in sparse column compressed form:
  plya->Ap[0] = 0;
  nz = 0;
  for (i=0; i<2*neq; i++){
    for (j=0; j<2*neq; j++){
      if (J[j][i] == 1){
	plya->Ai[nz] = j;
	nz++;
      }
    }
    plya->Ap[i+1] = nz;
  }
  plya->Ai = realloc(plya->Ai,sizeof(int)*nz);

  free(J[0]);
  free(J);

  /*
    Setup the workspace for the lyapunov calculation.
  */

  return _SUCCESS_;
};

int free_lya_param(lya_param *plya){
  int i;
  free(plya->xi);
  free(plya->ui);
  free(plya->vi);
  free(plya->Tvec);
  free(plya->duidT);
  free(plya->dvidT);
  free(plya->duidx);
  free(plya->dxidT);
  free(plya->a);
  free(plya->y_0);
  free(plya->maxstep);
  free(plya->x_grid);
  free(plya->u_grid);
  free(plya->v_grid);
  free(plya->dvdu_grid);
  free(plya->dudT_grid);
  for (i=0; i<(plya->Nres+2); i++) 
    free(plya->mat[i]);
  free(plya->mat);
  free(plya->vv);
  free(plya->indx);
  free(plya->Ap);
  free(plya->Ai);
  background_free_dof(&(plya->pbs));
       
  return _SUCCESS_;
};


int lya_initial_conditions(double Ti, double *y, lya_param *plya){
  /** Set initial conditions at temperature Ti: */
  int i;
  ErrorMsg error_message;
  double x, L, n_plus=2.0;
  double Vx,D,Vz,Vz_bar,Px,Py,Px_bar,Py_bar,v_length_sq;

  //Assuming y is calloc'ed -- dangerous, better to zero it.
  for (i=0; i<plya->neq*2; i++) y[i] = 0.0;

  L = plya->L_initial;
  //Set standard equilibrium initial conditions:
  y[plya->index_L] = L/_L_SCALE_;
  for (i=0; i<plya->vres; i++)
    y[plya->index_Pa_plus+i] = 4.0;
  
/** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  if (plya->is_electron==_TRUE_)
    plya->g_alpha = 1.0+4.0/((1.0-_SIN2_THETA_W_)*n_plus);
  else
    plya->g_alpha = 1.0;
  plya->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(Ti,3)/_PI_/_PI_*L;
  plya->Vx = plya->delta_m2/(2.0*Ti)*sin(2.0*plya->theta_zero);
  plya->V0 = -plya->delta_m2/(2.0*Ti)*cos(2.0*plya->theta_zero);
  plya->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
    _G_F_/_M_Z_/_M_Z_*pow(Ti,5)*n_plus*plya->g_alpha;

  lya_get_resonances_xi(Ti,L,plya);  
  //Get new x_grid
  lasagna_call(lya_get_parametrisation(Ti,
				   plya, 
				   error_message),
	       error_message,error_message);
  

  for (i=0; i<plya->vres; i++){
    x = plya->x_grid[i];
    Vx = plya->Vx/x;
    Vz = plya->V0/x + plya->V1*x + plya->VL;
    Vz_bar = plya->V0/x + plya->V1*x - plya->VL;
    D = 0.5*plya->C_alpha*_G_F_*_G_F_*x*pow(Ti,5);

    Px = Vx*Vz/(D*D+Vz*Vz);
    Px_bar = Vx*Vz_bar/(D*D+Vz_bar*Vz_bar);
    Py = -Vx*D/(D*D+Vz*Vz);
    Py_bar = -Vx*D/(D*D+Vz_bar*Vz_bar);

    //continue;
    y[plya->index_Px_plus+i] =  Px + Px_bar;
    y[plya->index_Px_minus+i] = Px - Px_bar;
    y[plya->index_Py_plus+i] = Py + Py_bar;
    y[plya->index_Py_minus+i] = Py - Py_bar;
 
  }

  //Pick a random vector using the provided seed and normalize it.
  srand(plya->lyapunov_seed);
  v_length_sq = 0;
  for(i=plya->neq; i<2*plya->neq; i++){
    y[i] = rand();
    v_length_sq += y[i]*y[i];
  }
  for(i=plya->neq; i<2*plya->neq; i++)
    y[i] /= sqrt(v_length_sq)*plya->v_scale;
  /*
  for(i=plya->neq; i<2*plya->neq; i++)
    y[i] = 0;
  */

  //initial conditions for lya_stop_at_divL
  plya->max_old = L;
  plya->max_cur = L;
  plya->should_break = _FALSE_;
  plya->breakpoint = 0;
  return _SUCCESS_;
};

int lya_init_output(lya_param *plya){
  int Tres=plya->Tres;
  int vres=plya->vres;
  int Nres=plya->Nres;
  int handle;
  double tmp_array[3];
  int32_t tmp_array_int[3];
  char *outf=plya->output_filename;
  //Initialises the output file and stores parameters
  mat_create_file(outf);
  //Add matrices that are defined on momentum grid:
  mat_add_matrix(outf,"Pa_plus",miDOUBLE,Tres,vres,&(plya->Pa_plus_handle));
  mat_add_matrix(outf,"Pa_minus",miDOUBLE,Tres,vres,&(plya->Pa_minus_handle));
  mat_add_matrix(outf,"Ps_plus",miDOUBLE,Tres,vres,&(plya->Ps_plus_handle));
  mat_add_matrix(outf,"Ps_minus",miDOUBLE,Tres,vres,&(plya->Ps_minus_handle));
  mat_add_matrix(outf,"Px_plus",miDOUBLE,Tres,vres,&(plya->Px_plus_handle));
  mat_add_matrix(outf,"Px_minus",miDOUBLE,Tres,vres,&(plya->Px_minus_handle));
  mat_add_matrix(outf,"Py_plus",miDOUBLE,Tres,vres,&(plya->Py_plus_handle));
  mat_add_matrix(outf,"Py_minus",miDOUBLE,Tres,vres,&(plya->Py_minus_handle));
  mat_add_matrix(outf,"x_grid",miDOUBLE,Tres,vres,&(plya->x_grid_handle));
  mat_add_matrix(outf,"u_grid",miDOUBLE,Tres,vres,&(plya->u_grid_handle));
  mat_add_matrix(outf,"v_grid",miDOUBLE,Tres,vres,&(plya->v_grid_handle));
  //Add resonance dependent matrices:
  mat_add_matrix(outf,"xi_vec",miDOUBLE,Tres,Nres,&(plya->xi_handle));
  mat_add_matrix(outf,"ui_vec",miDOUBLE,Tres,Nres,&(plya->ui_handle));
  mat_add_matrix(outf,"vi_vec",miDOUBLE,Tres,Nres,&(plya->vi_handle));
  mat_add_matrix(outf,"b_a_vec",miDOUBLE,Tres,(1+Nres),&(plya->b_a_vec_handle));
  //Add other matrices:
  mat_add_matrix(outf,"L_vec",miDOUBLE,Tres,1,&(plya->L_handle));
  mat_add_matrix(outf,"T_vec",miDOUBLE,Tres,1,&(plya->T_handle));
  mat_add_matrix(outf,"I_vec",miDOUBLE,Tres,1,&(plya->I_handle));
  mat_add_matrix(outf,"I_conserved",miDOUBLE,Tres,1,&(plya->I_conserved_handle));
  mat_add_matrix(outf,"V0_vec",miDOUBLE,Tres,1,&(plya->V0_handle));
  mat_add_matrix(outf,"V1_vec",miDOUBLE,Tres,1,&(plya->V1_handle));
  mat_add_matrix(outf,"Vx_vec",miDOUBLE,Tres,1,&(plya->Vx_handle));
  mat_add_matrix(outf,"VL_vec",miDOUBLE,Tres,1,&(plya->VL_handle));
  //Add lyapunov vector
  mat_add_matrix(outf,"lya_vector",miDOUBLE,Tres,plya->neq,&(plya->v_handle));
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
  //mat_write_data(outf,"T",plya->Tvec,0,Tres);
  tmp_array[0] = plya->delta_m2; tmp_array[1] = plya->theta_zero;
  mat_write_data(outf,"delta_m2_theta_zero",&(tmp_array),0,2);
  mat_write_data(outf,"is_electron",&(plya->is_electron),0,1);
  tmp_array_int[0]=plya->Tres; tmp_array_int[1]=plya->vres;
  mat_write_data(outf,"Tres_vres",tmp_array_int,0,2);
  tmp_array[0] = plya->xmin; tmp_array[1] = plya->xext; 
  tmp_array[2] = plya->xmax;
  mat_write_data(outf,"xmin_xext_xmax",tmp_array,0,3);
  tmp_array[0] = plya->alpha; tmp_array[1] = plya->rs;
  mat_write_data(outf,"alpha_rs",tmp_array,0,2);
  return _SUCCESS_;
}


int lya_store_output(double T,
			    double *y,
			    double *dy,
			    int index_t,
			    void *param,
			    ErrorMsg error_message){
  lya_param *plya=param;
  int vres=plya->vres;
  int Nres=plya->Nres;
  int i;
  double x,xp1,f0,f0p1,Pa_minus,Ps_minus,Ps_minusp1,Pa_minusp1,I_PaPs,L,Ilost,v_length_sq; 
  FILE *mat_file;
  /** Calculate integrated quantities for convenience: */
 
  I_PaPs = 0.0;
  for (i=0; i<plya->vres-1; i++){
    x = plya->x_grid[i];
    xp1 = plya->x_grid[i+1];
    f0 = 1.0/(1.0+exp(x));
    f0p1 = 1.0/(1.0+exp(xp1));
    Ps_minus = y[plya->index_Py_minus+i];
    Ps_minusp1 = y[plya->index_Py_minus+i+1];
    Pa_minus = y[plya->index_Pa_plus+i];
    Pa_minusp1 = y[plya->index_Pa_plus+i+1];
    I_PaPs += 0.5*(xp1-x)*(x*x*f0*(Ps_minus+Pa_minus)+
			   xp1*xp1*f0p1*(Ps_minusp1+Pa_minusp1));
  }
  //Calculate the information lost:
  v_length_sq = 0;
  for(i=plya->neq; i<2*plya->neq; i++)
    v_length_sq += y[i]*y[i];
  Ilost = log(sqrt(v_length_sq)*plya->v_scale)/log(2);

  //printf("Storing output at index: %d\n",index_t);
  mat_file = fopen(plya->output_filename,"r+b");
  //Write stuff from y_vector:
  L = y[plya->index_L]*_L_SCALE_;
  mat_write_fast(&L,plya->L_handle,8,1);
  mat_write_fast(y+plya->index_Pa_plus,plya->Pa_plus_handle,8,vres);
  mat_write_fast(y+plya->index_Pa_minus,plya->Pa_minus_handle,8,vres);
  mat_write_fast(y+plya->index_Ps_plus,plya->Ps_plus_handle,8,vres);
  mat_write_fast(y+plya->index_Ps_minus,plya->Ps_minus_handle,8,vres);
  mat_write_fast(y+plya->index_Px_plus,plya->Px_plus_handle,8,vres);
  mat_write_fast(y+plya->index_Px_minus,plya->Px_minus_handle,8,vres);
  mat_write_fast(y+plya->index_Py_plus,plya->Py_plus_handle,8,vres);
  mat_write_fast(y+plya->index_Py_minus,plya->Py_minus_handle,8,vres);
  //Write stuff from structure:
  mat_write_fast(plya->x_grid,plya->x_grid_handle,8,vres);
  mat_write_fast(plya->u_grid,plya->u_grid_handle,8,vres);
  mat_write_fast(plya->v_grid,plya->v_grid_handle,8,vres);
  mat_write_fast(plya->xi,plya->xi_handle,8,Nres);
  mat_write_fast(plya->ui,plya->ui_handle,8,Nres);
  mat_write_fast(plya->vi,plya->vi_handle,8,Nres);
  mat_write_fast(&(I_PaPs),plya->I_conserved_handle,8,1);
  mat_write_fast(&(plya->V0),plya->V0_handle,8,1);
  mat_write_fast(&(plya->V1),plya->V1_handle,8,1);
  mat_write_fast(&(plya->Vx),plya->Vx_handle,8,1);
  mat_write_fast(&(plya->VL),plya->VL_handle,8,1);
  mat_write_fast(&(plya->b),plya->b_a_vec_handle,8,1);
  mat_write_fast(&Ilost,plya->I_handle,8,1);
  mat_write_fast(y+plya->index_v,plya->v_handle,8,plya->neq);
  mat_write_fast(plya->vi,plya->b_a_vec_handle,8,Nres);//Wrong

  //Write temperature at last so we know that everything has been written:
  mat_write_fast(&T,plya->T_handle,8,1);
  fclose(mat_file);
  return _SUCCESS_;
};

int lya_print_variables(double T,
			double *y,
			double *dy,
			void *param,
			ErrorMsg error_message){
  
  lya_param *plya=param;
  int idx=53,i;
  double x,Vx,V0,V1,VL;
  double D,Gamma;
  double Pa_plus, Pa_minus, Ps_plus, Ps_minus, Px_plus, Px_minus;
  double Py_plus, Py_minus;
  double Ilost, v_length_sq;
  x = plya->x_grid[idx];
  Vx = plya->Vx/x;
  V0 = plya->V0/x;
  V1 = plya->V1*x;
  VL = plya->VL;
  Gamma = plya->C_alpha*_G_F_*_G_F_*x*pow(T,5);
  D = 0.5*Gamma;
  
  Pa_plus = y[plya->index_Pa_plus+idx];
  Pa_minus = y[plya->index_Pa_minus+idx];
  Ps_plus = y[plya->index_Ps_plus+idx];
  Ps_minus = y[plya->index_Ps_minus+idx];
  Px_plus = y[plya->index_Px_plus+idx];
  Px_minus = y[plya->index_Px_minus+idx];
  Py_plus = y[plya->index_Py_plus+idx];
  Py_minus = y[plya->index_Py_minus+idx];

  //Calculate the information lost:
  v_length_sq = 0;
  for(i=plya->neq; i<2*plya->neq; i++)
    v_length_sq += y[i]*y[i];
  Ilost = log(sqrt(v_length_sq)*plya->v_scale)/log(2);


  fprintf(stderr,
	  "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	  1e3*T,
	  y[plya->index_L]*_L_SCALE_,
	  Pa_plus,
	  Pa_minus,
	  Ps_plus,
	  Ps_minus,
	  Px_plus,
	  Px_minus,
	  Py_plus,
	  Py_minus,
	  dy[plya->index_L]/y[plya->index_L],
	  VL,
	  -(V0+V1)*Px_plus-VL*Px_minus, //15
	  0.5*Vx*(Pa_plus-Ps_plus),
	  D*Py_plus,
	  -(V0+V1)*Px_minus-VL*Px_plus,
	  0.5*Vx*(Pa_minus-Ps_minus),
	  D*Py_minus,
	  Ilost);

//lya_param *plya=param;
  
  return _SUCCESS_;
};

int lya_derivs(double T, 
	       double *y, 
	       double *dy, 
	       void *param,
	       ErrorMsg error_message){
  lya_param *plya=param;
  double L;
  int i, j;
  double *v_grid=plya->v_grid;
  double *dvdu_grid=plya->dvdu_grid;
  double *dudT_grid=plya->dudT_grid;
  double gentr,H;
  double Vx, VL;
  double x;
  double Gamma, D, V0, V1, Pa_plus, Pa_minus, Ps_plus, Ps_minus;
  double Px_plus, Px_minus, Py_plus, Py_minus, f0,  mu_div_T;
  double feq_plus, feq_minus, I_VxPy_minus, I_f0Pa_plus, I_rho_ss;
  double dudTdvdu, delta_v;
  double rs;
  int idx, stencil_method;
  double n_plus = 2.0;
  double dLdT;
  int nnz, i_Ap;
  SCCformat *J_SCC;
  DNRformat *J_DNR;
  if (plya->is_electron==_TRUE_)
    plya->g_alpha = 1.0+4.0/((1.0-_SIN2_THETA_W_)*n_plus);
  else
    plya->g_alpha = 1.0;
  L = y[plya->index_L]*_L_SCALE_;
  
  /** Calculate 'scalar' potentials Vx, V0, VL
      (not momentum dependent): */
  plya->VL = sqrt(2.0)*_G_F_*2.0*_ZETA3_*pow(T,3)/_PI_/_PI_*L;
  VL = plya->VL;
  plya->Vx = plya->delta_m2/(2.0*T)*sin(2.0*plya->theta_zero);
  plya->V0 = -plya->delta_m2/(2.0*T)*cos(2.0*plya->theta_zero);
  plya->V1 = -7.0*_PI_*_PI_/(45.0*sqrt(2.0))*
    _G_F_/_M_Z_/_M_Z_*pow(T,5)*n_plus*plya->g_alpha;

  lya_get_resonances_xi(T,L,plya);
  
  //Get new x_grid
  lasagna_call(lya_get_parametrisation(T,
				   plya, 
				   error_message),
	       error_message,error_message);
  
  //Get integrated quantities
  lasagna_call(lya_get_integrated_quantities(y,
					 plya,
					 &I_VxPy_minus,
					 &I_f0Pa_plus,
					 &I_rho_ss,
					 error_message),
	       error_message,error_message);
   
  //Get degrees of freedom from background:
  background_getdof(T,NULL,&gentr,&(plya->pbs));
  /** Use Friedmann equation in radiation dominated universe: 
      (The radiation approximation breaks down long before there 
      is a difference in g and gS) */  
  H = sqrt(8.0*pow(_PI_,3)*gentr/90.0)*T*T/_M_PL_;

  dLdT = -1.0/(8.0*H*T*_ZETA3_)*I_VxPy_minus;
  
  //Get partial derivatives:
  lasagna_call(lya_get_partial_derivatives(T,
				       L,
				       dLdT,
				       plya, 
				       error_message),
	       error_message,error_message);

  /** Calculate RHS: */
  dy[plya->index_L] = dLdT/_L_SCALE_;
  //Solving mu from L, using the chebyshev cubic root:
  mu_div_T = -2*_PI_/sqrt(3.0)*
    sinh(1.0/3.0*asinh(-18.0*sqrt(3.0)*_ZETA3_*L/pow(_PI_,3)));
  
  /** All quantities defined on the grid: */
  delta_v = v_grid[1]-v_grid[0];
  for (i=0; i<plya->vres; i++){
    x = plya->x_grid[i];
    Vx = plya->Vx/x;
    V0 = plya->V0/x;
    V1 = plya->V1*x;
    Gamma = plya->C_alpha*_G_F_*_G_F_*x*pow(T,5);
    D = 0.5*Gamma;
    
    Pa_plus = y[plya->index_Pa_plus+i];
    Pa_minus = y[plya->index_Pa_minus+i];
    Ps_plus = y[plya->index_Ps_plus+i];
    Ps_minus = y[plya->index_Ps_minus+i];
    Px_plus = y[plya->index_Px_plus+i];
    Px_minus = y[plya->index_Px_minus+i];
    Py_plus = y[plya->index_Py_plus+i];
    Py_minus = y[plya->index_Py_minus+i];
  
    //Regulator for sterile population:
    rs = plya->rs;
    //Distributions:
    feq_plus = 1.0/(1.0+exp(x-mu_div_T))+1.0/(1.0+exp(x+mu_div_T));
    //Use an expansion for feq_minus since mu_div_T is very small.
    feq_minus = exp(x)*2*mu_div_T/(1+exp(x-mu_div_T))/(1+exp(x+mu_div_T));
    //Use the unexpanded expression if the error grows too large.
    if(exp(x+mu_div_T)/6*pow(mu_div_T,3)/(1+exp(x-mu_div_T))/
       (1+exp(x+mu_div_T))/feq_minus > 1e-10)
      feq_minus = 1.0/(1.0+exp(x-mu_div_T))-1.0/(1.0+exp(x+mu_div_T));

    f0 = 1.0/(1.0+exp(x));

    dudTdvdu = dudT_grid[i]*dvdu_grid[i];


    //Define index steps for calculating derivatives:
    if (i==0)
      stencil_method = 12; //First order, forward   12
    else if (i==plya->vres-1)
      stencil_method = 10; //First order, backwards 10
    else if ((i==plya->vres-2)||(i==1))
      stencil_method = 21; //Second order, centered 21
    else
      stencil_method = 21; //Fifth order, centered  51
    

    idx = plya->index_Pa_plus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_plus+Gamma*(2.0*feq_plus/f0-Pa_plus))+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = plya->index_Pa_minus+i;
    dy[idx] = -1.0/(H*T)*(Vx*Py_minus+Gamma*(2.0*feq_minus/f0-Pa_minus))+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = plya->index_Ps_plus+i;
    dy[idx] = 1.0/(H*T)*(Vx*Py_plus - rs*Gamma*(1.0/(6.0*_ZETA3_)*I_rho_ss*
						feq_plus-0.5*f0*Ps_plus)) +
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = plya->index_Ps_minus+i;
    dy[idx] = 1.0/(H*T)*(Vx*Py_minus-rs*Gamma*f0*
			 (1.0-0.5*(Pa_minus+Ps_minus))) +
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = plya->index_Px_plus+i;
    dy[idx] = 1.0/(H*T)*((V0+V1)*Py_plus+VL*Py_minus+D*Px_plus)+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = plya->index_Px_minus+i;
    dy[idx] = 1.0/(H*T)*((V0+V1)*Py_minus+VL*Py_plus+D*Px_minus)+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);

    idx = plya->index_Py_plus+i;
    dy[idx] = 1.0/(H*T)*
      (-(V0+V1)*Px_plus-VL*Px_minus+0.5*Vx*(Pa_plus-Ps_plus)+D*Py_plus)+ 
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);
   
    idx = plya->index_Py_minus+i;
    dy[idx] = 1.0/(H*T)*
      (-(V0+V1)*Px_minus-VL*Px_plus+0.5*Vx*(Pa_minus-Ps_minus)+D*Py_minus)+
      dudTdvdu*drhodv(y, delta_v, idx, stencil_method);
  }

  
  /*
    Calculate dv/dT = J*v.
  */

  //Scale the L-component of v:
  y[plya->neq] /= _L_SCALE_;
  switch(plya->LinearAlgebraWrapper){
  case LINALG_WRAPPER_SPARSE:
  case LINALG_WRAPPER_SUPERLU:

    J_SCC = (SCCformat *) ((**(plya->J_pp)).Store);

    //Calculate J*v.

    for(i=0; i<plya->neq; i++) dy[i+plya->neq] = 0;
    for(i=0; i<plya->neq; i++){
      for(j=J_SCC->Ap[i]; j<J_SCC->Ap[i+1]; j++){
	if(J_SCC->Ai[j] < plya->neq)
	  dy[J_SCC->Ai[j]+plya->neq] += ((double*) J_SCC->Ax)[j]*y[i+plya->neq];
      }
    }
    break;

  case LINALG_WRAPPER_DENSE_NR:
    J_DNR = (DNRformat *) ((**(plya->J_pp)).Store);
    for(i=0; i<plya->neq; i++){
      dy[i+plya->neq]=0.0;
      for(j=0; j<plya->neq; j++){
	dy[i+plya->neq] += ((double **) J_DNR->Matrix)[i+1][j+1]*y[j+plya->neq];
      }
    }
    break;

  }

  // Scale the component of v that comes from L.
  y[plya->neq] *= _L_SCALE_;
  dy[plya->neq] *= _L_SCALE_;
  /*
  for(i=0; i<plya->neq; i++) printf("%6.0e",y[i]);
  printf("\n");
  for(i=0; i<plya->neq; i++) printf("%6.0e",dy[i]);
  printf("\n");
  for(i=0; i<plya->neq; i++) printf("%6.0e",y[i+plya->neq]);
  printf("\n");
  for(i=0; i<plya->neq; i++) printf("%6.0e",dy[i+plya->neq]);
  printf("\n");
  PrintMultiMatrix(*(plya->J_pp),"Jacobian");
  getchar();
  */
  return _SUCCESS_;
}


int lya_get_resonances_xi(double T, 
		      double L,
		      lya_param *plya){
  double *xi=plya->xi;
  double x0, x0b,A;
  int i;

  x0 = sqrt(fabs(plya->V0/plya->V1));
  
  xi[0] = x0;
  xi[1] = x0;
  A = fabs(0.5*plya->VL/sqrt(fabs(plya->V0*plya->V1)));
  if (plya->delta_m2<0.0){
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
  for (i=0; i<plya->Nres; i++){
    if(xi[i]>plya->xmax){
      xi[i] = plya->xmax;
    }
  }
  return _SUCCESS_;
};

int lya_get_resonances_dxidT(double T, 
			 double L,
			 double dLdT,
			 lya_param *plya){
  double x0, A,A2;
  double *dxidT=plya->dxidT;
  double F;
  double one_plus_dlogLdlogT;
  int i;

  x0 = sqrt(fabs(plya->V0/plya->V1));  
 
  dxidT[0] = -3.0*x0/T;
  dxidT[1] = -3.0*x0/T;

  if (fabs(L)>1e-100){
 
    A = fabs(0.5*plya->VL/sqrt(fabs(plya->V0*plya->V1)));
    A2 = A*A;
    one_plus_dlogLdlogT = 1.0+T/L*dLdT;
 
    if (plya->delta_m2<0.0){
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
  for (i=0; i<plya->Nres; i++){
    if(plya->xi[i]>=(plya->xmax-1e-12)){
      dxidT[i] = 0.0;
    }
  }
  return _SUCCESS_;
};


int lya_u_of_x(double x, double *u, double *dudx, lya_param *plya){
  double K,xmin,xmax,xext;
  xmin = plya->xmin;
  xmax = plya->xmax;
  xext = plya->xext;
  K = (xext+xmax)/(xmax-xmin);
  *u = K*(x-xmin)/(x+xext);
  *dudx = K*(xext+xmin)/pow(x+xext,2);

  return _SUCCESS_;
};

int lya_x_of_u(double u, double *x, lya_param *param){
  *x = param->xext*(u+param->eps1)/(1+param->eps2-u);
  return _SUCCESS_;
};

int lya_nonlinear_rhs(double *y, double *Fy, void *param){
  lya_param *plya=param;
  double b=y[0];
  double *vi=y+1;
  double alpha = plya->alpha;
  double *ui=plya->ui;
  double *a=plya->a;
  int i,n=plya->Nres;

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

int lya_nonlinear_rhs_jac(double *y, double **jac, void *param){
  lya_param *plya=param;
  double b=y[0];
  double *vi=y+1;
  double alpha = plya->alpha;
  int i,j,n=plya->Nres;
  
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

int lya_stop_at_divL(double t,
		     double *y,
		     double *dy,
		     void *param,
		     ErrorMsg error_message){
  lya_param *plya=param;
  double L, T, stopfactor,v_length_sq;
  int i;
  stopfactor = 2;
  L = y[plya->index_L]*_L_SCALE_;
  T = t;
  // Check for sign change:
  if(L*plya->max_cur < 0){
    plya->max_old = max(plya->max_old,fabs(plya->max_cur));
    plya->max_cur = 0;
    plya->should_break = _FALSE_;
  }
  // Is the value larger than the current max:
  if((fabs(L)>fabs(plya->max_cur)) && (plya->should_break==_FALSE_)){
    plya->max_cur = L;
    if((fabs(plya->max_cur) > plya->max_old*stopfactor)){
      plya->breakpoint = T - plya->T_wait;
      plya->should_break = _TRUE_;
    }
  }
  else if(plya->should_break==_TRUE_){
    if(fabs(L) < plya->max_old)
      plya->should_break = _FALSE_;
    else if(T<plya->breakpoint)
      return _TRUE_;
  }
  v_length_sq = 0;
  for(i=plya->neq; i<2*plya->neq; i++)
    v_length_sq += y[i]*y[i];
  if(sqrt(v_length_sq)*plya->v_scale > pow(2,plya->I_stop)){
    printf("The information loss is larger than %g bits. Calculation stopped.\n",plya->I_stop);
    return _TRUE_;
  }
  return _FALSE_;
};


int lya_print_L(double T,
			double *y,
			double *dy,
			void *param,
			ErrorMsg error_message){

  lya_param *plya=param;
  fprintf(stderr,"%.16e %.16e\n",1e3*T,y[plya->index_L]*_L_SCALE_);
  return _SUCCESS_;
}

int lya_get_integrated_quantities(double *y,
			      lya_param *plya,
			      double *I_VxPy_minus,
			      double *I_f0Pa_plus,
			      double *I_rho_ss,
			      ErrorMsg error_message){
  int i;
  double w_trapz, x, x2, f0, Vx;
  double Py_minus, Pa_plus, Ps_plus_Ps_minus;
			      
  /** Integrated quantities needed. We integrate in x space: */
  *I_VxPy_minus = 0.0;
  *I_f0Pa_plus = 0.0;
  *I_rho_ss = 0.0;
  for (i=0; i<plya->vres; i++){
    if (i==0)
      w_trapz = 0.5*(plya->x_grid[i+1]-plya->x_grid[i]);
    else if (i==plya->vres-1)
      w_trapz = 0.5*(plya->x_grid[i]-plya->x_grid[i-1]);
    else
      w_trapz = 0.5*(plya->x_grid[i+1]-plya->x_grid[i-1]);
    x = plya->x_grid[i];
    x2 = x*x;
    f0 = 1.0/(1.0+exp(x));
    Vx = plya->Vx/x;
    Py_minus = y[plya->index_Py_minus+i];
    Pa_plus = y[plya->index_Pa_plus+i];
    Ps_plus_Ps_minus = (y[plya->index_Ps_plus+i] + y[plya->index_Ps_minus+i]);
      
    *I_VxPy_minus += w_trapz*(x2*f0*Vx*Py_minus);
    *I_f0Pa_plus += w_trapz*(x2*f0*Pa_plus);
    *I_rho_ss += w_trapz*(x2*f0*Ps_plus_Ps_minus); //From equation 2.18 in KS01.
  }
  return _SUCCESS_;
}

int lya_get_parametrisation(double T,lya_param *plya, ErrorMsg error_message){
  double alpha=plya->alpha;
  double *maxstep=plya->maxstep;;
  double *y_0=plya->y_0;
  double *ui=plya->ui;
  double *vi=plya->vi;
  double *x_grid=plya->x_grid;
  double *u_grid=plya->u_grid;
  double *v_grid=plya->v_grid;
  double tol_newton=1e-12;
  double wi;
  int i,j;
  int niter;

  for (i=0; i<plya->Nres; i++){
    lya_u_of_x(plya->xi[i],ui+i,plya->duidx+i,plya);
  }
  
  //Establish guess and set maximum steps for Newton method:
  /** This has been a very hard 'bug' to trace: Some of the evolvers rely
      on the derivative function to give exactly the same output for the same
      inputs, so making a non-deterministic guess here spoils this!
  */
  if (1==1){//(plya->guess_exists == _FALSE_){
    y_0[0] = 1.0 - alpha;
    maxstep[0] = 100.0;
    for (i=1; i<=plya->Nres; i++){
      y_0[i] = ui[i-1];
      maxstep[i] = 0.1;
    }
  }
  else{
    //We can make a more realistic guess:
    y_0[0] = plya->b;
    maxstep[0] = 100.0;
    for (i=1; i<=plya->Nres; i++){
      y_0[i] = vi[i-1];
      maxstep[i] = 0.1;
    }
  }
  //Find parametrisation parameters vi, a and b using Newton:
  lasagna_call(Newton(lya_nonlinear_rhs,
		      lya_nonlinear_rhs_jac,
		      y_0,
		      plya,
		      maxstep,
		      tol_newton,
		      &niter,
		      100,
		      plya->Nres+1,
		      error_message),
	       error_message,error_message);
  plya->guess_exists = _TRUE_;
  plya->b = y_0[0];
  for (i=0; i<plya->Nres; i++){
    vi[i] = y_0[i+1];
    plya->a[i] = ui[i]-alpha*vi[i];
  }

  //Now update grids:
  /** Find the splitting of v and store it temporarily in plya->indx.
      We use the fact that v is uniform.
  */
  plya->indx[0] = 0;
  plya->indx[plya->Nres] = plya->vres;
  for (i=1; i<plya->Nres; i++){
    wi = 0.5*(vi[i-1]+vi[i]); //Weighted average
    plya->indx[i] = (int)((wi-plya->v_left)/(v_grid[1]-v_grid[0]));
    if (plya->indx[i]<plya->indx[i-1])
      plya->indx[i] = plya->indx[i-1];
  }
  // Now loop over resonances:
  for (i=0; i<plya->Nres; i++){
    if (plya->indx[i]==plya->indx[i+1])
      continue;
    //Loop over each segment:
    for (j=plya->indx[i]; j<plya->indx[i+1]; j++){
      u_grid[j] = alpha*v_grid[j]+plya->a[i]+plya->b*pow(v_grid[j]-vi[i],3);
      lya_x_of_u(u_grid[j],&(x_grid[j]),plya);
    }
  }

  return _SUCCESS_;
}

int lya_get_partial_derivatives(double T,
			    double L,
			    double dLdT,
			    lya_param *plya, 
			    ErrorMsg error_message){
  int i,j;
  double u1, v1, vN;
  double gamma_j, beta_j, wi;
  double alpha=plya->alpha;
  double *dvidT=plya->dvidT;
  double *duidT=plya->duidT;
  double *v_grid=plya->v_grid;
  double *ui=plya->ui;
  double *vi=plya->vi;
  double *dvdu_grid=plya->dvdu_grid;
  double *dudT_grid=plya->dudT_grid;
  double dbdT,daidT;
  double lu_sgn;

/** Get partial derivatives of vi with respect to T by solving a 
      linear system A*dvidT = B(duidT):
  */
  lya_get_resonances_dxidT(T,L,dLdT,plya);

  //Set partial derivatives of ui with respect to T:
  for (i=0; i<plya->Nres; i++){ 
    plya->duidT[i]=plya->duidx[i]*plya->dxidT[i];
  }
  for(j=0; j<plya->Nres; j++){
    for(i=0; i<plya->Nres; i++){
        plya->mat[j+1][i+1] = 0.0;
    }
  }
  u1 = ui[0];
  v1 = vi[0];
  vN = vi[plya->Nres-1];
  
  for (j=0; j<plya->Nres-1; j++){
    /** Enter non-zero entries in A: */
    gamma_j = 0.25*pow(vi[j+1]/v1-vi[j]/v1,3)*(2.0*alpha-3.0*u1/v1);
    beta_j = alpha + 0.75*plya->b*pow(vi[j+1]-vi[j],2);
    plya->mat[j+1][1] = gamma_j;
    plya->mat[j+1][j+1] -= beta_j;
    plya->mat[j+1][j+2] += beta_j;
    /** Setup RHS: */
    dvidT[j] = duidT[j+1]-duidT[j]-0.25*pow(vi[j+1]/v1-vi[j]/v1,3)*duidT[0];
  }
  plya->mat[plya->Nres][1] = -pow((1.0-vN)/v1,3)*(2.0*alpha-3.0*u1/v1);
  plya->mat[plya->Nres][plya->Nres] = alpha+3.0*plya->b*pow(1.0-vN,2);
  dvidT[plya->Nres-1] = duidT[plya->Nres-1]+pow((1.0-vN)/v1,3)*duidT[0];

  //LU decomposition of matrix:
  lasagna_call(ludcmp(plya->mat,plya->Nres,plya->indx,&lu_sgn,plya->vv),
	       error_message,error_message);
  lasagna_call(lubksb(plya->mat,plya->Nres,plya->indx,dvidT-1),
	       error_message,error_message);

  /** Calculate x and u on the v-grid along with derivatives dvdu and dudT:*/
  dbdT = (duidT[0]+(2.0*alpha-3.0*u1/v1)*dvidT[0])/pow(v1,3);
  /** Find the splitting of v and store it temporarily in plya->indx.
      We use the fact that v is uniform.
  */
  plya->indx[0] = 0;
  plya->indx[plya->Nres] = plya->vres;
  for (i=1; i<plya->Nres; i++){
    wi = 0.5*(vi[i-1]+vi[i]); //Weighted average
    plya->indx[i] = (int)((wi-plya->v_left)/(v_grid[1]-v_grid[0]));
    if (plya->indx[i]<plya->indx[i-1])
      plya->indx[i] = plya->indx[i-1];
  }

  // Now loop over resonances:
  for (i=0; i<plya->Nres; i++){
    if (plya->indx[i]==plya->indx[i+1])
      continue;
    //Loop over each segment:
    for (j=plya->indx[i]; j<plya->indx[i+1]; j++){
      daidT = duidT[i]-alpha*dvidT[i];
      dvdu_grid[j] = 1.0/(alpha+3.0*plya->b*pow(v_grid[j]-vi[i],2));
      dudT_grid[j] = daidT+dbdT*pow(v_grid[j]-vi[i],3)-
	3.0*plya->b*pow(v_grid[j]-vi[i],2)*dvidT[i];
    }
  }
  /** Debug area*/ 
  /**
  fprintf(plya->tmp,"%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	  T,
	  plya->a[0],
	  plya->a[1],
	  duidT[0]-alpha*dvidT[0],
	  duidT[1]-alpha*dvidT[1], 
	  plya->b,
	  dbdT,
	  plya->xi[0],
	  plya->xi[1],
	  plya->dxidT[0], 
	  plya->dxidT[1],
	  plya->ui[0],
	  plya->ui[1],
	  plya->duidT[0],
	  plya->duidT[1],
	  plya->vi[0],
	  plya->vi[1],
	  plya->dvidT[0],
	  plya->dvidT[1]);
  */
return _SUCCESS_;
}
