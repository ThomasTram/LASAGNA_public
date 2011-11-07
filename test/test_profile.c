/** @file lasagna.c 
 * Thomas Tram Bülow, 25.08.2011    
 */
#include "test_profile.h"
int main() {
  double tmp1,tmp2;
  struct background_structure bg_struct;
  struct qke_param qke_struct;
  double *y_inout;
  double *dy;
  int *interp_idx;
  double Ti,Tf,Vl,Vr,sinsq2theta;
  ErrorMsg error_message;
  int i;
  FILE *datfile;
  time_t start,stop;
  double spend_time;

  printf("Hello World! Garfield is coming.. \n");


  strncpy(bg_struct.dof_filename,"dsdofHP_B.dat",_FILENAMESIZE_);
  printf("%s",bg_struct.dof_filename);
  background_init_dof(&bg_struct);
  background_getdof(5.0e-3,&tmp1,&tmp2,&bg_struct);
  printf("From file: %s, found sqrtg=%g and gentr=%g\n",
	 bg_struct.dof_filename,tmp1,tmp2);
  qke_struct.pbs = &bg_struct;

  qke_struct.nproc = 1; //omp_get_num_procs();
  printf("Number of cores: %d\n",qke_struct.nproc);
  qke_struct.xext = 3.1;
  qke_struct.xmin = 0.1; //1e-4;
  qke_struct.xmax = 20.0; //100.0;
  init_qke_param(&qke_struct,2,100,2560); //vres Tres
  /** Do stuff */
  y_inout = calloc(qke_struct.neq,sizeof(double));
  dy = malloc(sizeof(double)*qke_struct.neq);
  interp_idx = malloc(sizeof(int)*qke_struct.neq);
  for(i=0; i<qke_struct.neq; i++)
    interp_idx[i] = 1;
  Ti = 50e-3;
  Tf = 30e-3;
  
  for(i=0; i<qke_struct.Tres; i++)
    qke_struct.Tvec[i] = Ti+i*(Tf-Ti)/(qke_struct.Tres-1);
  Vr = 1.0;
  Vl = 0.0;
  for(i=0; i<qke_struct.vres; i++){
    qke_struct.v_grid[i] = Vl+i*(Vr-Vl)/(qke_struct.vres-1);
    printf("%g ",qke_struct.v_grid[i]);
  }
  printf("\n");
  /** We must have non-zero alpha, otherwise the matrix for 
      solving for dvidT becomes singular.
  */
  strcpy(qke_struct.output_filename,"output/dump.mat");
  qke_struct.alpha = 1.0;
  qke_struct.rs = 0.0;
  qke_struct.L_initial = 5e-10;
  qke_struct.delta_m2 = -1e-18; //-10e-18;
  sinsq2theta = 1e-7;
  qke_struct.theta_zero = 0.5*asin(sqrt(sinsq2theta));
  
  qke_struct.is_electron = _FALSE_;
  if (qke_struct.is_electron == _TRUE_)
    qke_struct.C_alpha = 1.27;
  else
    qke_struct.C_alpha = 0.92;

  qke_initial_conditions(Ti, y_inout, &qke_struct);

  start = time(NULL);
  for (i=0; i<20000; i++){
    qke_derivs(Ti, y_inout, dy, &qke_struct, error_message);
  }
  stop = time(NULL);
  printf("%ld seconds spend.\n",stop-start);
  int func_return = _SUCCESS_;
  free(dy);
    free(y_inout);
    free(interp_idx);
    free_qke_param(&qke_struct);
    background_free_dof(&bg_struct);

    if (func_return == _FAILURE_){
      printf("%s\n",error_message);
      return _FAILURE_;
    }
    else{
      return _SUCCESS_;
    }
}
