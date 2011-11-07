/** @file lasagna.c 
 * Thomas Tram Bülow, 25.08.2011    
 */
 
#include "lasagna.h"

int main() {
  double tmp1,tmp2;
  struct background_structure bg_struct;
  struct qke_param qke_struct;
  double *y_inout;
  double Tmid;
  int *interp_idx;
  double Ti,Tf,Vl,Vr,sinsq2theta;
  ErrorMsg error_message;
  int i;
  int func_return;
  printf("Hello World! Garfield is coming.. \n");


  strncpy(bg_struct.dof_filename,"dsdofHP_B.dat",_FILENAMESIZE_);
  printf("%s",bg_struct.dof_filename);
  background_init_dof(&bg_struct);
  background_getdof(5.0e-3,&tmp1,&tmp2,&bg_struct);
  printf("From file: %s, found sqrtg=%g and gentr=%g\n",
	 bg_struct.dof_filename,tmp1,tmp2);
  qke_struct.pbs = &bg_struct;

  /** Parameters for the u(x) transformation must be known
      before calling init_qke_param.*/
  qke_struct.xext = 3.1;
  qke_struct.xmin =  0.0001;// 0.0001; //1e-4;
  qke_struct.xmax = 100.0; //100.0;

#ifdef _OPENMP
  qke_struct.nproc = omp_get_num_procs();
  printf("Number of cores: %d\n",qke_struct.nproc);
#endif

  qke_struct.evolve_vi = _FALSE_;
  init_qke_param(&qke_struct,2,256,10000); //vres Tres
  /** Do stuff */
  y_inout = calloc(qke_struct.neq,sizeof(double));
  interp_idx = malloc(sizeof(int)*qke_struct.neq);
  for(i=0; i<qke_struct.neq; i++)
    interp_idx[i] = 1;
  Ti = 0.035;
  Tf = 0.010;
  Tmid = 0.027;
  for(i=0; i<qke_struct.Tres/2; i++)
    qke_struct.Tvec[i] = Ti+i*(Tmid-Ti)/(qke_struct.Tres/2-1);
  for(i=0; i<qke_struct.Tres/2; i++)
    qke_struct.Tvec[i+qke_struct.Tres/2] =Tmid+(i+1)*(Tf-Tmid)/(qke_struct.Tres/2);
  Vr = 1.0;
  Vl = 0.0;
  for(i=0; i<qke_struct.vres; i++){
    qke_struct.v_grid[i] = Vl+i*(Vr-Vl)/(qke_struct.vres-1);
  }
  /** We must have non-zero alpha, otherwise the matrix for 
      solving for dvidT becomes singular.
  */
  strcpy(qke_struct.output_filename,"output/dump.mat");
  qke_struct.alpha = 0.3;
  qke_struct.rs = 0.0;
  qke_struct.L_initial = 2e-10;
  qke_struct.delta_m2 = -1e-17; //-10e-18;
  sinsq2theta = pow(10.0,-6.1);
  qke_struct.theta_zero = 0.5*asin(sqrt(sinsq2theta));
  qke_struct.is_electron = _FALSE_;
  if (qke_struct.is_electron == _TRUE_)
    qke_struct.C_alpha = 1.27;
  else
    qke_struct.C_alpha = 0.92;

  qke_initial_conditions(Ti, y_inout, &qke_struct);
  
  qke_init_output(&qke_struct);
  /**
     printf("First resonance at x = %g, last resonance at %g.\n",
	 sqrt(1.64e8*fabs(qke_struct.delta_m2)*cos(2.0*qke_struct.theta_zero)/pow(Ti,6)), sqrt(1.64e8*fabs(qke_struct.delta_m2)*cos(2.0*qke_struct.theta_zero)/pow(Tf,6)));
  return _SUCCESS_;
  */
  func_return = evolver_ndf15(qke_derivs,
			      Ti,
			      Tf,
			      y_inout, 
			      interp_idx,
			      qke_struct.neq, 
			      &qke_struct,
			      1e-2, 
			      1e-14, 
			      qke_struct.Tvec, 
			      qke_struct.Tres,
			      qke_store_output,
			      NULL,
			      error_message);
  
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
