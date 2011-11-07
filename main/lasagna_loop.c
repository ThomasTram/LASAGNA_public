/** @file lasagna_loop.c 
 * Thomas Tram Bülow, 25.08.2011    
 */
 
#include "lasagna_loop.h"

int main() {
  double tmp1,tmp2;
  struct background_structure bg_struct;
  struct qke_param qke_struct;
  double *y_inout;
  double Tmid;
  int *interp_idx;
  double Ti,Tf,Vl,Vr,sinsq2theta;
  ErrorMsg error_message;
  int i,j,k;
  int func_return;

  int nt=5, nm=5;
  double sinsq2theta_exp[nt];
  double deltamsq[nm];
  double sinsq2theta_exp_min = -8.0, sinsq2theta_exp_max=-5.0;
  double deltamsq_exp_min=-20.0, deltamsq_exp_max=-17.0;
  double sinsq2theta_exp2, deltamsq2;
  
  for (i=0; i<nt; i++)
    sinsq2theta_exp[i] = sinsq2theta_exp_min + 
      i*(sinsq2theta_exp_max-sinsq2theta_exp_min)/(nt-1.0);
  for (i=0; i<nm; i++)
    deltamsq[i] = pow(10.0,deltamsq_exp_min+i*deltamsq_exp_max/(nm-1.0));

  strncpy(bg_struct.dof_filename,"dsdofHP_B.dat",_FILENAMESIZE_);
  printf("%s",bg_struct.dof_filename);
  background_init_dof(&bg_struct);


#ifdef _OPENMP
  qke_struct.nproc = omp_get_num_procs();
  printf("Number of cores: %d\n",qke_struct.nproc);
  omp_set_num_threads(qke_struct.nproc);
#pragma omp parallel for default(shared) private(j,i,y_inout,Ti,Tf,Tmid,Vr,Vl,func_return,error_message,interp_idx,qke_struct) schedule(static)
#endif
  for (j=0; j<(nt*nm); j++){    
    sinsq2theta_exp2 = sinsq2theta_exp[j%5];
    deltamsq2 = deltamsq[j/5];
    qke_struct.theta_zero = 0.5*asin(pow(10.0,0.5*sinsq2theta_exp2));
    qke_struct.delta_m2 = deltamsq2; //-10e-18;    
    qke_struct.pbs = &bg_struct;

    /** Parameters for the u(x) transformation must be known
	before calling init_qke_param.*/
    qke_struct.xext = 3.1;
    qke_struct.xmin =  0.0001;// 0.0001; //1e-4;
    qke_struct.xmax = 100.0; //100.0;
    qke_struct.nproc = 1;
    qke_struct.evolve_vi = _FALSE_;
    init_qke_param(&qke_struct,2,256,100); //vres Tres
    /** Do stuff */
    y_inout = calloc(qke_struct.neq,sizeof(double));
  
    Ti = 0.070;
    Tf = 0.001;
    Tmid = 0.035;

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

    sprintf(qke_struct.output_filename,"output/dump_%g_%g.mat",sinsq2theta_exp2,deltamsq2);
    printf("%s\n",qke_struct.output_filename);
    qke_struct.alpha = 0.3;
    qke_struct.rs = 0.0;
    qke_struct.L_initial = 1e-10;

    qke_struct.is_electron = _FALSE_;
    if (qke_struct.is_electron == _TRUE_)
      qke_struct.C_alpha = 1.27;
    else
      qke_struct.C_alpha = 0.92;
    
    interp_idx = malloc(sizeof(int)*qke_struct.neq);
    for(i=0; i<qke_struct.neq; i++)
      interp_idx[i] = 1;

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
				1e-3, 
				1e-14, 
				qke_struct.Tvec, 
				qke_struct.Tres,
				qke_store_output,
				NULL,
				error_message);
    if (func_return == _FAILURE_)
      printf("Failure at j=%d, message: %s.\n",j,error_message);
    
    free(y_inout);
    free(interp_idx);
    free_qke_param(&qke_struct);  
  }
  background_free_dof(&bg_struct);
  return _SUCCESS_;
}
