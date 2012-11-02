/** @file lasagna_lya.c 
 * Thomas Tram Bülow, 25.08.2011    
 */
 
#include "lasagna_lya.h"
#include <time.h>
int main(int argc, char **argv) {
  lya_param lya_struct;
  double *y_inout;
  int *interp_idx;
  ErrorMsg error_message;
  int i;
  int func_return;
  clock_t start, end;
  double cpu_time_used, elapsed;
  time_t wtime1, wtime2;

  EvolverOptions options;
  extern int evolver_radau5();
  extern int evolver_ndf15(); 	
  extern int evolver_rkdp45(); 	
  int (*generic_evolver)();  

  
  func_return = lya_input_init_from_arguments(argc, 
					      argv,
					      &lya_struct,
					      error_message);
  if (func_return == _FAILURE_){
    printf("\n\nError running input_init_from_arguments\n=>%s\n",error_message);
    return _FAILURE_;
  }

  //Dump jacobian pattern:
  FILE *jacfile=fopen("jac_anal_Ap.dat","w");
  for (i=0; i<=lya_struct.neq; i++) fprintf(jacfile,"%d ",lya_struct.Ap[i]);
  fclose(jacfile);
  jacfile=fopen("jac_anal_Ai.dat","w");
  for (i=0; i<lya_struct.Ap[lya_struct.neq]; i++) fprintf(jacfile,"%d ",lya_struct.Ai[i]);
  fclose(jacfile);

  //  return 0;

  /** Do stuff */
  y_inout = calloc(2*lya_struct.neq,sizeof(double));
  interp_idx = malloc(sizeof(int)*2*lya_struct.neq);
  for(i=0; i<2*lya_struct.neq; i++)
    interp_idx[i] = 1;

  if(lya_struct.evolver == 0){
    generic_evolver = evolver_radau5;
  }
  else if (lya_struct.evolver == 2){
    printf("Runge-Kutta evolver\n");
    generic_evolver = evolver_rkdp45;
  }
  else{
    generic_evolver = evolver_ndf15;
  }
  
  lya_initial_conditions(lya_struct.T_initial, 
			    y_inout, 
			    &lya_struct);


  lya_init_output(&lya_struct);

  //Handle options:
  DefaultEvolverOptions(&options,lya_struct.LinearAlgebraWrapper);
  options.used_in_output=interp_idx;
  options.RelTol = lya_struct.rtol;
  options.AbsTol = lya_struct.abstol;
  options.t_vec = lya_struct.Tvec; 
  options.tres = lya_struct.Tres;
  options.Ap = lya_struct.Ap;
  options.Ai = lya_struct.Ai;
  options.output = lya_store_output;
  //  options.print_variables = qke_print_L;
  //  options.stop_function = qke_stop_at_L;
  if(lya_struct.T_wait >= 0)  options.stop_function = qke_stop_at_divL;
  options.EvolverVerbose=lya_struct.verbose;
  options.Cores = lya_struct.nproc;
  options.J_pointer_flag = _TRUE_;
  lya_struct.J_pp = &(options.J_pointer);

  printf("theta: %g\n",lya_struct.theta_zero);
  start = clock();  
  time(&wtime1);
  func_return = generic_evolver(lya_derivs,
				&lya_struct,
				lya_struct.T_initial,
				lya_struct.T_final,
				y_inout, 
				2*lya_struct.neq, 
				&(options),
				error_message);

  end = clock();
  time(&wtime2);  
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("CPU time used: %g minutes.\n",cpu_time_used/60);
  elapsed = difftime(wtime2,wtime1);
  printf("Wall clock time used: %g minutes.\n",elapsed/60);
      
  
  free(y_inout);
  free(interp_idx);
  free_lya_param(&lya_struct);


  if (func_return == _FAILURE_){
    printf("%s\n",error_message);
    return _FAILURE_;
  }
  else{
    return _SUCCESS_;
  }
}
