/** @file lasagna.c 
 * Thomas Tram Bülow, 25.08.2011    
 */
 
#include "lasagna.h"
#include <time.h>
int main(int argc, char **argv) {
  qke_param qke_struct;
  double *y_inout;
  int *interp_idx;
  ErrorMsg error_message;
  int i;
  int func_return;
  clock_t start, end;
  double cpu_time_used;

  extern int evolver_radau5();
  extern int evolver_ndf15(); 	
  extern int evolver_rkdp45(); 	
  int (*generic_evolver)();  

  func_return = input_init_from_arguments(argc, 
					  argv,
					  &qke_struct,
					  error_message);
  if (func_return == _FAILURE_){
    printf("\n\nError running input_init_from_arguments\n=>%s\n",error_message);
    return _FAILURE_;
  }

  /** Do stuff */
  y_inout = calloc(qke_struct.neq,sizeof(double));
  interp_idx = malloc(sizeof(int)*qke_struct.neq);
  for(i=0; i<qke_struct.neq; i++)
    interp_idx[i] = 1;

  if(qke_struct.evolver == 0){
    generic_evolver = evolver_radau5;
  }
  else if (qke_struct.evolver == 2){
    printf("Runge-Kutta evolver\n");
    generic_evolver = evolver_rkdp45;
  }
  else{
    generic_evolver = evolver_ndf15;
  }
  
  qke_initial_conditions(qke_struct.T_initial, 
			 y_inout, 
			 &qke_struct);


  qke_init_output(&qke_struct);

  printf("theta: %g\n",qke_struct.theta_zero);
  start = clock();  
    func_return = generic_evolver(qke_derivs,
				  qke_struct.T_initial,
				  qke_struct.T_final,
				  y_inout, 
				  interp_idx,
				  qke_struct.neq, 
				  &qke_struct,
				  qke_struct.rtol, 
				  qke_struct.abstol, 
				  qke_struct.Tvec, 
				  qke_struct.Tres,
				  qke_struct.Ap,
				  qke_struct.Ai,
				  qke_store_output,
				  qke_print_L,//NULL,
				  NULL,//qke_stop_at_L,
				  error_message);
  
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("CPU time used: %g minutes.\n",cpu_time_used/60.0);    
  
  free(y_inout);
  free(interp_idx);
  free_qke_param(&qke_struct);


  if (func_return == _FAILURE_){
    printf("%s\n",error_message);
    return _FAILURE_;
  }
  else{
    return _SUCCESS_;
  }
}
