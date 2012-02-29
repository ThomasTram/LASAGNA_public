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
  double T_lowtemp;

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
  else if (qke_struct.evolver == 1){
    generic_evolver = evolver_ndf15;
  }
  else {
    generic_evolver = evolver_rkdp45;
  }
  
  qke_initial_conditions_fixed_grid(qke_struct.T_initial, 
				    y_inout, 
				    &qke_struct);
  
  qke_init_output(&qke_struct);

  T_lowtemp = 4e-3;
    //qke_struct.T_initial-0.1*(qke_struct.T_initial-qke_struct.T_final);
  start = clock();  
  func_return = generic_evolver(qke_derivs_fixed_grid,
				qke_struct.T_initial,
				T_lowtemp,
				y_inout, 
				interp_idx,
				qke_struct.neq, 
				&qke_struct,
				qke_struct.rtol, 
				qke_struct.abstol, 
				qke_struct.Tvec, 
				qke_struct.Tres,
				NULL,//qke_struct.Ap,
				NULL,//qke_struct.Ai,
				qke_store_output,
				qke_print_variables,
				NULL,//qke_stop_at_L,
				error_message);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("CPU time used: %g minutes.\n",cpu_time_used/60.0);    

  qke_struct.L_initial = 0.0;
  //y_inout[qke_struct.index_L] = 0.0;
  start = clock();  
  func_return = generic_evolver(qke_derivs_fixed_grid,
				T_lowtemp,
				qke_struct.T_final,
				y_inout, 
				interp_idx,
				qke_struct.neq, 
				&qke_struct,
				qke_struct.rtol, 
				qke_struct.abstol, 
				qke_struct.Tvec, 
				qke_struct.Tres,
				NULL,//qke_struct.Ap,
				NULL,//qke_struct.Ai,
				qke_store_output,
				qke_print_variables,
				NULL,//qke_stop_at_L,
				error_message);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("CPU time used: %g minutes.\n",cpu_time_used/60.0);    


  //Proceed using low-temperature basis:
  /**
  init_qke_param_fixed_grid_low_temp(&qke_struct);
  double *y_inout2=malloc(sizeof(double)*qke_struct.neq);
  qke_copy_IC_fixed_grid_to_fixed_grid_low_temp(T_lowtemp,
						y_inout,
						y_inout2,
						&qke_struct);
  free(y_inout);
  free(interp_idx);
  interp_idx = malloc(sizeof(int)*qke_struct.neq);
  for(i=0; i<qke_struct.neq; i++){
    interp_idx[i] = 1;
    //printf("%g ",y_inout2[i]);
  }
  printf("\n");
  
  start = clock();  
  func_return = generic_evolver(qke_derivs_fixed_grid_low_temp,
				T_lowtemp,
				qke_struct.T_final,
				y_inout2, 
				interp_idx,
				qke_struct.neq, 
				&qke_struct,
				qke_struct.rtol, 
				qke_struct.abstol, 
				qke_struct.Tvec, 
				qke_struct.Tres,
				NULL,//qke_struct.Ap,
				NULL,//qke_struct.Ai,
				qke_store_output_fixed_grid_low_temp,
				qke_print_variables_fixed_grid_low_temp,
				NULL,//qke_stop_at_L,
				error_message);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("CPU time used: %g minutes.\n",cpu_time_used/60.0);    
  

  free(y_inout2);
  */
  free(interp_idx);
  free_qke_param(&qke_struct);
  background_free_dof(&qke_struct.pbs);


  if (func_return == _FAILURE_){
    printf("%s\n",error_message);
    return _FAILURE_;
  }
  else{
    return _SUCCESS_;
  }
}
