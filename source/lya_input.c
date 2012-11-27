#include "lya_input.h"
int lya_input_init_from_arguments(int argc, 
				  char **argv,
				  lya_param *plya,
				  ErrorMsg errmsg){
  struct file_content fc;
  fc.size = 0;

  lasagna_test(argc>2, 
	       errmsg,
	       "Too many input arguments to LASAGNA!");

  if (argc==2){
    lasagna_call(parser_read_file(argv[1],&fc,errmsg),
		 errmsg,
		 errmsg);
  }
  
  lasagna_call(lya_input_init(&fc,
			      plya,
			      errmsg),
	       errmsg,
	       errmsg);
  
  lasagna_call(parser_free(&fc),errmsg,errmsg);


  return _SUCCESS_;
}

int lya_input_init(struct file_content *pfc,
		   lya_param *plya,
		   ErrorMsg errmsg){
  int flag1,flag2;
  double param1,param2;
  int int1;
  char string1[_ARGUMENT_LENGTH_MAX_];

  //Read default values:
  lasagna_call(lya_input_default_params(plya),
	       errmsg,
	       errmsg);
  
  lasagna_read_string("dof_filename",plya->pbs.dof_filename);
  lasagna_read_double("xext",plya->xext);
  lasagna_read_double("xmin",plya->xmin);
  lasagna_read_double("xmax",plya->xmax);
  lasagna_read_int("evolver",plya->evolver);	
  lasagna_read_int("linalg_wrapper",plya->LinearAlgebraWrapper);
  lasagna_read_int("nproc",plya->nproc);
  lasagna_read_int("verbose",plya->verbose);
  lasagna_read_double("T_initial",plya->T_initial);
  lasagna_read_double("T_final",plya->T_final);
  lasagna_read_double("Tres",plya->Tres);
  lasagna_read_double("v_left",plya->v_left);
  lasagna_read_double("v_right",plya->v_right);
  lasagna_read_double("vres",plya->vres);
  lasagna_read_double("rtol",plya->rtol);
  lasagna_read_double("abstol",plya->abstol);
  lasagna_read_string("output_filename",plya->output_filename);
  lasagna_read_double("rs",plya->rs);
  lasagna_read_double("alpha",plya->alpha);
  lasagna_read_double("L_initial",plya->L_initial);
  lasagna_read_double("L_final",plya->L_final);
  lasagna_read_double("delta_m2",plya->delta_m2);
  lasagna_read_double_one_of_two("sinsq2theta","theta_zero",plya->theta_zero);
  if (flag1 == _TRUE_)
    plya->theta_zero = 0.5*asin(sqrt(plya->theta_zero));
  lasagna_read_double("trigger_dLdT_over_L",plya->trigger_dLdT_over_L);
  lasagna_read_int("is_electron",plya->is_electron);
  if (plya->is_electron == _TRUE_)
    printf("Flavour of active species: Electron\n");
  else
    printf("Flavour of active species: Muon/Tau\n");
  lasagna_read_int("fixed_grid", plya->fixed_grid);
  lasagna_read_int("lyapunov_seed",plya->lyapunov_seed);
  lasagna_read_double("v_scale",plya->v_scale);
  lasagna_read_double("T_wait",plya->T_wait);
  lasagna_read_double("I_stop",plya->I_stop);

  //Initialise background somewhere
  background_init_dof(&(plya->pbs));
  //Initialise rest of plya 
  init_lya_param(plya);


  return _SUCCESS_;
}

int lya_input_default_params(lya_param *plya){
  double sinsq2theta;
  strncpy(plya->pbs.dof_filename,"dsdofHP_B.dat",_FILENAMESIZE_);
  plya->xext = 3.1;
  plya->xmin =  0.0001;// 0.0001; //1e-4;
  plya->xmax = 100.0; //100.0;
  plya->evolver = 1;
  plya->LinearAlgebraWrapper = LINALG_WRAPPER_SPARSE;
  plya->nproc = 1;
  plya->verbose = 4;
  plya->fixed_grid = 0;
  plya->Nres = 2;
  plya->T_initial = 0.025;
  plya->T_final = 0.010;
  plya->Tres = 500;

  plya->v_right = 1.0;
  plya->v_left = 0.0;
  plya->vres = 200;
  strcpy(plya->output_filename,"output/dump.mat");
  /** We must have non-zero alpha, otherwise the matrix for 
      solving for dvidT becomes singular.
  */
  plya->alpha = 0.1;
  plya->rs = 0.0;
  plya->L_initial = 2e-10;
  plya->L_final = 5e-5;
  plya->delta_m2 = -1e-19;
  sinsq2theta = pow(10.0,-9.0);
  plya->theta_zero = 0.5*asin(sqrt(sinsq2theta));
  plya->is_electron = _FALSE_;
  plya->rtol = 1e-3;
  plya->abstol = 1e-6;

  plya->trigger_dLdT_over_L = 1e100;
  plya->lyapunov_seed = 1;
  plya->v_scale = 1e10;
  plya->T_wait = -1; //Deactivate stop_at_divL
  plya->I_stop = 100;
  return _SUCCESS_;
}
