#include "input.h"
int input_init_from_arguments(int argc, 
			      char **argv,
			      qke_param *pqke,
			      ErrorMsg errmsg){
  int i;
  struct file_content fc;
  char input_file[_ARGUMENT_LENGTH_MAX_];
  fc.size = 0;
  input_file[0]='\0';

  lasagna_test(argc>2, 
	       errmsg,
	       "Too many input arguments to LASAGNA!");

  if (argc==2){
    lasagna_call(parser_read_file(argv[1],&fc,errmsg),
		 errmsg,
		 errmsg);
  }
  
  lasagna_call(input_init(&fc,
			  pqke,
			  errmsg),
	       errmsg,
	       errmsg);
  
  lasagna_call(parser_free(&fc),errmsg,errmsg);


  return _SUCCESS_;
}

int input_init(struct file_content *pfc,
	       qke_param *pqke,
	       ErrorMsg errmsg){
  int flag1,flag2,flag3;
  double param1,param2,param3;
  int entries_read;
  int int1, fileentries;
   double * pointer1;
  char string1[_ARGUMENT_LENGTH_MAX_];

  double sinsq2theta;

  //Read default values:
  lasagna_call(input_default_params(pqke),
	       errmsg,
	       errmsg);
  
  lasagna_read_string("dof_filename",pqke->pbs.dof_filename);
  lasagna_read_double("xext",pqke->xext);
  lasagna_read_double("xmin",pqke->xmin);
  lasagna_read_double("xmax",pqke->xmax);
  lasagna_read_int("evolve_vi",pqke->evolve_vi);	
  lasagna_read_int("evolver",pqke->evolver);	
  lasagna_read_double("T_initial",pqke->T_initial);
  lasagna_read_double("T_final",pqke->T_final);
  lasagna_read_double("Tres",pqke->Tres);
  lasagna_read_double("v_left",pqke->v_left);
  lasagna_read_double("v_right",pqke->v_right);
  lasagna_read_double("vres",pqke->vres);
  lasagna_read_double("rtol",pqke->rtol);
  lasagna_read_double("abstol",pqke->abstol);
  lasagna_read_string("output_filename",pqke->output_filename);
  lasagna_read_double("rs",pqke->rs);
  lasagna_read_double("alpha",pqke->alpha);
  lasagna_read_double("L_initial",pqke->L_initial);
  lasagna_read_double("L_final",pqke->L_final);
  lasagna_read_double("delta_m2",pqke->delta_m2);
  lasagna_read_double_one_of_two("sinsq2theta","theta_zero",pqke->theta_zero);
  if (flag1 == _TRUE_)
    pqke->theta_zero = 0.5*asin(sqrt(pqke->theta_zero));
  lasagna_read_int("is_electron",pqke->is_electron);

  //Initialise background somewhere
  background_init_dof(&(pqke->pbs));
  //Initialise rest of pqke somewhere
  init_qke_param(pqke);
  return _SUCCESS_;
}

int input_default_params(qke_param *pqke){
  double sinsq2theta;
  strncpy(pqke->pbs.dof_filename,"dsdofHP_B.dat",_FILENAMESIZE_);
  pqke->xext = 3.1;
  pqke->xmin =  0.0001;// 0.0001; //1e-4;
  pqke->xmax = 100.0; //100.0;
  pqke->evolve_vi = _FALSE_;
  pqke->evolver = 0;
  pqke->Nres = 2;
  pqke->T_initial = 0.040;
  pqke->T_final = 0.010;
  pqke->Tres = 10000;

  pqke->v_right = 1.0;
  pqke->v_left = 0.0;
  pqke->vres = 500;
  strcpy(pqke->output_filename,"output/dump.mat");
  /** We must have non-zero alpha, otherwise the matrix for 
      solving for dvidT becomes singular.
  */
  pqke->alpha = 0.2;
  pqke->rs = 0.0;
  pqke->L_initial = 2e-10;
  pqke->L_final = 5e-5;
  pqke->delta_m2 = -1e-17; //-10e-18;
  sinsq2theta = pow(10.0,-6.0);
  pqke->theta_zero = 0.5*asin(sqrt(sinsq2theta));
  pqke->is_electron = _FALSE_;
  pqke->rtol = 1e-2;
  pqke->abstol = 1e-5;
  return _SUCCESS_;
}
