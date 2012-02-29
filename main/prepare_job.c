#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

int main(int argc, char **argv){
  char parameter_filename[150];
  char script_filename[150];
  char jobscript_filename[150];
  char line_output[256];
  char line_cpin[256];
  char line_cpout[256];
  char line_sinsq[256];
  char line_deltam2[256];
  char line_T_initial[256];
  char line_dispatch[65536];
  char line_lasagna[256];
  char wrkdir[256];
  
     double sinsq_expl = -4.0, sinsq_expr = -0.8, sinsq_exp;
     double deltam2_expl = -3.0, deltam2_expr = 1.0, deltam2_exp;
  
  /**
  double sinsq_expl = -2.0, sinsq_expr = -1.0, sinsq_exp;
  double deltam2_expl = -1.0, deltam2_expr = 0.5, deltam2_exp;
  */  
  double deltam2, sinsq2theta;
  double L_initial=1e-2,T_initial;
  int sign_deltam2 = 1;
  int i,j,deltares=4,sinsqres=4;
  FILE *parameter_file;
  FILE *aux_file;
  //delta_m2 = 1e-17
  //sinsq2theta = 1e-6
  //output_filename = output/dump.mat
  char *base_parameters = 
    "--- Background parameters ------------\n"
    "dof_filename = dsdofHP_B.dat\n"
    "--- Oscillation parameters -----------\n"
    "is_electron = 0\n"
    "--- Evolution parameters -----------\n"
    "L_initial = 1e-2\n"
    "T_final = 0.001\n"
    "L_final = -1.0\n"
    "--- Output parameters ----------------\n"
    "Tres = 100\n"
    "--- Precision parameters -------------\n"
    "evolver = 2\n"
    "rtol = 1e-3\n"
    "abstol = 1e-6\n"
    "vres = 300\n"
    "alpha = 1.0\n"
    "xext = 3.1\n"
    "xmin = 1e-4\n"
    "xmax = 100.0\n"
    "evolve_vi = 0\n"
    "v_left = 0.0\n"
    "v_right = 1.0\n"
    "rs = 0.0\n";
  char *master_string = 
    "#!/bin/bash\n"
    "cd `dirname $0`\n"
    "./$1 &> $1.out\n";
  char *slave_string = 
    "#!/bin/bash\n"
    "mkdir /scratch/$PBS_JOBID/$0\n"
    "cp dsdofHP_B.dat /scratch/$PBS_JOBID/$0\n"
    "cp lasagna /scratch/$PBS_JOBID/$0\n";

  char *jobscript_string =
    "#!/bin/bash\n"
    "#PBS -q q8\n" 
    "#PBS -l nodes=2:ppn=8\n"
    "#PBS -l walltime=40:30:00\n"
    "echo \"========= Job started  at `date` ==========\"";
   char *jobscript_string_end = 
     "echo \"========= Job finished at `date` ==========\"";
  
  getcwd(wrkdir,sizeof(wrkdir));
  sprintf(line_dispatch,"dispatch -s %s/master.sh ",wrkdir);
  printf("%s\n", wrkdir);
  for (i=0; i<deltares; i++){
    deltam2_exp = deltam2_expl + i*(deltam2_expr-deltam2_expl)/(deltares-1.0);
    deltam2 = sign_deltam2*pow(10,deltam2_exp)*1e-18;
    T_initial = 1e-3*pow(0.1242*fabs(deltam2*1e18)/(fabs(L_initial)*0.1),0.25);
    T_initial +=1e-3*10.0;
    if (T_initial>0.040)
      T_initial = 0.040;
    if (T_initial<0.0011)
      T_initial = 0.0011;
    //T_initial=0.040;
    sprintf(line_deltam2,"delta_m2 = %.14e\n",deltam2);
    sprintf(line_T_initial,"T_initial = %.14e\n",T_initial);
    for (j=0; j<sinsqres; j++){
      sinsq_exp = sinsq_expl + j*(sinsq_expr-sinsq_expl)/(sinsqres-1.0);
      sinsq2theta = pow(10.0,sinsq_exp);
      sprintf(line_sinsq,"sinsq2theta = %.14e\n",sinsq2theta);
      //Write filename:
      sprintf(parameter_filename,"param_%03d_%03d.ini",i,j);
      sprintf(line_output,"output_filename = dump_%03d_%03d.mat\n",i,j);
      parameter_file = fopen(parameter_filename,"w");
      fprintf(parameter_file,"%s\n%s\n%s\n%s\n%s\n",
	      base_parameters,
	      line_T_initial,
	      line_deltam2,
	      line_sinsq,
	      line_output);
      fclose(parameter_file);
      sprintf(script_filename,"runlasagna_%03d_%03d",i,j);
      aux_file = fopen(script_filename,"w");
      sprintf(line_cpin,
	      "cp %s /scratch/$PBS_JOBID/$0\n",
	      parameter_filename);
      sprintf(line_lasagna,
	      "cd /scratch/$PBS_JOBID/$0\n./lasagna %s > $0.txt\n",
	      parameter_filename);
      sprintf(line_cpout,
	      "cp $0.txt $PBS_O_WORKDIR\ncp dump_%03d_%03d.mat $PBS_O_WORKDIR\n",
	      i,j);
      fprintf(aux_file,
	      "%s\n%s\n%s\n%s\n",
	      slave_string,
	      line_cpin,
	      line_lasagna,
	      line_cpout);
      fclose(aux_file);
      strcat(line_dispatch, script_filename);
      strcat(line_dispatch," ");
    }
  }
  aux_file = fopen("master.sh","w");
  fprintf(aux_file,"%s\n",master_string);
  fclose(aux_file);
  sprintf(jobscript_filename,"%s.js",argv[1]);
  aux_file = fopen(jobscript_filename,"w");
  if (aux_file == NULL){
    printf("Could not open %s\n",jobscript_filename);
    return 1;
  }
  fprintf(aux_file,"%s\n%s\n%s\n",
	  jobscript_string,
	  line_dispatch,
	  jobscript_string_end);
  fclose(aux_file);
  return 0;
}
