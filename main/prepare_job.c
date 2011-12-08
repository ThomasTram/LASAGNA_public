#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

int main(int argc, char **argv){
  char parameter_filename[150];
  char script_filename[150];
  char line_output[256];
  char line_cpin[256];
  char line_cpout[256];
  char line_sinsq[256];
  char line_deltam2[256];
  char line_dispatch[32768];
  char line_lasagna[256];
  char wrkdir[256];
  double sinsq_expl = -6.5, sinsq_expr = -5.5, sinsq_exp;
  double deltam2_expl = 0.5, deltam2_expr = 1.5, deltam2_exp;
  double deltam2, sinsq2theta;
  int i,j,deltares=36,sinsqres=36;
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
    "T_initial = 0.040\n"
    "L_initial = 2e-10\n"
    "T_final = 0.010\n"
    "L_final = 5e-5\n"
    "--- Output parameters ----------------\n"
    "Tres = 100\n"
    "--- Precision parameters -------------\n"
    "evolver = 0\n"
    "rtol = 1e-2\n"
    "abstol = 1e-5\n"
    "vres = 500\n"
    "alpha = 0.2\n"
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
    "#PBS -q qexp\n" 
    "#PBS -l nodes=1:ppn=8\n"
    "##PBS -l walltime=500:30:00\n";
   
  
  getcwd(wrkdir,sizeof(wrkdir));
  sprintf(line_dispatch,"dispatch -s %s/master.sh ",wrkdir);
  printf("%s\n", wrkdir);
  for (i=9; i<deltares; i++){
    deltam2_exp = deltam2_expl + i*(deltam2_expr-deltam2_expl)/(deltares-1.0);
    deltam2 = -pow(10,deltam2_exp)*1e-18;
    sprintf(line_deltam2,"delta_m2 = %.14e\n",deltam2);
    for (j=0; j<sinsqres; j++){
      sinsq_exp = sinsq_expl + j*(sinsq_expr-sinsq_expl)/(sinsqres-1.0);
      sinsq2theta = pow(10.0,sinsq_exp);
      sprintf(line_sinsq,"sinsq2theta = %.14e\n",sinsq2theta);
      //Write filename:
      sprintf(parameter_filename,"param_%d_%d.ini",i,j);
      sprintf(line_output,"output_filename = dump_%d_%d.mat\n",i,j);
      parameter_file = fopen(parameter_filename,"w");
      fprintf(parameter_file,"%s\n%s\n%s\n%s\n",
	      base_parameters,
	      line_deltam2,
	      line_sinsq,
	      line_output);
      fclose(parameter_file);
      sprintf(script_filename,"runlasagna_%d_%d",i,j);
      aux_file = fopen(script_filename,"w");
      sprintf(line_cpin,
	      "cp %s /scratch/$PBS_JOBID/$0\n",
	      parameter_filename);
      sprintf(line_lasagna,
	      "cd /scratch/$PBS_JOBID/$0\n./lasagna %s > $0.txt\n",
	      parameter_filename);
      sprintf(line_cpout,
	      "cp $0.txt $PBS_O_WORKDIR\ncp dump_%d_%d.mat $PBS_O_WORKDIR\n",
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
  printf("%s\n",line_dispatch);
  aux_file = fopen("master.sh","w");
  fprintf(aux_file,"%s\n",master_string);
  fclose(aux_file);
  aux_file = fopen("runloop.js","w");
  fprintf(aux_file,"%s\n%s\n",jobscript_string,line_dispatch);
  fclose(aux_file);
  return 0;
}
