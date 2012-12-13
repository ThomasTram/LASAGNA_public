#ifndef __LYA_EQUATIONS__
#define __LYA_EQUATIONS__
#include "common.h"
#include "evolver_ndf15.h"
#include "evolver_common.h"
#include "newton.h"
#include "background.h"
#include "qke_equations.h"
#include "mat_io.h"
#include "multimatrix.h"

typedef struct lya_param_structure{
  FILE *tmp;
  char output_filename[_FILENAMESIZE_]; //Where to write output.
  int evolver;   //Which time integrator to use
  LinAlgWrapper LinearAlgebraWrapper; //Wrapper for LU-decompositions
  int nproc;     //Number of cores available.
  int verbose;   //Level of output
  int fixed_grid;//Fixed grid?
  int Nres;      //Number of resinances
  size_t neq;       //Number of equations
  int Tres;      //Entries in time/Temperature vector.
  int is_electron; //True if we have electron neutrino, False otherwise.
  int guess_exists;
  double rtol;   //Relative tolerance of integrator
  double abstol; //Absolute tolerance of integrator
  double alpha;  //Sampling density aound resonances (0=most dense, 1=uniform)
  double rs;
  double T_initial; //Initial temperature
  double T_final; //Final temperature
  double *Tvec;  //Temperature vector (Tvec[Tres])
  double *xi;    //Resonances in v-space.
  double *ui;    //Resonances in u-space, ui[Nres]
  double *vi;
  double *duidT; //Partial derivative of resonances in u-space.
  double *dvidT; //Partial derivative of resonances in v-space.
  double *duidx; //Partial derivative of u wrt x
  double *dxidT; //Partial derivative of resonances in x space.
  double *a;     //a and b parameters controlling the parametrisation.
  double b;
  double *y_0;   //Workspace for lya_derivs, for Newton method.
  double *maxstep;
  int vres;      //Number of momentum bins in v-space. (Resolution)
  double v_left;    //Boundaries of v, usually just 0 and 1.
  double v_right;
  double *v_grid;  //v_grid[vres]
  double *u_grid;
  double *x_grid;
  double *dvdu_grid;
  double *dudT_grid;
  double n_plus;
  double xmin;     //minimum value of x=p/T considered
  double xmax;     //Maximum value of x=p/T considered
  double xext;     //extremum value of distribution, xext \simeg 3.1 for FD
  double eps1;     //Function of xmin,xmax and xext
  double eps2;     //Function of xmin,xmax and xext
  //Common indices:
  int index_L;       //Index of lepton asymetry in y-vector
  int index_n_plus; 
 //Indices for normal method:
  int index_Pa_plus; //Start index of P_a+ in y-vector
  int index_Pa_minus;
  int index_Ps_plus;
  int index_Ps_minus;
  int index_Px_plus;
  int index_Px_minus;
  int index_Py_plus;
  int index_Py_minus;
  //Model parameters:
  double theta_zero; //Active sterile mixing angle
  double delta_m2;		//squared mass difference
  double L_initial;
  //Trigger parameters:
  double L_final; //Final value of abs(L)
  double trigger_dLdT_over_L;
  double C_alpha;
  double **mat;      //(Nres+2)x(Nres+2) matrix used in more than one occasion for solving linear systems.
  double *vv;        //Workarray for LU decomposition.
  int *indx;         //Permutation vector for LU decomposition.
  int *Ai;  //Row indices of jacobian
  int *Ap;  //Column indices of jacobian
  struct background_structure pbs;
  double Vx;
  double V0;
  double V1;
  double VL;
  double g_alpha; //1 for muon/tau, 1+4sec²(theta_w)/n_plus for e.
  int Pa_plus_handle; //Position of given matrix in output file
  int Pa_minus_handle; //Position of given matrix in output file
  int Ps_plus_handle; //Position of given matrix in output file
  int Ps_minus_handle; //Position of given matrix in output file
  int Px_plus_handle; //Position of given matrix in output file
  int Px_minus_handle; //Position of given matrix in output file
  int Py_plus_handle; //Position of given matrix in output file
  int Py_minus_handle; //Position of given matrix in output file
  int x_grid_handle; //Position of given matrix in output file
  int u_grid_handle; //Position of given matrix in output file
  int v_grid_handle; //Position of given matrix in output file
  int b_a_vec_handle; //Position of given matrix in output file 
  int xi_handle; //Position of given matrix in output file
  int ui_handle; //Position of given matrix in output file
  int vi_handle; //Position of given matrix in output file
  int L_handle; //Position of given matrix in output file
  int T_handle; //Position of given matrix in output file
  int I_conserved_handle; //Position of given matrix in output file
  int V0_handle; //Position of given matrix in output file
  int V1_handle; //Position of given matrix in output file
  int Vx_handle; //Position of given matrix in output file
  int VL_handle; //Position of given matrix in output file
  // For qke_stop_at_divL.
  double T_wait;
  double max_old;
  double max_cur;
  int should_break;
  double breakpoint;
  double I_stop;
  // For lyapunov calculation only:
  double v_scale; //Scale to suppres the lyapunov vectors.
  int lyapunov_seed; //For the random start vector if running lasagna_lya.
  int index_v;
  int v_handle;
  int I_handle;
  MultiMatrix **J_pp;
} lya_param;

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  //Initialise:
  int init_lya_param(lya_param *plya);
  //Free:
  int free_lya_param(lya_param *plya);
  //Set or copy initial conditions:
  int lya_initial_conditions(double Ti, double *y, lya_param *plya);
  //Handle binary output:
  int lya_init_output(lya_param *plya);
  int lya_store_output(double t,
		       double *y,
		       double *dy,
		       int index_t,
		       void *plya,
		       ErrorMsg error_message);

  //Handle ASCII output for verbose mode:
  int lya_print_variables(double t,
			  double *y,
			  double *dy,
			  void *plya,
			  ErrorMsg error_message);
  //LYA derivs:
  int lya_derivs(double T, 
		 double *y, 
		 double *dy, 
		 void *plya, 
		 ErrorMsg error_message);

  int lya_print_L(double t,
		  double *y,
		  double *dy,
		  void *plya,
		  ErrorMsg error_message);
  //Stop conditions:
  int lya_stop_at_L(double t,
		    double *y,
		    double *dy,
		    void *param,
		    ErrorMsg error_message);
  int lya_stop_at_divL(double t,
		       double *y,
		       double *dy,
		       void *param,
		       ErrorMsg error_message);
  int lya_stop_at_trigger(double t,
			  double *y,
			  double *dy,
			  void *param,
			  ErrorMsg error_message);
  //Misc.:
  int lya_get_resonances_xi(double T, 
			    double L,
			    lya_param *plya);
  int lya_get_resonances_dxidT(double T, 
			       double L,
			       double dLdT,
			       lya_param *plya);
  int lya_get_integrated_quantities(double *y,
				    lya_param *plya,
				    double *I_VxPy_minus,
				    double *I_f0Pa_plus,
				    double *I_rho_ss,
				    ErrorMsg error_message);
  int lya_get_parametrisation(double T,lya_param *plya, ErrorMsg error_message);
  int lya_get_partial_derivatives(double T,
				  double L,
				  double dLdT,
				  lya_param *plya, 
				  ErrorMsg error_message);
  int lya_u_of_x(double x, double *u, double *dudx, lya_param *param);
  int lya_x_of_u(double u, double *x, lya_param *param);
  int lya_nonlinear_rhs(double *y, double *Fy, void *param);
  int lya_derivs_test_partial(double T, 
			      double *y, 
			      double *dy, 
			      void *param,
			      ErrorMsg error_message);
  int lya_test_partial_output(double T,
			      double *y,
			      double *dy,
			      int index_t,
			      void *param,
			      ErrorMsg error_message);


#ifdef __cplusplus
}
#endif

#endif //__LYA_EQUATIONS__
