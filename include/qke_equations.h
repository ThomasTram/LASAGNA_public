#ifndef __QKE__
#define __QKE__
#include "common.h"
#include "evolver_ndf15.h"
#include "newton.h"
#include "background.h"
#include "mat_io.h"
/**************************************************************/

typedef struct qke_param_structure{
  char output_filename[_FILENAMESIZE_]; //Where to write output.
  int evolver;   //Which time integrator to use
  int nproc;     //Number of cores available.
  int Nres;      //Number of resinances
  int neq;       //Number of equations
  int Tres;      //Entries in time/Temperature vector.
  int is_electron; //True if we have electron neutrino, False otherwise.
  int evolve_vi;
  double rtol;   //Relative tolerance of integrator
  double abstol; //Absolute tolerance of integrator
  double alpha;  //Sampling density aound resonances (0=most dense, 1=uniform)
  double rs;     //Regulator for sterile sector: (0 turns off regulation)
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
  double *y_0;   //Workspace for qke_derivs, for Newton method.
  double *maxstep;
  int vres;      //Number of momentum bins in v-space. (Resolution)
  int v_left;    //Boundaries of v, usually just 0 and 1.
  int v_right;
  double *v_grid;  //v_grid[vres]
  double *u_grid;
  double *x_grid;
  double *dvdu_grid;
  double *dudT_grid;
  double xmin;     //minimum value of x=p/T considered
  double xmax;     //maximum value of x=p/T considered
  double xext;     //extremum value of distribution, xext \simeg 3.1 for FD
  double eps1;     //Function of xmin,xmax and xext
  double eps2;     //Function of xmin,xmax and xext
  int index_Pa_plus; //Start index of P_a+ in y-vector
  int index_Pa_minus;
  int index_Ps_plus;
  int index_Ps_minus;
  int index_Px_plus;
  int index_Px_minus;
  int index_Py_plus;
  int index_Py_minus;
  int index_L;       //Index of lepton asymetry in y-vector
  int index_b;
  int index_vi;
  double theta_zero; //Active sterile mixing angle
  double delta_m2;		//squared mass difference
  double L_initial;
  double L_final; //Final value of abs(L)
  double mu_div_T_initial;
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
} qke_param;

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
  int u_of_x(double x, double *u, double *dudx, qke_param *param);
  int init_qke_param(qke_param *pqke);
  int free_qke_param(qke_param *pqke);
  int get_resonances_xi(double T, double L, qke_param *param);
  int x_of_u(double u, double *x, qke_param *param);
  int nonlinear_rhs(double *y, double *Fy, void *param);
  int qke_initial_conditions(double Ti, double *y, qke_param *pqke);
  int qke_init_output(qke_param *pqke);
  int qke_store_output(double t,
		       double *y,
		       double *dy,
		       int index_t,
		       void *pqke,
		       ErrorMsg error_message);
  int qke_print_variables(double t,
			double *y,
			double *dy,
			void *pqke,
			  ErrorMsg error_message);
  int qke_derivs(double T, 
		 double *y, 
		 double *dy, 
		 void *pqke, 
		 ErrorMsg error_message);

  



#ifdef __cplusplus
}
#endif

/**************************************************************/
//Some physical quantities in units where c=hbar=k_b=1
#define _G_F_ 1.16637e-5		//GeV^-2 The Fermi coupling constant
#define _SIN2_THETA_W_ 0.23120	//sin(theta_W)^2, the square of the sine of the Weinberg angle
#define _M_W_ 80.398			//GeV. Mass of the W boson
#define _M_Z_ 91.1876 			//GeV. Mass of the Z boson
#define _M_PL_ 1.2209e19                //GeV. Planck mass sqrt(hbar*c/G)
//Code parameters:
#define _L_SCALE_ 1e-15                //The scale of L
#define _ZETA3_ 1.20205690315959428539973816151
#endif
