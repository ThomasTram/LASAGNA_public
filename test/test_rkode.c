#include "test_rkode.h"
#include "evolver_rk45.h"
#include "evolver_ndf15.h"
#include "evolver_radau5.h"
int main(){
  
  int i;
  double x0,xf;
  double y[2];
  int tres=600;
  double t_vec[tres];
  double nu;
  int NU;
  int useless[2];
  ErrorMsg error_message;

  useless[0] = 1;
  useless[1] = 1;
  nu=2.0;
  NU = 2;
  x0 = 1.0;
  xf = 50;
  y[0] = jn(NU,x0);
  y[1] = 0.5*(jn(NU-1,x0)-jn(NU+1,x0));
  
  for (i=0; i<tres; i++)
    t_vec[i] = x0+i*(xf-x0)/(tres-1);

  evolver_rkdp45(besderivs,
		 x0,
		 xf,
		 y,
		 useless,
		 2,
		 &nu,
		 1e-4,
		 1e-15,
		 NULL,//t_vec,
		 30,//tres,
		 NULL,
		 NULL,
		 besoutput,
		 NULL,
		 NULL,
		 error_message);
  return _SUCCESS_;
}

int besoutput(double t,
	      double *y,
	      double *dy, 
	      int index_t, 
	      void *param, 
	      ErrorMsg error_message){
  fprintf(stderr,"%.16e %.16e %.16e\n",t,y[0],y[1]);
  return _SUCCESS_;
}

int besderivs(
	      double t,
	      double *y,
	      double *dy,
	      void * ppaw,
	      ErrorMsg error_message){
  double *pnu = ppaw;
  double nu = *pnu;
  dy[0] = y[1];
  dy[1] = -y[1]/t-(1.0-nu*nu/t/t)*y[0];
  return _SUCCESS_;
}
