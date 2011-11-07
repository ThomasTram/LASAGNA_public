#include "common.h"

int besoutput(double t,
	      double *y,
	      double *dy, 
	      int index_t, 
	      void *param, 
	      ErrorMsg error_message);

int besderivs(
	      double t,
	      double *y,
	      double *dy,
	      void * ppaw,
	      ErrorMsg error_message);
