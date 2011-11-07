#ifndef _BG_
#define _BG_


#include "common.h"
#include "arrays.h"
struct background_structure{
  double *T;       //Temperature array
  double *dof_array;     //DoF spline array
  int ndof;        //Entries in DoF spline array
  char dof_filename[_FILENAMESIZE_]; //File where data is read..
  int klo;         //Last interpolation indices in dof_array
  int khi;
  int Tdir;        //1 for increasing, -1 for decreasing.
};

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int background_free_dof(struct background_structure *pbs);
  int background_init_dof(struct background_structure *pbs);
  int background_getdof(double T, 
			double *sqrth, 
			double *gentr, 
			struct background_structure *pbs);


#ifdef __cplusplus
}
#endif


#endif
