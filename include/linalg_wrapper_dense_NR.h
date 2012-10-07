#include "common.h"
#include "multimatrix.h"

typedef struct {
  MultiMatrix *LU; //The LU decomposition.
  void *LUw;       //Workspace for LU
  int *luidx;      //Permutation vector for LU.
  int neq;
} DNR_structure;
