#ifndef __MULTIMATRIX /* allow multiple inclusions */
#define __MULTIMATRIX
#include "common.h"
/********************************************
 * The matrix types are defined as follows. *
 ********************************************/
typedef enum {
  L_SCC,   /* Sparse Compressed Column format */
  L_DNR,   /* Dense Numerical Recipes style */ 
} StoreType;

typedef enum {
  L_DBL,   /* double */
  L_DBL_CX,/* double complex */
} DataType;

typedef struct {
  StoreType Stype; /* Storage type: interprets the storage structure 
		      pointed to by *Store. */
  DataType Dtype; /* Data type. */
  int_t  nrow;   /* number of rows */
  int_t  ncol;   /* number of columns */
  void *Store;   /* pointer to the actual storage of the matrix */
} MultiMatrix;

/***********************************************
 * The storage schemes are defined as follows. *
 ***********************************************/

/* StoreType == L_SCC (Also known as Harwell-Boeing sparse matrix format) */
typedef struct {
  int  nnz;	    /* number of nonzeros in the matrix */
  int  *Ai;         /* pointer to array of row indices of the nonzeros */
  int  *Ap;         /* pointer to array of beginning of columns in nzval[] 
		       and rowind[]  */
  void *Ax;         /* pointer to array of nonzero values, packed by column */
  
 /* Note:
     Zero-based indexing is used;
     colptr[] has ncol+1 entries, the last one pointing
     beyond the last column, so that colptr[ncol] = nnz. */
} SCCformat;

/* StoreType == L_DNR (This is the storage format advocated in Numerical 
   Recipes, with 1-based indexing.) */
typedef struct {
  void * Matrix;   /* This is an array of pointers to the 
		      rows of the matrix. Because of the choice 
		      of 1-based indexing, it is (nrow+1) long.*/
  void * Data;     /* This is a pointer to (ncol*nrow+1) values,
		      and the matrix itself is stored in 
		      row-major order starting from 
		      pointer_to_data + 1. */
} DNRformat;
