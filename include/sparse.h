#ifndef __SPA__
#define __SPA__
/****************************************/
/* Sparse Matrix algorithms for CLASS   */
/* 15/11 2010                           */
/* Thomas Tram                          */
/****************************************/
#include "common.h"
#include "complex.h"

/* Structures: */
typedef struct sparse_matrix{
	/* Sparse matrix in compressed column form: */
	int ncols;		/* Number of columns */
	int nrows;		/* Number of rows */
	int maxnz;		/* Maximum number of non-zero entries*/
	int *Ap;		/* Ap[0..ncols]. Ap[k+1]-Ap[k] is the number of entries in the k'th column. */
	int *Ai;		/* Ai[0..(maxnz-1)]. Contains the row indices of the entries. */
	double *Ax;		/* Ax[0..(maxnz-1)]. Contains the values of the entries. */
} sp_mat;

typedef struct sparse_numerical{
	/* Sparse LU decomposition along with enough information to do a fast refactorization: */
	int n;			/*Matrix assumed square, [nxn] */
	sp_mat *L;		/*L and U is the factors of the decomposed matrix.*/
	sp_mat *U;
	int **xi;		/*xi[k] points to a row of xi, which holds the topological ordered indices.*/
	int *topvec;	/*topvec[k] holds the first index in xi[k].*/
	int *pinv;		/*Inverse row permutation. */
	int *p;			/*Row permutation. */
	int *q;			/* Column permutation */
	int *wamd;		/* Work array for sp_amd */
	double *w;		/* Work array for sp_lu */
} sp_num;

typedef struct sparse_matrix_complex{
	/* Sparse matrix in compressed column form: */
	int ncols;		/* Number of columns */
	int nrows;		/* Number of rows */
	int maxnz;		/* Maximum number of non-zero entries*/
	int *Ap;		/* Ap[0..ncols]. Ap[k+1]-Ap[k] is the number of entries in the k'th column. */
	int *Ai;		/* Ai[0..(maxnz-1)]. Contains the row indices of the entries. */
	double complex *Ax;		/* Ax[0..(maxnz-1)]. Contains the values of the entries. */
} sp_mat_cx;

typedef struct sparse_numerical_complex{
	/* Sparse LU decomposition along with enough information to do a fast refactorization: */
	int n;			/*Matrix assumed square, [nxn] */
	sp_mat_cx *L;		/*L and U is the factors of the decomposed matrix.*/
	sp_mat_cx *U;
	int **xi;		/*xi[k] points to a row of xi, which holds the topological ordered indices.*/
	int *topvec;	/*topvec[k] holds the first index in xi[k].*/
	int *pinv;		/*Inverse row permutation. */
	int *p;			/*Row permutation. */
	int *q;			/* Column permutation */
	int *wamd;		/* Work array for sp_amd */
	double complex *w;		/* Work array for sp_lu */
} sp_num_cx;


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
/* Routines and macros: */
  int sp_mat_alloc(sp_mat** A, int ncols, int nrows, int maxnz, ErrorMsg error_message);
  int sp_mat_free(sp_mat *A);
  int sp_num_alloc(sp_num** N, int n,ErrorMsg error_message);
  int sp_num_free(sp_num *N);
  int reachr(int Gncol, 
	     int *Bp, 
	     int *Bi, 
	     int *Gp, 
	     int *Gi, int k, int *xik,int *pinv);
  void dfsr(int j, int *Gp, int *Gi, int *top, int *xik, int *pinv);
  int sp_splsolve(sp_mat *G, sp_mat *B, int k, int*xik, int top, double *x, int *pinv);
  int sp_ludcmp(sp_num *N, sp_mat *A, double pivtol);
  int sp_lusolve(sp_num *N, double *b, double *x);
  int sp_refactor(sp_num *N, sp_mat *A);
  int column_grouping(sp_mat *G, int *col_g, int *col_wi);
  int column_grouping2(sp_mat *G, int *col_g, int *col_wi);
  int sp_amd(int *Cp, int *Ci, int n, int cnzmax, int *P, int *W);
  int sp_wclear(int mark, int lemax, int *w, int n);
  int sp_tdfs(int j, int k, int *head, const int *next, int *post, int *stack);

  int sp_mat_alloc_cx(sp_mat_cx** A, int ncols, int nrows, int maxnz, ErrorMsg error_message);
  int sp_mat_free_cx(sp_mat_cx *A);
  int sp_num_alloc_cx(sp_num_cx** N, int n,ErrorMsg error_message);
  int sp_num_free_cx(sp_num_cx *N);
  int sp_splsolve_cx(sp_mat_cx *G, sp_mat_cx *B, int k, int*xik, int top, double complex *x, int *pinv);
  int sp_ludcmp_cx(sp_num_cx *N, sp_mat_cx *A, double pivtol);
  int sp_lusolve_cx(sp_num_cx *N, double complex *b, double complex *x);
  int sp_refactor_cx(sp_num_cx *N, sp_mat_cx *A);
  int get_pattern_A_plus_AT(int *Ap, 
			    int *Ai, 
			    int n, 
			    int **Cp, 
			    int **Ci, 
			    ErrorMsg error_message);
  
#define SPFLIP(i) (-(i)-2)
#define SPUNFLIP(i) (((i)<0) ? SPFLIP(i) : (i))
#define SPMARKED(w,j) (w[j] < 0)
#define SPMARK(w,j) {w[j] = SPFLIP(w[j]);}

#ifdef __cplusplus
}
#endif


#endif
