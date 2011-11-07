#include "common.h"
int arrays_spline_natural(double *x, double *y, double *d2y, int n){
  double *u;
  double p,sig;
  int i,k;

  u = malloc(sizeof(double)*n);
  d2y[0] = 0.0;
  u[0] = 0.0;
  for (i=1; i<n-1; i++){
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*d2y[i-1]+2.0;
    d2y[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  d2y[n-1] = 0.0;
  for(k=n-2; k>=0; k--){
    d2y[k] = d2y[k]*d2y[k+1]+u[k];
  }
  free(u);
  return _SUCCESS_;
}

int arrays_spline_interpolate(double *xa, 
			      double *ya, 
			      double *d2ya, 
			      int n,
			      int xdir,
			      double x,
			      double *y,
			      int klo,
			      int khi){
  
  int k,kk;
  double h,b,a;

  //  printf("Looking for x=%g in range [%g; %g]\n",x,xa[0],xa[n-1]);

  if ((x-xa[0])*xdir<=0.0){
    *y = ya[0];
    //  printf("Was before\n");
    return _SUCCESS_;
  }
  if ((xa[n-1]-x)*xdir<=0.0){
    *y = ya[n-1];
    // printf("Was after\n");
    return _SUCCESS_;
  }
 
  if ((klo<0)||(klo>n)){
    //klo out of bounds
    klo = 0;
  }
  else{
    //move klo such that it x is after xa[klo] in table:
    for (kk=1; klo==0;kk *=2){
      if ((x-xa[klo])*xdir<0.0)
	klo = max(klo-kk,0);
      else
	break;
    }
  }
  if ((khi<0)||(khi>n)){
    //khi out of bounds
    khi = n-1;
  }
  else{
    //move khi such that x is before x[khi] in table:
    for (kk=1; khi==n-1;kk *=2){
      if ((xa[khi]-x)*xdir<0.0)
	khi = min(khi+kk,n-1);
      else
	break;
    }
  }
  
  //printf("Is x=%g in x[%d]=%g -> x[%d]=%g.",
  // x,klo,xa[klo],khi,xa[khi]);
    //Do bisection to find bounds:
  while ((khi-klo)>1){
    k = (khi+klo) >> 1; //Divide by 2 using bit shift
    if ((xa[k] - x)*xdir>0.0) khi = k;
    else klo = k;
  }
    h = xa[khi]-xa[klo];
    a = (xa[khi]-x)/h;
    b = (x-xa[klo])/h;
    *y = a*ya[klo]+b*ya[khi] + ((a*a*a-a)*d2ya[klo]+(b*b*b-b)*d2ya[khi])*(h*h)/6.0;

  return _SUCCESS_;
}

