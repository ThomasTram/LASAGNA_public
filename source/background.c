#include "background.h"
int background_free_dof(struct background_structure *pbs){
  free(pbs->T);
  free(pbs->dof_array);
  return _SUCCESS_;
}

int background_init_dof(struct background_structure *pbs){
  double tmp1,tmp2,tmp3;
  FILE *datafile;
  int row,status,tablesize,firstrow=0,header_found;
  char tmpstring[256];
  char *ptr_string=tmpstring;
  datafile = fopen(pbs->dof_filename,"r");
  if (datafile==NULL) printf("Null...\n");
  // Find size of table:
  for (row=0,header_found = 0; ; row++){
    ptr_string = fgets(tmpstring,256,datafile);
    if (ptr_string==NULL){
      tablesize = row;
      break;
    }
    else if(header_found == 0){
      //If header has not been found, search for it:
      status = sscanf(tmpstring, "%lf %lf %lf",&tmp1,&tmp2,&tmp3);
      if (status==3){
	//Header found.
	header_found = 1;
	firstrow = row;
      }
    }
  }
  rewind(datafile);
  for (row=0; row<firstrow; row++){
    ptr_string = fgets(tmpstring,256,datafile);
    printf("%s\n",tmpstring);
  }
  tablesize -= firstrow;
  pbs->ndof = tablesize;
  pbs->klo = 0;
  pbs->khi = pbs->ndof-1;
  //Allocate T andy:
  pbs->T = malloc(sizeof(double)*tablesize);
  pbs->dof_array = malloc(sizeof(double)*tablesize*2*2);

  for (row=0; row<tablesize; row++){
    status = fscanf(datafile,"%lf %lf %lf",
		    pbs->T+row,
		    pbs->dof_array+row,
		    pbs->dof_array+2*tablesize+row);
  }
  
  arrays_spline_natural(pbs->T,
			pbs->dof_array,
			pbs->dof_array+tablesize,
			tablesize);
  arrays_spline_natural(pbs->T,
			pbs->dof_array+2*tablesize,
			pbs->dof_array+3*tablesize,
			tablesize);
  if (pbs->T[tablesize-1]>=pbs->T[0])
    pbs->Tdir=1;
  else
    pbs->Tdir=-1;    
  fclose(datafile);
  return _SUCCESS_;
}


int background_getdof(double T, 
		      double *sqrtg, 
		      double *gentr, 
		      struct background_structure *pbs){
  /** Interpolation table for the square root of radiation dof and entropy dof as a funtion of the temperature T in GeV */
  int n=pbs->ndof;
  if (sqrtg!=NULL){
    arrays_spline_interpolate(pbs->T,
			      pbs->dof_array,
			      pbs->dof_array+n,
			      n,
			      pbs->Tdir,
			      T,
			      sqrtg,
			      &(pbs->klo),
			      &(pbs->khi));
  }
  if (gentr!=NULL){
    arrays_spline_interpolate(pbs->T,
			      pbs->dof_array+2*n,
			      pbs->dof_array+3*n,
			      n,
			      pbs->Tdir,
			      T,
			      gentr,
			      &(pbs->klo),
			      &(pbs->khi));
  }
 return _SUCCESS_;
}


