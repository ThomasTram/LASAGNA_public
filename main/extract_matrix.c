#include "extract_matrix.h"
int main(int argc, char *argv[]){
  int cols,rows,data_type;
  int i,j,k,elements;
  void *generic_pointer;
  double *dbl_ptr;
  float *float_ptr;
  int32_t *int32_ptr;
  int funcreturn;
  FILE *datfile;
  char *tmpstring;
  int name_length;

  //Get temperature to decide how much to plot:
  funcreturn = mat_read_data(argv[1],
			     "T",
			     &generic_pointer,
			     &cols,
			     &rows,
			     &data_type);
  if (funcreturn == _FAILURE_){
      printf("Couldn't extract temperature T from file %s..\n",argv[1]);
      return _FAILURE_;
  }
  dbl_ptr = generic_pointer;
  for (i=0; (i<cols)&&(dbl_ptr[i]!=0.0); i++);
  elements = i;
  free(generic_pointer);
  printf("%d out of %d output values computed.\n",elements,cols);
  for (i=2; i<argc; i++){
    //Loop over input arguments containing names of matrices:
    funcreturn = mat_read_data(argv[1],
			       argv[i],
			       &generic_pointer,
			       &cols,
			       &rows,
			       &data_type);
    if (funcreturn == _FAILURE_){
      printf("Couldn't extract matrix %s from file..\n",argv[i]);
      continue;
    }
    //Type cast pointer: (I could switch on data_type, but not neccesary.)
    dbl_ptr = generic_pointer;
    int32_ptr = generic_pointer;
    float_ptr = generic_pointer;

    name_length = strlen(argv[i]);
    tmpstring = malloc(sizeof(char)*(name_length+12));
    strcpy(tmpstring,"output/");
    strcat(tmpstring,argv[i]);
    strcat(tmpstring,".dat");
    printf("Writing [%d x (%d)%d] matrix %s to file %s.\n",
	   rows,min(elements,cols),cols,argv[i],tmpstring);
    datfile = fopen(tmpstring,"w");
    for (j=0; j<min(cols,elements); j++){
      for (k=0; k<rows; k++){
	switch (data_type){
	case miSINGLE:
	  fprintf(datfile,"%e ",float_ptr[j*rows+k]);
	  break;
	case miINT32:
	  fprintf(datfile,"%d ",int32_ptr[j*rows+k]);
	  break;
	case miDOUBLE:
	  fprintf(datfile,"%e ",dbl_ptr[j*rows+k]);
	  break;
	}
      }
      fprintf(datfile,"\n");
    }
    fclose(datfile);
    free(tmpstring);
    free(generic_pointer);
  }
  return _SUCCESS_;
}
