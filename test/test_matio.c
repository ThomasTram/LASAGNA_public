#include "test_matio.h"
int main(){
  double data[6];
  float fl_data[6];
  int i, idata[4];
  char my_filename[100]="testmat.mat";
  for (i=0; i<6; i++){
    data[i] = (i+1.0)*(i+2.0);
    fl_data[i] = data[i];
    printf("%g, ",data[i]);
  }
  for (i=0; i<4; i++)
    idata[i] = 2*i+1;
  printf("Size of double: %d, Size of single: %d.\n",
	 sizeof(data[0]),sizeof(fl_data[0]));
  printf("String length is: %d\n",strlen(my_filename));
  mat_create_file(my_filename);
  mat_add_matrix("testmat.mat","matrix1",miDOUBLE,2,3);
  mat_add_matrix("testmat.mat","matrix2_I_am_soo_good",miDOUBLE,5,1);
  mat_add_matrix("testmat.mat","matrix3",miDOUBLE,8,2);
  mat_add_matrix("testmat.mat","intmatrix",miINT32,2,5); 
  mat_add_matrix("testmat.mat","matrix4",miSINGLE,6,6);
 
  mat_write_data("testmat.mat","matrix3",data,6,6);
  mat_write_data("testmat.mat","intmatrix",idata,4,4);
  mat_write_data("testmat.mat","matrix4",fl_data,30,6);
  
  return _SUCCESS_;
}
