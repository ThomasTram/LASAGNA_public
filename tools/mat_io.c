#include "mat_io.h"
int mat_create_file(char *filename){
  FILE *mat_file;
  int16_t version=256;
  int32_t zero=0;
  int16_t endian=0x4d49;
  uint16_t zero8=0;
  int current_pos,pad_left;

  mat_file = fopen(filename,"wb");
  fprintf(mat_file,"%s",
	 "LASAGNA binary output file compatible with MATLAB.\n");
  current_pos = ftell(mat_file);
  for (pad_left=116-current_pos; pad_left>0; pad_left--)
    fwrite(&zero8,1,1,mat_file);
  if (pad_left<0)
    printf("Too large header, file corrupt!");

  //size_t fwrite ( const void * ptr, size_t size, size_t count, FILE * stream );
  //Write subsys data offset:
  fwrite(&zero, 4, 1, mat_file);
  fwrite(&zero, 4, 1, mat_file);
  //Write version and endian indicator:
  fwrite(&version, 2, 1, mat_file);
  fwrite(&endian, 2, 1, mat_file);
  fclose(mat_file);
  return _SUCCESS_;
}
  
int mat_add_matrix(char *filename, 
		   char *matrix_name,
		   int data_type,
		   int cols, 
		   int rows,
		   int *handle){
  
  FILE *mat_file;
  int32_t MAT_F_CLASS_T=0x00ff;
  int32_t array_class;
  int type_byte_size, data_byte_size, total_byte_size;
  int name_length,name_length_padded, data_byte_size_padded;
  int i;
  int8_t int8_tmp;
  int32_t int32_tmp;
  uint32_t uint32_tmp;
  uint64_t uint64_tmp;

  //Open file in append mode:
  mat_file = fopen(filename,"ab");
  switch(data_type){
  case miSINGLE:
    array_class = mxSINGLE_CLASS;
    type_byte_size = 4;
    break;
  case miINT32:
    array_class = mxINT32_CLASS;
    type_byte_size = 4;
    break;
  case miDOUBLE:
    array_class = mxDOUBLE_CLASS;
    type_byte_size = 8;
    break;
  default:
    array_class = mxDOUBLE_CLASS;
    type_byte_size = 8;
    printf("Warning: type not recognized!\n");
    break;
  }
  data_byte_size = cols*rows*type_byte_size;
  if ((data_byte_size%8)==0)
    data_byte_size_padded = data_byte_size;
  else
    data_byte_size_padded = data_byte_size+8-data_byte_size%8;
  
  //How much storage to set aside for the matrix name?
  name_length = strlen(matrix_name);
  if ((name_length%8)==0)
    name_length_padded = name_length;
  else
    name_length_padded = name_length+8-name_length%8;
  
  //Now calculate the total size including all subelements:
  total_byte_size = data_byte_size_padded + name_length_padded+6*8;
  /**  printf("name_length=%d, padded to %d, data_size = %d, padded to %d\n",
	 name_length,name_length_padded,data_byte_size,data_byte_size_padded);
  */
  //Write the tag of the MATLAB array:
  uint32_tmp = miMATRIX;
  fwrite(&uint32_tmp,4,1,mat_file);
  //Total size has been calculated:
  uint32_tmp = total_byte_size;
  fwrite(&uint32_tmp,4,1,mat_file);
  
  //Set subelement: array flags.
  //Tag
  uint32_tmp = miUINT32;
  fwrite(&uint32_tmp,4,1,mat_file);
  uint32_tmp = 8;
  fwrite(&uint32_tmp,4,1,mat_file);
  //Data
  int32_tmp = array_class & MAT_F_CLASS_T;
  fwrite(&int32_tmp,4,1,mat_file);
  int32_tmp=0;
  fwrite(&int32_tmp,4,1,mat_file);
  
  //Set subelement: Dimensions
  //Tag
  uint32_tmp = miINT32;
  fwrite(&uint32_tmp,4,1,mat_file);
  uint32_tmp = 8;
  fwrite(&uint32_tmp,4,1,mat_file);
  //Data
  int32_tmp = rows;
  fwrite(&int32_tmp,4,1,mat_file);
  int32_tmp = cols;
  fwrite(&int32_tmp,4,1,mat_file);

  //Set subelement: Array name
  //Tag
  uint32_tmp = miINT8;
  fwrite(&uint32_tmp,4,1,mat_file);
  uint32_tmp = name_length;
  fwrite(&uint32_tmp,4,1,mat_file);
  //Data
  for(i=0; i<name_length_padded; i++){
    if (i<name_length)
      int8_tmp = matrix_name[i];
    else
      int8_tmp = 0;
    fwrite(&int8_tmp,1,1,mat_file);
  }

  //Set subelement: Real part
  //Tag
  uint32_tmp = data_type;
  fwrite(&uint32_tmp, 4, 1, mat_file);
  uint32_tmp = data_byte_size;
  fwrite(&uint32_tmp, 4, 1, mat_file);
  //Data
  //Set handle do data-part:
  *handle = ftell(mat_file);
  //Fill with zeros:
  uint64_tmp = 0;
  for (i=0; i<(data_byte_size+7)/8; i++)
    fwrite(&uint64_tmp, 8, 1, mat_file);
  
  fclose(mat_file);
  return _SUCCESS_;
}

int mat_goto_matrix(FILE *mat_file,
		    char *matrix_name,
		    int *beginning_of_matrix){
  uint32_t matrix_size, name_size;
  int bytes_from_matrix_size_to_name_size = 36;
  int i,matrix_found=_FALSE_,start_data_element;
  int8_t int8_tmp;
  int read_status;

  /** This small helper routine finds the start of a matrix with a
      specific name */
  //Skip directly to the size of the first data element:
  fseek(mat_file,128,SEEK_SET);
  start_data_element = ftell(mat_file);

  while (matrix_found==_FALSE_){
    fseek(mat_file,4,SEEK_CUR);
    read_status = fread(&matrix_size, 4, 1, mat_file);

    //Skip to the size of array name:
    fseek(mat_file, bytes_from_matrix_size_to_name_size, SEEK_CUR);
    read_status = fread(&name_size, 4, 1, mat_file);
    if (read_status!=1){
      //There is no more matrices in file!
      printf("Error: Could not find matrix!!\n");
      if (feof(mat_file)==1)
	printf("Program has reached End-of-file!\n");
      return _FAILURE_;
    }

    //Is this the correct matrix?
    //Test for correct length, and if true, correct name:
    if (strlen(matrix_name)!=name_size){
      matrix_found = _FALSE_;
    }
    else{
      matrix_found=_TRUE_;
      for (i=0; i<name_size; i++){
	read_status = fread(&int8_tmp,1,1,mat_file);
	if (int8_tmp!=((int8_t) matrix_name[i])){
	  matrix_found = _FALSE_;
	  break;
	}
      }
    }
    if (matrix_found == _FALSE_){
      //Try to jump to next matrix:
      start_data_element +=(8+matrix_size);
      fseek(mat_file, start_data_element, SEEK_SET);
    }
    else{
      *beginning_of_matrix = start_data_element;
      return _SUCCESS_;
    }
  }
  return _FAILURE_;
}

int mat_write_data(char *filename, 
		   char *matrix_name, 
		   void *data_ptr,
		   int start_entry,
		   int entries){
  uint32_t array_size, data_type, name_size, name_size_padded;
  int type_byte_size;
  int start_data_element;
  FILE *mat_file;
  int func_return, read_status;
  
  //Open file for both reading and writing:
  mat_file = fopen(filename,"r+b");
  
  /** If a handle is not passed, check if the matrix name 
      exists in file, and return handle:
  */
  func_return = mat_goto_matrix(mat_file,
				  matrix_name,
				  &start_data_element);
  if (func_return == _FAILURE_){
    fclose(mat_file);
    return _FAILURE_;
  }
  name_size = strlen(matrix_name);
  if ((name_size%8)==0)
    name_size_padded = name_size;
  else
    name_size_padded = name_size + 8-name_size%8;
  //Store data:
  
  fseek(mat_file,start_data_element+6*8+name_size_padded,SEEK_SET);
  //Here is the beginning of the array. Read data type and array size:
  read_status = fread(&data_type,4,1,mat_file);
  read_status = fread(&array_size,4,1,mat_file);
  //Determine byte requirements for data-type:
  switch(data_type){
  case miSINGLE:
  case miINT32:
    type_byte_size = 4;
    break;
  case miDOUBLE:
    type_byte_size = 8;
    break;
  default:
    type_byte_size = 8;
    printf("Warning: type not recognized!\n");
    break;
  }
  //Test if matrix is large enough:
  if ((start_entry+entries)*type_byte_size>array_size){
    printf("Error: Not enough space in matrix! \n");
  }
  else{
    //Jump to where we want to write:
    fseek(mat_file,start_entry*type_byte_size,SEEK_CUR);
    fwrite(data_ptr,type_byte_size,entries,mat_file);
  }
  fclose(mat_file);
  return _SUCCESS_;
}

int mat_read_data(char *filename, 
		  char *matrix_name, 
		  void **data_ptr,
		  int *cols,
		  int *rows,
		  int *data_type){
  uint32_t array_size, name_size, name_size_padded;
  int32_t int32_tmp;
  int type_byte_size;
  int start_data_element;
  FILE *mat_file;
  int func_return, read_status;
  int entries;
  //Open file for reading:
  mat_file = fopen(filename,"rb");
  
  //If matrix name exists in file, get the beginning of the matrix data element:
  func_return = mat_goto_matrix(mat_file,
				matrix_name,
				&start_data_element);
  if (func_return == _FAILURE_){
    fclose(mat_file);
    return _FAILURE_;
  }
  name_size = strlen(matrix_name);
  if ((name_size%8)==0)
    name_size_padded = name_size;
  else
    name_size_padded = name_size + 8-name_size%8;
  //Get dimension:
  fseek(mat_file,start_data_element+4*8,SEEK_SET);
  read_status = fread(&int32_tmp,4,1,mat_file);
  *rows = int32_tmp;
  read_status = fread(&int32_tmp,4,1,mat_file);
  *cols = int32_tmp;
  //Goto the beginning of the array:
  fseek(mat_file,8+name_size_padded,SEEK_CUR);
  //Here is the beginning of the array. Read data type and array size:
  read_status = fread(data_type,4,1,mat_file);
  read_status = fread(&array_size,4,1,mat_file);
  //Determine byte requirements for data-type:
  switch(*data_type){
  case miSINGLE:
  case miINT32:
    type_byte_size = 4;
    break;
  case miDOUBLE:
    type_byte_size = 8;
    break;
  default:
    type_byte_size = 8;
    printf("Warning: type not recognized!\n");
    break;
  }
  //Allocate pointer:
  entries = (*rows)*(*cols);
  *data_ptr = malloc(type_byte_size*entries);
  //Read data:
  read_status = fread(*data_ptr,type_byte_size,entries,mat_file);
  fclose(mat_file);
  return _SUCCESS_;
}
