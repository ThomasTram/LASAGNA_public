#ifndef _IO_
#define _IO_

#include "common.h"
#include <stdint.h>

#define MAT_NO_HANDLE -1
/** MATLAB array classes: */
#define mxSINGLE_CLASS 7
#define mxDOUBLE_CLASS 6
#define mxINT32_CLASS 12
/** Storage data types: */
#define miINT8 1
#define miUINT8 2
#define miINT32 5
#define miUINT32 6
#define miSINGLE 7
#define miDOUBLE 9
#define miMATRIX 14

/** Macro for updating output file: */
#define mat_write_fast(ptr,handle,size_element,entries)       		\
  do {									\
    fseek(mat_file,handle,SEEK_SET);					\
    fwrite(ptr,size_element,entries,mat_file);				\
    handle += size_element*entries;					\
  } while(0);

/** Macro for reading */
#define lasagna_fread(destination_ptr,					\
		      size,						\
		      count,						\
		      stream)						\
  do {									\
    read_status = fread(destination_ptr, size, count, stream);		\
    if (read_status != count) {						\
      if (feof(mat_file)==1)						\
	printf("Error: %s(L:%d). Read failed, End-of-File reached.\n",	\
	       __func__,__LINE__);					\
      else								\
	printf("Error: %s(L:%d). Read failed, unknown reason.\n",	\
	       __func__,__LINE__);					\
      return _FAILURE_;							\
    }									\
  } while(0);								\
  
/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
  int mat_create_file(char *filename);
  int mat_add_matrix(char *filename, 
		     char *matrix_name,
		     int data_type,
		     int cols, 
		     int rows,
		     int *handle);
  int mat_write_data(char *filename, 
		     char *matrix_name, 
		     void *data_ptr, 
		     int start_entry,
		     int entries);
  int mat_read_data(char *filename, 
		    char *matrix_name, 
		    void **data_ptr,
		    int *cols,
		    int *rows,
		    int *data_type);
  int mat_goto_matrix(FILE *mat_file,
		      char *matrix_name,
		      int *beginning_of_matrix);


  /** Remember that data is stored in subsequent columns. 
      So mat_write_double is easy for overwriting a column
      in the existing matrix, but overwriting a row is a non-trivial
      operation.
  */

#ifdef __cplusplus
}
#endif


#endif
