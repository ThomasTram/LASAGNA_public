/** @file input.h Documented includes for input module */

#ifndef __INPUT__
#define __INPUT__

#include "common.h"
#include "parser.h"
#include "qke_equations.h"

/* macro for reading parameter values with routines from the parser */
#define lasagna_read_double(name,destination)				\
  do {									\
    lasagna_call(parser_read_double(pfc,name,&param1,&flag1,errmsg),      \
	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_)						\
      destination = param1;						\
  } while(0);


#define lasagna_read_int(name,destination)				\
  do {									\
    lasagna_call(parser_read_int(pfc,name,&int1,&flag1,errmsg),		\
 	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_)						\
      destination = int1;						\
  } while(0);

#define lasagna_read_string(name,destination)				\
  do {									\
    lasagna_call(parser_read_string(pfc,name,&string1,&flag1,errmsg),	\
 	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_)						\
      strcpy(destination,string1);					\
  } while(0);

#define lasagna_read_double_one_of_two(name1,name2,destination)		\
  do {									\
    lasagna_call(parser_read_double(pfc,name1,&param1,&flag1,errmsg),	\
	       errmsg,							\
	       errmsg);							\
    lasagna_call(parser_read_double(pfc,name2,&param2,&flag2,errmsg),	\
	       errmsg,							\
	       errmsg);							\
    lasagna_test((flag1 == _TRUE_) && (flag2 == _TRUE_),			\
	       errmsg,							\
	       "In input file, you can only enter one of %s, %s, choose one", \
	       name1,name2);						\
    if (flag1 == _TRUE_)						\
      destination = param1;						\
    if (flag2 == _TRUE_)						\
      destination = param2;						\
  } while(0);

#define lasagna_at_least_two_of_three(a,b,c)		\
  ((a == _TRUE_) && (b == _TRUE_)) ||		\
  ((a == _TRUE_) && (c == _TRUE_)) ||		\
  ((b == _TRUE_) && (c == _TRUE_))

#define lasagna_none_of_three(a,b,c)				\
  (a == _FALSE_) && (b == _FALSE_) && (c == _FALSE_)

/* macro for reading parameter values with routines from the parser */
#define lasagna_read_list_of_doubles_or_default(name,destination,default,siz)	\
  do {									\
    lasagna_call(parser_read_list_of_doubles(pfc,name,			\
	&entries_read,&(destination),&flag1,errmsg),			\
	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_){						\
        lasagna_test(entries_read != siz,errmsg,			\
             "Number of entries in %s, %d, does not match number of indistinguishable ncdm species, %d.", \
		name,entries_read,siz);				\
    }else{								\
	lasagna_alloc(destination,siz*sizeof(double),errmsg);		\
	for(n=0; n<siz; n++) destination[n] = default;		\
    }									\
  } while(0);

#define lasagna_read_list_of_integers_or_default(name,destination,default,siz) \
  do {									\
    lasagna_call(parser_read_list_of_integers(pfc,name,			\
	&entries_read,&(destination),&flag1,errmsg),			\
	       errmsg,							\
	       errmsg);							\
    if (flag1 == _TRUE_){						\
        lasagna_test(entries_read != siz,errmsg,			\
             "Number of entries in %s, %d, does not match number of indistinguishable ncdm species, %d.", \
		name,entries_read,siz);				\
    }else{								\
	lasagna_alloc(destination,siz*sizeof(int),errmsg);		\
	for(n=0; n<siz; n++) destination[n] = default;		\
    }									\
  } while(0);

#define lasagna_read_list_of_doubles(name,destination,siz)			\
  do {									\
    lasagna_call(parser_read_list_of_doubles(pfc,name,			\
	&entries_read,&(destination),&flag1,errmsg),			\
	       errmsg,							\
	       errmsg);							\
    lasagna_test(flag1 == _FALSE_,errmsg,					\
	"Entry %s is required but not found!",name)			\
        lasagna_test(entries_read != siz,errmsg,			\
             "Number of entries in %s, %d, does not match number of indistinguishable ncdm species, %d.", \
		name,entries_read,siz);				\
  } while(0);

#define lasagna_read_list_of_integers(name,destination,siz)			\
  do {									\
    lasagna_call(parser_read_list_of_integers(pfc,name,			\
	&entries_read,&(destination),&flag1,errmsg),			\
	       errmsg,							\
	       errmsg);							\
    lasagna_test(flag1 == _FALSE_,errmsg,					\
	"Entry %s is required but not found!",name)			\
        lasagna_test(entries_read != siz,errmsg,			\
             "Number of entries in %s, %d, does not match number of indistinguishable ncdm species, %d.", \
		name,entries_read,siz);				\
  } while(0);


/**************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int input_init_from_arguments(
				int argc, 
				char **argv,
				qke_param * pqke,
				ErrorMsg errmsg);

  int input_init(
		 struct file_content * pfc,
		 qke_param *pqke,
		 ErrorMsg errmsg
		 );

  int input_default_params(
			   qke_param * pqke
			   );


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
