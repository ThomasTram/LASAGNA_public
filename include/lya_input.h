/** @file lya_input.h Documented includes for input module */

#ifndef __LYA_INPUT__
#define __LYA_INPUT__

#include "common.h"
#include "parser.h"
#include "lya_equations.h"
#include "input.h"

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int lya_input_init_from_arguments(
				    int argc, 
				    char **argv,
				    lya_param * pqke,
				    ErrorMsg errmsg);

  int lya_input_init(
		     struct file_content * pfc,
		     lya_param *pqke,
		     ErrorMsg errmsg
		     );

  int lya_input_default_params(
			       lya_param * pqke
			       );

  int lya_input(
		int argc, 
		char **argv,
		lya_param *plya,
		ErrorMsg errmsg
		);


#ifdef __cplusplus
}
#endif

#endif
