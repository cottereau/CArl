/*
 * carl_mult_coupling_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef CARL_MULT_COUPLING_INPUT_PARSER_H_
#define CARL_MULT_COUPLING_INPUT_PARSER_H_

#include "carl_headers.h"

namespace carl
{
struct carl_mult_coupling_input_params {
	std::string coupl_matrix_file;	///< Path to the coupling matrix file.
	std::string input_vec_file;	///< Path to the input vector file.
	std::string output_base; 	///< Output filename base.
};

/**	\brief Parser function for the coupling matrix multiplication operation.
 *	
 *	Required parameters:
 *	  - `CouplMatrix` or `-mC`: path to the coupling matrix file.
 *    - `InputVector` or `-vI`: path to the input vector file.
 *    - `OutputBase` or `-o`: output filename base.
 *
 */
void get_input_params(GetPot& field_parser,
		carl_mult_coupling_input_params& input_params);
};
#endif /* CARL_MULT_COUPLING_INPUT_PARSER_H_ */
