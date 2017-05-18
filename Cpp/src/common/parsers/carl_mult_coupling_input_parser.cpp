/*
 * carl_mult_coupling_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "carl_mult_coupling_input_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
		carl_mult_coupling_input_params& input_params) {

	if (field_parser.search(2, "CouplMatrix", "-mC")) {
		input_params.coupl_matrix_file = field_parser.next(
				input_params.coupl_matrix_file);
	} else {
		homemade_error_msg("Missing the coupling matrix file!");
	}

	if (field_parser.search(2, "InputVector", "-vI")) {
		input_params.input_vec_file = field_parser.next(
				input_params.input_vec_file);
	} else {
		homemade_error_msg("Missing the input vector file!");
	}

	if (field_parser.search(2, "OutputBase", "-o")) {
		input_params.output_base = field_parser.next(
				input_params.output_base);
	} else {
		homemade_error_msg("Missing the output filename base!");
	}
};
};
