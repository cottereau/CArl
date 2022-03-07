/*
 * libmesh_apply_solution_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */
#include "libmesh_apply_solution_input_parser.h"

namespace carl
{
void get_input_params(GetPot& field_parser,
		libmesh_apply_solution_input_params& input_params) {

	if (field_parser.search(2, "InputVector", "--input-vec" )) {
		input_params.input_vector = field_parser.next(
				input_params.input_vector);
	} else {
		homemade_error_msg("Missing the displacement field vector!");
	}

	if (field_parser.search(2, "InputMesh", "--input-mesh")) {
		input_params.input_mesh = field_parser.next(
				input_params.input_mesh);
	} else {
		homemade_error_msg("Missing the input mesh!");
	}

	if (field_parser.search(2, "PhysicalParameters", "--physical-params")) {
		input_params.physical_params_file = field_parser.next(
				input_params.physical_params_file);
	} else {
		homemade_error_msg("Missing the physical params file!");
	}

	if (field_parser.search(2, "OutputMesh", "--output-mesh")) {
		input_params.output_mesh = field_parser.next(
				input_params.output_mesh);
	} else {
		homemade_error_msg("Missing the output mesh path!");
	}
};

};
