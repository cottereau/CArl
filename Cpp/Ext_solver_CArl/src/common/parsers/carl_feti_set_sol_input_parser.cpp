/*
 * carl_feti_set_sol_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */


#include "carl_feti_set_sol_input_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
		feti_set_sol_params& input_params) {

	if (field_parser.search(1, "ScratchFolderPath")) {
		input_params.scratch_folder_path = field_parser.next(
				input_params.scratch_folder_path);
	} else {
		homemade_error_msg("Missing the external scratch folder path!");
	}

	if (field_parser.search(1,"UseRigidBodyModesB"))
	{
		input_params.bUseRigidBodyModes = true;
	} else {
		input_params.bUseRigidBodyModes = false;
	}

	if (field_parser.search(1,"OutputBase")) {
		input_params.output_base = field_parser.next(
				input_params.output_base);
	} else {
		homemade_error_msg("Missing the output filename base!");
	}
};

};
