/*
 * carl_feti_set_sol_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef CARL_FETI_SET_SOL_INPUT_PARSER_H_
#define CARL_FETI_SET_SOL_INPUT_PARSER_H_

#include "carl_headers.h"

namespace carl
{
/// Structure containing the parameters for the setup initialization of the FETI solver.
struct feti_set_sol_params {
	// --- Parameters used directly by the CArl_FETI_set_sol program (some are also used by the other CArl_FETI programs)

	// Path to "scratch" folder
	std::string scratch_folder_path;	///< Path to the folder which will be used to save the temporary files during the solve operation

	// Rigid body mode options for the micro system
	bool bUseRigidBodyModes;			///< [RB] Use the rigid body modes for the micro system?

	// Path to the "final" output base
	std::string output_base;			///< Base of the final output path.
};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *	  - `ScratchFolderPath` : path to the folder where the temporary files used by the coupled solver will be saved.
 *	  - `OutputBase` : base of the final output files.
 *
 *  Boolean flags:
 *    - `UseRigidBodyModesB` : use the rigid body modes for system B.
 */
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
#endif /* CARL_FETI_SET_SOL_INPUT_PARSER_H_ */
