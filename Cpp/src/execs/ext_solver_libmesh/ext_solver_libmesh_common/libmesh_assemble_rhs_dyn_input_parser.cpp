/*
 * libmesh_assemble_rhs_dyn_input_parser.h
 *
 *  Created on: July 27, 2026
 *      Author: Romain Ruyssen
 */

#include "libmesh_assemble_rhs_dyn_input_parser.h"
void get_input_params(GetPot& field_parser,
		libmesh_assemble_rhs_dyn_input_params& input_params) {

	if (field_parser.search(1, "ScratchFolderPath")) {
		input_params.scratch_folder_path = field_parser.next(
				input_params.scratch_folder_path);
		std::cout << input_params.scratch_folder_path << std::endl;
	} else {
		homemade_error_msg("Missing the scratch folder path!");
	}		

	if (field_parser.search(1, "TimeMeshFile")) {
		input_params.time_mesh_file = field_parser.next(
				input_params.time_mesh_file);
		std::cout << input_params.time_mesh_file << std::endl;
	} else {
		homemade_error_msg("Missing the time mesh file !");
	}

	if (field_parser.search(1, "Mesh")) {
		input_params.space_mesh_file = field_parser.next(
				input_params.space_mesh_file);
		std::cout << input_params.space_mesh_file << std::endl;
	} else {
		homemade_error_msg("Missing the space mesh file !");
	}

    if (field_parser.search(1, "NewmarkParameters")) {
        input_params.newmark_params_file = field_parser.next(
                input_params.newmark_params_file);
        std::cout << input_params.newmark_params_file << std::endl;
    } else {
        homemade_error_msg("Missing the Newmark parameters file !");
    }

    //-------------------------------------------------------
	// Fichiers contenant les différentes matrices du domaine
	if (field_parser.search(1, "SysMatrix")) {
		input_params.domain_matrix_file = field_parser.next(
				input_params.domain_matrix_file);
		std::cout << input_params.domain_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_M")) {
		input_params.domain_M_matrix_file = field_parser.next(
				input_params.domain_M_matrix_file);
		std::cout << input_params.domain_M_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the M matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_K")) {
		input_params.domain_K_matrix_file = field_parser.next(
				input_params.domain_K_matrix_file);
		std::cout << input_params.domain_K_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the K matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_C")) {
		input_params.domain_C_matrix_file = field_parser.next(
				input_params.domain_C_matrix_file);
		std::cout << input_params.domain_C_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the C matrix file !");
	}

};
