/*
 * libmesh_assemble_rhs_dyn_input_parser.h
 *
 *  Created on: July 27, 2026
 *      Author: Romain Ruyssen
 */

#ifndef LIBMESH_ASSEMBLE_RHS_DYN_INPUT_PARSER_H_
#define LIBMESH_ASSEMBLE_RHS_DYN_INPUT_PARSER_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"

struct libmesh_assemble_rhs_dyn_input_params {
	
    // Path to "scratch" folder
	std::string scratch_folder_path;	///< Path to the folder which will be used to save the temporary files during the solve operation
	
	// Path to the time mesh file
	std::string time_mesh_file;	        ///< Path to the file containing the time mesh. 

	// Path to the space mesh file
	std::string space_mesh_file;		///< Path to the file containing the space mesh.

    // Path to the newmark parameters file
    std::string newmark_params_file;      ///< Newmark parameters.
	
	// Domains matrices files
	std::string domain_matrix_file;	    ///< Path to the file containing the matrix.
	std::string domain_M_matrix_file;	///< Path to the file containing the matrix.
	std::string domain_K_matrix_file;	///< Path to the file containing the matrix.
	std::string domain_C_matrix_file;	///< Path to the file containing the matrix.

};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *	  - `Mesh` : path to the mesh.
 *    - 
 *
 *  Optional parameter:
 *    
 *  Boolean flags:
 *    
 */
void get_input_params(GetPot& field_parser,
		libmesh_assemble_rhs_dyn_input_params& input_params);
#endif /* LIBMESH_ASSEMBLE_RHS_DYN_INPUT_PARSER_H_ */
