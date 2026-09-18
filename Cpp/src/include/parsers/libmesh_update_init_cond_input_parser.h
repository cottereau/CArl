/*
 * libmesh_update_init_cond_input_parser.h
 *
 *  Created on: August 5, 2026
 *      Author: Romain Ruyssen
 */

#ifndef LIBMESH_UPDATE_INIT_COND_INPUT_PARSER_H_
#define LIBMESH_UPDATE_INIT_COND_INPUT_PARSER_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"

struct libmesh_update_init_cond_input_params {
	std::string time_mesh_file;		      ///< Path to the time mesh file.
	std::string space_mesh_file;		  ///< Path to the file containing the space mesh.
    std::string newmark_params_file;	  ///< Path to the Newmark parameters file.
	std::string solution_file;            ///< Path to the solution file.
    std::string init_cond_disp_file;      ///< Path to the initial displacement file.
	std::string init_cond_vel_file;       ///< Path to the initial velocity file.
	std::string init_cond_acc_file;       ///< Path to the initial acceleration file.
};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *	  - `Mesh` : path to the mesh.
*/
void get_input_params(GetPot& field_parser,
		libmesh_update_init_cond_input_params& input_params);
#endif /* LIBMESH_UPDATE_INIT_COND_INPUT_PARSER_H_ */
