/*
 * libmesh_apply_solution_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef LIBMESH_APPLY_SOLUTION_INPUT_PARSER_H_
#define LIBMESH_APPLY_SOLUTION_INPUT_PARSER_H_

#include "carl_headers.h"
#include "ext_solver_libmesh_enums.h"

namespace carl
{

struct libmesh_apply_solution_input_params {
	std::string input_vector; 	///< Path to the input vector containing the displacements
	std::string input_mesh;     ///< Path to the mesh that will be deformed
	std::string physical_params_file;	  ///< Physical parameters.
	std::string output_mesh;    ///< Path to the output mesh
};

/**	\brief Parser function for mesh deformation (apply solution) programs
 *	
 *	Required parameters:
 *	  - `InputVector` or `--input-vec`: path to the input vector containing the displacements.
 *    - `InputMesh` or `--input-mesh`: path to the mesh that will be deformed.
 *    - `PhysicalParameters` or `--physical-params`: physical parameters file.
 *    - `OutputMesh` or `--output-mesh`: path to the output mesh.
 */
void get_input_params(GetPot& field_parser,
		libmesh_apply_solution_input_params& input_params);
}

#endif /* LIBMESH_APPLY_SOLUTION_INPUT_PARSER_H_ */
