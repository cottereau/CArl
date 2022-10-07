/*
 * \file libmesh_apply_solution_dyn_input_parser.h
 *
 *  Created on: Feb 23, 2022
 *      Author: Chensheng Luo
 */

#ifndef LIBMESH_APPLY_SOLUTION_DYN_INPUT_PARSER_H_
#define LIBMESH_APPLY_SOLUTION_DYN_INPUT_PARSER_H_

#include "carl_headers.h"
#include "ext_solver_libmesh_enums.h"
#include "newmark_param_parser.h"

namespace carl
{

struct libmesh_apply_solution_dyn_input_params {
	std::string input_vector_folder; 	///< Path to the input vector folder containing the displacements
	std::string input_mesh;     ///< Path to the mesh that will be deformed
	std::string physical_params_file;	  ///< Physical parameters.
	std::string output_mesh_folder;    ///< Path to the output mesh
    int total_loop_times;       ///< Total loop times=(max(index))
    int step_loop_times;       ///< Step loop time=difference between each loop step
    carl::NewmarkParams Newmark;    ///< Newmark parameter of the model

};

/**	\brief **DYN-DI/DYN-CG** Parser function for mesh deformation (apply solution) programs
 *	
 *	Required parameters:
 *	  - `InputVectorFolder` or `--input-vec`: path to the input vector folder containing the displacements.
 *    - `InputMesh` or `--input-mesh`: path to the mesh that will be deformed.
 *    - `PhysicalParameters` or `--physical-params`: physical parameters file.
 *    - `OutputMesh` or `--output-mesh`: path to the output mesh.
 *    - `TotalLoopTimes` or `--total-loop`: total loop times.
 *    - `StepLoopTimes` or `--step-loop`: step loop times.
 *    - `TimeStep` or `--time-step`: time step of B
 * 
 */
void get_input_params(GetPot& field_parser,
		libmesh_apply_solution_dyn_input_params& input_params);
}

#endif /* LIBMESH_APPLY_SOLUTION_DYN_INPUT_PARSER_H_ */
