/*
 * \file libmesh_apply_solution_dyn_input_parser.cpp
 *
 *  Created on: Feb 23, 2022
 *      Author: Chensheng Luo
 *
 * \brief **DYN-DI/DYN-CG** Parser function for mesh deformation (apply solution) programs
 */
#include "libmesh_apply_solution_dyn_input_parser.h"

namespace carl
{
void get_input_params(GetPot& field_parser,
		libmesh_apply_solution_dyn_input_params& input_params) {

	if (field_parser.search(2, "InputVectorFolder", "--input-vec" )) {
		input_params.input_vector_folder = field_parser.next(
				input_params.input_vector_folder);
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
		input_params.output_mesh_folder = field_parser.next(
				input_params.output_mesh_folder);
	} else {
		homemade_error_msg("Missing the output mesh path!");
	}

	if (field_parser.search(2, "TotalLoopTimes", "--total-loop")) {
    	input_params.total_loop_times = field_parser.next(input_params.total_loop_times);
  	} else {
    	input_params.total_loop_times = 10;
  	}

  	if (field_parser.search(2, "StepLoopTimes", "--step-loop")) {
    	input_params.step_loop_times = field_parser.next(input_params.step_loop_times);
  	} else {
    	input_params.step_loop_times = 1;
  	}

    if (field_parser.search(1,"NewmarkParameters")){
        std::string filename;
        filename = field_parser.next(filename);
        GetPot newmark_parser;
        newmark_parser.parse_input_file(filename, "#", "\n", " \t\n");
        carl::get_newmark_params(newmark_parser, input_params.Newmark);
      } else{
        homemade_error_msg("[CArl Parameters]Missing Newmark parameters file!");
      }

};

};
