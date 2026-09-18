/*
 * libmesh_init_dyn_input_parser.h
 *
 *  Created on: Apr 3, 2026
 *      Author: Romain Ruyssen
 */

#include "libmesh_init_dyn_input_parser.h"
void get_input_params(GetPot& field_parser,
		libmesh_init_dyn_input_params& input_params) {

	// Set mesh files
	if (field_parser.search(1, "Mesh")) {
		input_params.mesh_file = field_parser.next(
				input_params.mesh_file);
	} else {
		homemade_error_msg("Missing the system mesh file!");
	}

	// Set constant parameters
	if ( field_parser.search(1, "PhysicalParameters") )
	{
		input_params.physical_params_file = field_parser.next(input_params.physical_params_file);
	}
	else
	{
		homemade_error_msg("Missing the physical parameters file!");
	}

	// Set constant parameters
	if ( field_parser.search(1, "TimeDiscretization") )
	{
		input_params.time_discretization_file = field_parser.next(input_params.time_discretization_file);
	}
	else
	{
		homemade_error_msg("Missing the Time Discretization parameters file!");
	}

	// Set constant parameters
	if ( field_parser.search(1, "NewmarkParameters") )
	{
		input_params.newmark_params_file = field_parser.next(input_params.newmark_params_file);
	}
	else
	{
		homemade_error_msg("Missing the Newmark parameters file!");
	}

	// Set weight function
	std::string sys_type;
	if ( field_parser.search(1, "SystemType") )
	{
		sys_type = field_parser.next(sys_type);
		if(sys_type == "Macro" || sys_type == "MACRO" || sys_type == "macro")
			input_params.system_type = WeightFunctionSystemType::MACRO;
		else if(sys_type == "Micro" || sys_type == "MICRO" || sys_type == "micro")
			input_params.system_type = WeightFunctionSystemType::MICRO;
		else if(sys_type == "NoWeight" || sys_type == "NOWEIGHT" || sys_type == "noweight")
		{
			input_params.system_type = WeightFunctionSystemType::NO_WEIGHT;
			std::cout << " >> Warning: Will not use the weight parameters!" << std::endl;
		}
		else
			homemade_error_msg("Invalid system type (must be either Macro, Micro or NoWeight)!");
	}
	else
	{
		homemade_error_msg("Missing the system type (must be either Macro, Micro or NoWeight)!");
	}

	if ( field_parser.search(1, "MeshWeight") )
	{
		input_params.mesh_weight_file = field_parser.next(input_params.mesh_weight_file);
	}
	else
	{
		homemade_error_msg("Missing the weight mesh file!");
	}

	if( field_parser.search(1, "WeightIndexes") )
	{
		input_params.weight_domain_idx_file = field_parser.next(input_params.weight_domain_idx_file);
	}
	else
	{
		homemade_error_msg("Missing the weight value file!");
	}

	// Output
	if (field_parser.search(2, "--output", "OutputBaseMatrix"))
	{
		input_params.output_base_matrix = field_parser.next(
			input_params.output_base_matrix);
	} else {
		input_params.output_base_matrix = "test_system";
	}

	if (field_parser.search(2, "--output", "OutputBaseInitCond"))
	{
		input_params.output_base_init_cond = field_parser.next(
			input_params.output_base_init_cond);
	} else {
		input_params.output_base_init_cond = "test_system";
	}

	if (field_parser.search(1, "ExportRBVectors")) {
		input_params.bCalculateRBVectors = true;
	} else {
		input_params.bCalculateRBVectors = false;
	}
};
