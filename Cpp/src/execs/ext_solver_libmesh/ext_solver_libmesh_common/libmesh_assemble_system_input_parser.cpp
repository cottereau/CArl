/*
 * libmesh_assemble_system_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 *  Added dynamics on: May 31, 2018
 *      Author: Filippo Gatti
 */
// This file contains the subroutines to read the input files
#include "libmesh_assemble_system_input_parser.h"
void get_input_params(GetPot& field_parser,
		libmesh_assemble_input_params& input_params) {

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
        std::cout << " >> INPUT FILE: " << input_params.physical_params_file; 
	}
	else
	{
		homemade_error_msg("Missing the physical parameters file!");
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
	if (field_parser.search(2, "--output", "OutputBase"))
	{
		input_params.output_base = field_parser.next(
			input_params.output_base);
	} else {
		input_params.output_base = "test_system";
	}

	if (field_parser.search(1, "ExportRBVectors")) {
		input_params.bCalculateRBVectors = true;
	} else {
		input_params.bCalculateRBVectors = false;
	}
};
