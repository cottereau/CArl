/*
 * libmesh_assemble_system_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef LIBMESH_ASSEMBLE_SYSTEM_INPUT_PARSER_H_
#define LIBMESH_ASSEMBLE_SYSTEM_INPUT_PARSER_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"

struct libmesh_assemble_input_params {
	std::string mesh_file;				  ///< Path to the system mesh.
	std::string physical_params_file;	  ///< Physical parameters.
	WeightFunctionSystemType system_type; ///< Indicates if the system to be assembled is a micro or a macro system (used to choose the proper weight function).

	std::string mesh_weight_file;		///< Path to the mesh containing the weight region indices.
	std::string weight_domain_idx_file; ///< Path to the file identifying the weight function regions.

	std::string output_base; 	///< Output filename base.
};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *	  - `Mesh` : path to the mesh.
 *    - `PhysicalParameters` : physical parameters.
 *    - `SystemType` : parameter used to tell the assembler which weight functions must be used. *Values*: `Micro` or `Macro`.
 *	  - `MeshWeight` : path to the mesh defining the domains of the Arlequin weight parameters.
 *    - `WeightIndexes` : path to the indices of the domains of the Arlequin weight parameters.
 *    - `OutputBase` or `--output` : base of the output files (including folders). *Default*: `test_system`.
 */
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
	if (field_parser.search(2, "--output", "OutputBase"))
	{
		input_params.output_base = field_parser.next(
			input_params.output_base);
	} else {
		input_params.output_base = "test_system";
	}
};
#endif /* LIBMESH_ASSEMBLE_SYSTEM_INPUT_PARSER_H_ */
