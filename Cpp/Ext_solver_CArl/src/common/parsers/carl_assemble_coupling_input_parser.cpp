/*
 * carl_assemble_coupling_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "carl_assemble_coupling_input_parser.h"

namespace carl
{

void get_assemble_coupling_input_params(GetPot& field_parser,
		coupling_assemble_coupling_input_params& input_params) {

	// Set mesh files
	if (field_parser.search(3, "--meshA", "-mA", "MeshA")) {
		input_params.mesh_BIG_file = field_parser.next(
				input_params.mesh_BIG_file);
	} else {
		homemade_error_msg("Missing the A mesh file!");
	}

	if (field_parser.search(3, "--meshB", "-mB", "MeshB")) {
		input_params.mesh_micro_file = field_parser.next(
				input_params.mesh_micro_file);
	} else {
		homemade_error_msg("Missing the B mesh file!");
	}

	if (field_parser.search(3, "--meshI", "-mI", "InterBase")) {
		input_params.common_inter_file = field_parser.next(
				input_params.common_inter_file);
	} else {
		homemade_error_msg("Missing the path to the intersection files!");
	}

	// Set coupling parameters
	if( field_parser.search(2, "--ce","CouplingWidth") )
	{
		input_params.coupling_width = field_parser.next(input_params.coupling_width);
	} else {
		homemade_error_msg("Missing the coupling region width!");
	}

	if( field_parser.search(2, "--ck","CouplingRigidity") )
	{
		input_params.coupling_rigidity = field_parser.next(input_params.coupling_rigidity);
	} else {
		homemade_error_msg("Missing the coupling rigidity!");
	}

	// Output
	if (field_parser.search(3, "--output", "-mO", "OutputBase"))
	{
		input_params.output_base = field_parser.next(
			input_params.output_base);
	} else {
		input_params.output_base = "coupling_matrix";
	}

	if (field_parser.search(3, "--meshAR", "-mAR", "Mesh_A_Restriction")) {
		input_params.mesh_restrict_BIG_file = field_parser.next(
				input_params.mesh_restrict_BIG_file);
	} else {
		input_params.mesh_restrict_BIG_file = input_params.common_inter_file + "_A_restriction.msh";
	}

	if (field_parser.search(3, "--meshBR", "-mBR", "Mesh_B_Restriction")) {
		input_params.mesh_restrict_micro_file = field_parser.next(
				input_params.mesh_restrict_micro_file);
	} else {
		input_params.mesh_restrict_micro_file = input_params.common_inter_file + "_B_restriction.msh";
	}

	// Set the equivalence tables
	if (field_parser.search(2, "--tableRA", "Mesh_A_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_BIG_file = field_parser.next(
				input_params.equivalence_table_restrict_BIG_file);
	} else {
		input_params.equivalence_table_restrict_BIG_file = input_params.common_inter_file + "_A_restriction_restrict.dat";
	}

	if (field_parser.search(2, "--tableRB", "Mesh_B_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_micro_file = field_parser.next(
				input_params.equivalence_table_restrict_micro_file);
	} else {
		input_params.equivalence_table_restrict_micro_file = input_params.common_inter_file + "_B_restriction_restrict.dat";
	}

	// Set the mediator mesh
	input_params.mediator_type = carl::MediatorType::USE_MACRO;
	input_params.mesh_mediator_file = input_params.mesh_restrict_BIG_file;

	std::string mediator_type;
	if (field_parser.search(1,"MediatorMesh"))
	{
		mediator_type = field_parser.next(mediator_type);

		if(mediator_type == "UseRestricted_A")
		{
			input_params.mediator_type = carl::MediatorType::USE_MACRO;
			input_params.mesh_mediator_file = input_params.mesh_restrict_BIG_file;
		}
		else if(mediator_type == "UseRestricted_B")
		{
			input_params.mediator_type = carl::MediatorType::USE_MICRO;
			input_params.mesh_mediator_file = input_params.mesh_restrict_micro_file;
		}
	}
};

};