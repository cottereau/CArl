/*
 * carl_assemble_coupling_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef CARL_ASSEMBLE_COUPLING_INPUT_PARSER_H_
#define CARL_ASSEMBLE_COUPLING_INPUT_PARSER_H_

#include "carl_headers.h"

namespace carl
{
/**	\brief Structure containing the parameters for the construction of the coupling matrices.
 *	
 *		Details on the parameters setup are found in the documentation of carl::get_assemble_coupling_input_params(GetPot& field_parser,
 *		coupling_assemble_coupling_input_params& input_params).
 */
struct coupling_assemble_coupling_input_params {
	std::string mesh_BIG_file;		///< Path to the macro (BIG) system mesh.
	std::string mesh_micro_file;	///< Path to the micro system mesh.

	std::string mesh_restrict_BIG_file;		///< Path to the restricted macro (BIG) system mesh.
	std::string mesh_restrict_micro_file;	///< Path to the restricted micro system mesh.

	std::string mesh_mediator_file;		///< Path to the mediator system mesh.
	std::string common_inter_file;		///< Common path to the intersection meshes and tables.

	/// Equivalence table between the macro (BIG) system mesh and its restriction mesh
	std::string equivalence_table_restrict_BIG_file;
	/// Equivalence table between the micro system mesh and its restriction mesh
	std::string equivalence_table_restrict_micro_file;

	carl::MediatorType   mediator_type;		///< Use which mesh as the mediator mesh? *Values*: carl::MediatorType::USE_MACRO, carl::MediatorType::USE_MICRO or carl::MediatorType::USE_EXTERNAL (LAST CASE NOT IMPLEMENTED YET!).
	bool b_UseMesh_micro_AsMediator;	///< Use the micro system's restricted mesh as the mediator?

	double coupling_width;		///< Width  of the coupling region.
	double coupling_rigidity; 	///< Rigidity used for the coupling matrix.

	std::string output_base; 	///< Output filename base.

};

/**	\brief Parser function for the construction of the coupling matrices.
 *	
 *	Required parameters:
 *  + System and intersection meshes:
 *	  - `MeshA`, `-mA` or `--meshA` : path to the mesh A.
 *	  - `MeshB`, `-mB` or `--meshB` : path to the mesh B.
 *	  - `InterBase`, `-mI` or `--meshI` : common path to the intersection meshes and tables.
 *  + Coupling parameters:
 *    - `CouplingWidth` or `--ce` : width  of the coupling region (same unit as the meshes, \f$e\f$ in the \f$L_2\f$ coupling term).
 *    - `CouplingRigidity` or `--ck` : rigidity used for the coupling matrix (in MPa, if mm was used for the meshes, \f$\kappa\f$ in both the \f$L_2\f$ and \f$H_1\f$ terms).
 *
 *  Optional parameters:
 *  + Output:
 *    - `OutputBase` or `--output` : base of the output files (including folders). *Default*: `coupling_matrix`.
 *  + Restriction meshes and tables:
 *	  - `Mesh_A_Restriction`, `-mAR` or `--meshAR` : path to the restricted mesh A (formed by elements of the mesh A intersecting the coupling region). *Default*: `[InterBase]_A_restriction.msh`.
 *	  - `Mesh_B_Restriction`, `-mBR` or `--meshBR` : path to the restricted mesh B (formed by elements of the mesh A intersecting the coupling region). *Default*: `[InterBase]_B_restriction.msh`.
 *    - `Mesh_A_RestrictionEquivalenceTable` or `--tableRA` : path to the equivalence table between the mesh A and its restriction. *Default*: `[InterBase]_A_restriction_restrict.dat`.
 *    - `Mesh_B_RestrictionEquivalenceTable` or `--tableRB` : path to the equivalence table between the mesh B and its restriction. *Default*: `[InterBase]_B_restriction_restrict.dat`.
 *  + Mediator mesh:
 *    - `MediatorMesh` : choice of the mediator mesh. *Values*: `UseRestricted_A` or `UseRestricted_B`. *Default*: `UseRestricted_A`.
 */
void get_assemble_coupling_input_params(GetPot& field_parser,
		coupling_assemble_coupling_input_params& input_params) {

	// Set mesh files
	if (field_parser.search(3, "--meshA", "-mA", "MeshA")) {
		input_params.mesh_BIG_file = field_parser.next(
				input_params.mesh_BIG_file);
		std::cout << input_params.mesh_BIG_file << std::endl;
	} else {
		homemade_error_msg("Missing the A mesh file!");
	}

	if (field_parser.search(3, "--meshB", "-mB", "MeshB")) {
		input_params.mesh_micro_file = field_parser.next(
				input_params.mesh_micro_file);
		std::cout << input_params.mesh_micro_file << std::endl;
	} else {
		homemade_error_msg("Missing the B mesh file!");
	}

	if (field_parser.search(3, "--meshI", "-mI", "InterBase")) {
		input_params.common_inter_file = field_parser.next(
				input_params.common_inter_file);
		std::cout << input_params.common_inter_file << std::endl;
	} else {
		homemade_error_msg("Missing the path to the intersection files!");
	}

	// Set coupling parameters
	if( field_parser.search(2, "--ce","CouplingWidth") )
	{
		input_params.coupling_width = field_parser.next(input_params.coupling_width);
		std::cout << input_params.coupling_width << std::endl;
	} else {
		homemade_error_msg("Missing the coupling region width!");
	}

	if( field_parser.search(2, "--ck","CouplingRigidity") )
	{
		input_params.coupling_rigidity = field_parser.next(input_params.coupling_rigidity);
		std::cout << input_params.coupling_rigidity << std::endl;
	} else {
		homemade_error_msg("Missing the coupling rigidity!");
	}

	// Output
	if (field_parser.search(3, "--output", "-mO", "OutputBase"))
	{
		input_params.output_base = field_parser.next(
			input_params.output_base);
		std::cout << input_params.output_base << std::endl;
	} else {
		input_params.output_base = "coupling_matrix";
	}

	if (field_parser.search(3, "--meshAR", "-mAR", "Mesh_A_Restriction")) {
		input_params.mesh_restrict_BIG_file = field_parser.next(
				input_params.mesh_restrict_BIG_file);
	} else {
		input_params.mesh_restrict_BIG_file = input_params.common_inter_file + "_A_restriction.msh";
	}
		std::cout << input_params.mesh_restrict_BIG_file << std::endl;

	if (field_parser.search(3, "--meshBR", "-mBR", "Mesh_B_Restriction")) {
		input_params.mesh_restrict_micro_file = field_parser.next(
				input_params.mesh_restrict_micro_file);
	} else {
		input_params.mesh_restrict_micro_file = input_params.common_inter_file + "_B_restriction.msh";
	}
		std::cout << input_params.mesh_restrict_micro_file << std::endl;

	// Set the equivalence tables
	if (field_parser.search(2, "--tableRA", "Mesh_A_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_BIG_file = field_parser.next(
				input_params.equivalence_table_restrict_BIG_file);
	} else {
		input_params.equivalence_table_restrict_BIG_file = input_params.common_inter_file + "_A_restriction_restrict.dat";
	}
		std::cout << input_params.equivalence_table_restrict_BIG_file << std::endl;

	if (field_parser.search(2, "--tableRB", "Mesh_B_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_micro_file = field_parser.next(
				input_params.equivalence_table_restrict_micro_file);
	} else {
		input_params.equivalence_table_restrict_micro_file = input_params.common_inter_file + "_B_restriction_restrict.dat";
	}
		std::cout << input_params.equivalence_table_restrict_micro_file << std::endl;

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
		std::cout << input_params.mediator_type << std::endl;
		std::cout << input_params.mesh_mediator_file << std::endl;
};

};
#endif /* CARL_ASSEMBLE_COUPLING_INPUT_PARSER_H_ */
