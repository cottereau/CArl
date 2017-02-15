/*
 * carl_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef CARL_INPUT_PARSER_H_
#define CARL_INPUT_PARSER_H_

#include "carl_headers.h"

namespace carl
{
struct coupling_generation_input_params {
	std::string physical_params_file;

	std::string mesh_BIG_file;
	std::string mesh_micro_file;

	std::string mesh_restrict_BIG_file;
	std::string mesh_restrict_micro_file;

	std::string mesh_mediator_file;
	std::string mesh_inter_file;

	std::string mesh_weight_file;

	std::string equivalence_table_restrict_BIG_file;
	std::string equivalence_table_restrict_micro_file;
	std::string equivalence_table_mediator;

	std::string intersection_table_full;

	std::string weight_domain_idx_file;

	std::string scaling_data_file;

	bool b_UseMesh_BIG_AsMediator;
	bool b_UseMesh_micro_AsMediator;
	bool b_UseMesh_extra_AsMediator;
	bool b_Repartition_micro;
	bool LATIN_b_UseRestartFiles;
	bool LATIN_b_PrintRestartFiles;
	bool b_PrintOutput;
	bool b_ExportScalingData;

	double mean_distance;

	double k_dA;
	double k_dB;
	double k_cA;
	double k_cB;

	double LATIN_eps;
	int LATIN_conv_max;
	double LATIN_relax;

	double coupled_conv_abs;
	double coupled_conv_rel;
	double coupled_div;
	int coupled_iter_max;

	std::string coupled_convergence_output;
	std::string coupled_restart_file_base;

	std::string output_file_BIG;
	std::string output_file_micro;

	std::string solver_type_string;
	carl::CoupledSolverType solver_type;
	std::string CG_precond_type_string;
	carl::BaseCGPrecondType CG_precond_type;
};

void get_input_params(GetPot& field_parser,
		coupling_generation_input_params& input_params) {

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

	if (field_parser.search(3, "--meshAR", "-mAR", "Mesh_A_Restriction")) {
		input_params.mesh_restrict_BIG_file = field_parser.next(
				input_params.mesh_restrict_BIG_file);
	} else {
		homemade_error_msg("Missing the restricted A mesh file!");
	}

	if (field_parser.search(3, "--meshBR", "-mBR", "Mesh_B_Restriction")) {
		input_params.mesh_restrict_micro_file = field_parser.next(
				input_params.mesh_restrict_micro_file);
	} else {
		homemade_error_msg("Missing the restricted B mesh file!");
	}

	if (field_parser.search(3, "--meshI", "-mI", "MeshInter")) {
		input_params.mesh_inter_file = field_parser.next(
				input_params.mesh_inter_file);
	} else {
		homemade_error_msg("Missing the intersection mesh file!");
	}

	if ( field_parser.search(3, "--meshWeight", "-mW", "MeshWeight") )
	{
		input_params.mesh_weight_file = field_parser.next(input_params.mesh_weight_file);
	}
	else
	{
		homemade_error_msg("Missing the weight mesh file!");
	}

	// Set the equivalence and intersection tables
	if (field_parser.search(2, "--tableRA", "Mesh_A_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_BIG_file = field_parser.next(
				input_params.equivalence_table_restrict_BIG_file);
	} else {
		homemade_error_msg("Missing the equivalence table for the mesh A!");
	}

	if (field_parser.search(2, "--tableRB", "Mesh_B_RestrictionEquivalenceTable")) {
		input_params.equivalence_table_restrict_micro_file = field_parser.next(
				input_params.equivalence_table_restrict_micro_file);
	} else {
		homemade_error_msg("Missing the equivalence table for the mesh B!");
	}

	if (field_parser.search(2, "--tableFullI", "FullIntersectionElementsTable")) {
		input_params.intersection_table_full = field_parser.next(
				input_params.intersection_table_full);
	} else {
		homemade_error_msg("Missing the full intersection elements file!");
	}

	// Set table files
	if( field_parser.search(2, "--weightIdx", "WeightIndexes") )
	{
		input_params.weight_domain_idx_file = field_parser.next(input_params.weight_domain_idx_file);
	}
	else
	{
		homemade_error_msg("Missing the weight value file!");
	}

	// Set the mediator mesh
	input_params.b_UseMesh_BIG_AsMediator = false;
	input_params.b_UseMesh_micro_AsMediator = false;
	input_params.b_UseMesh_extra_AsMediator = false;
	input_params.b_Repartition_micro = true;
	if (field_parser.search(1,"Use_A_AsMediator"))
	{
		input_params.b_UseMesh_BIG_AsMediator = true;
	}
	if (field_parser.search(1,"Use_B_AsMediator"))
	{
		input_params.b_UseMesh_micro_AsMediator = true;
	}
	if (field_parser.search(1,"Use_extra_AsMediator"))
	{
		input_params.b_UseMesh_extra_AsMediator = true;
	}
	if(input_params.b_UseMesh_BIG_AsMediator
			+ input_params.b_UseMesh_micro_AsMediator
			+ input_params.b_UseMesh_extra_AsMediator > 1)
	{
		homemade_error_msg("Choose only one mesh as mediator!");
	}
	if (field_parser.search(1,"Do_notRepartition_Micro"))
	{
		input_params.b_Repartition_micro = false;
	}

	if(input_params.b_UseMesh_BIG_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_BIG_file;
		input_params.equivalence_table_mediator = input_params.equivalence_table_restrict_BIG_file;
	}
	if(input_params.b_UseMesh_micro_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_micro_file;
		input_params.equivalence_table_mediator = input_params.equivalence_table_restrict_micro_file;
	}
	if(input_params.b_UseMesh_extra_AsMediator)
	{
//		if (field_parser.search(3, "--meshM", "-mM", "MeshMediator")) {
//			input_params.mesh_mediator_file = field_parser.next(
//					input_params.mesh_mediator_file);
//		} else {
//			homemade_error_msg("Missing the mediator mesh file!");
//		}
		libmesh_not_implemented_msg("Still implementing the external mesh case!");
	}

	// Set constant parameters
	if ( field_parser.search(3, "-p","--parameters","PhysicalParameters") )
	{
		input_params.physical_params_file = field_parser.next(input_params.physical_params_file);
	}
	else
	{
		homemade_error_msg("Missing the physical parameters file!");
	}

	input_params.mean_distance = 0.2;

	if (field_parser.search(2, "--dist", "CouplingMeshScale")) {
		input_params.mean_distance = field_parser.next(
				input_params.mean_distance);
	}

	// Set coupling parameters
	bool bUseSameSearchCoeff = false;

	input_params.mean_distance = 0.2;

	input_params.k_dA = 2.5;
	input_params.k_dB = 2.5;
	input_params.k_cA = 2.5;
	input_params.k_cB = 2.5;

	if( field_parser.search(1,"LATINUseSameSearchParameters") )
	{
		bUseSameSearchCoeff = field_parser.next(bUseSameSearchCoeff);
	}

	if( bUseSameSearchCoeff )
	{
		if( field_parser.search(2, "-k","LATINk") )
		{
			input_params.k_dA = field_parser.next(input_params.k_dA);
			input_params.k_dB = input_params.k_dA;
			input_params.k_cA = input_params.k_dA;
			input_params.k_cB = input_params.k_dA;
		}
	}
	else
	{
		if( field_parser.search(2, "-kdA","LATINkDecoupledA") )
		{
			input_params.k_dA = field_parser.next(input_params.k_dA);
		}

		if( field_parser.search(2, "-kdB","LATINkDecoupledB") )
		{
			input_params.k_dB = field_parser.next(input_params.k_dB);
		}

		if( field_parser.search(2, "-kcA","LATINkCoupledA") )
		{
			input_params.k_cA = field_parser.next(input_params.k_cA);
		}

		if( field_parser.search(2, "-kcB","LATINkCoupledB") )
		{
			input_params.k_cB = field_parser.next(input_params.k_cB);
		}
	}

	if( field_parser.search(2, "--dist","CouplingMeshScale") )
	{
		input_params.mean_distance = field_parser.next(input_params.mean_distance);
	}

	// Set LATIN parameters
	input_params.LATIN_b_UseRestartFiles = false;
	input_params.LATIN_b_PrintRestartFiles = false;
	input_params.LATIN_eps = 1E-2;
	input_params.LATIN_conv_max = 10000;
	input_params.LATIN_relax = 0.8;

	input_params.coupled_convergence_output = "coupled_convergence.dat";
	input_params.coupled_restart_file_base = "coupled_restart";

	if( field_parser.search(2, "--LATINeps","LATINEps") )
	{
		input_params.LATIN_eps = field_parser.next(input_params.LATIN_eps);
	}

	if( field_parser.search(2, "--LATINconv","LATINConvergneceLimit") )
	{
		input_params.LATIN_conv_max = field_parser.next(input_params.LATIN_conv_max);
	}

	if( field_parser.search(2, "--LATINrelax","LATINRelax") )
	{
		input_params.LATIN_relax = field_parser.next(input_params.LATIN_relax);
	}

	if( field_parser.search(2, "--CoupledConvOutput","CoupledConvergenceOutput") )
	{
		input_params.coupled_convergence_output = field_parser.next(input_params.coupled_convergence_output);
	}

	// Set coupling solver convergence
	input_params.coupled_conv_abs = 1e-20;
	input_params.coupled_conv_rel = 1e-5;
	input_params.coupled_div = 1e5;
	input_params.coupled_iter_max = 1e4;

	if( field_parser.search(1,"CoupledConvAbs") )
	{
		input_params.coupled_conv_abs = field_parser.next(input_params.coupled_conv_abs);
	}
	if( field_parser.search(1,"CoupledConvRel") )
	{
		input_params.coupled_conv_rel = field_parser.next(input_params.coupled_conv_rel);
	}
	if( field_parser.search(1,"CoupledDiv") )
	{
		input_params.coupled_div = field_parser.next(input_params.coupled_div);
	}
	if( field_parser.search(1,"CoupledIterMax") )
	{
		input_params.coupled_iter_max = field_parser.next(input_params.coupled_iter_max);
	}

	if( field_parser.search(1,"Use_restart_data") )
	{
		input_params.LATIN_b_UseRestartFiles = true;
	}
	if( field_parser.search(1,"Save_restart_data") )
	{
		input_params.LATIN_b_PrintRestartFiles = true;
	}

	if( field_parser.search(1,"CoupledRestartDataFiles") )
	{
		input_params.coupled_restart_file_base = field_parser.next(input_params.coupled_restart_file_base);
	}

	if( field_parser.search(1,"SkipOutput") )
	{
		input_params.b_PrintOutput = false;
	}
	else
	{
		input_params.b_PrintOutput = true;
	}

	// Set output files
	input_params.output_file_BIG = "meshes/3D/output/carl_multi_crystal_test_micro.exo";
	input_params.output_file_micro = "meshes/3D/output/carl_multi_crystal_test_macro.exo";
	if ( field_parser.search(3, "-oA","--outputA", "OutputEXOFileA") )
	{
		input_params.output_file_BIG = field_parser.next(input_params.output_file_BIG);
	}
	if ( field_parser.search(3, "-oB","--outputB", "OutputEXOFileB") )
	{
		input_params.output_file_micro = field_parser.next(input_params.output_file_micro);
	}

	input_params.solver_type = carl::CoupledSolverType::LATIN_MODIFIED_STIFFNESS;
	input_params.CG_precond_type = carl::BaseCGPrecondType::NO_PRECONDITIONER;
	if ( field_parser.search(1, "SolverType") )
	{
		input_params.solver_type_string = field_parser.next(input_params.solver_type_string);
		if(input_params.solver_type_string == "LATIN_Modified")
			input_params.solver_type = carl::CoupledSolverType::LATIN_MODIFIED_STIFFNESS;
		if(input_params.solver_type_string == "LATIN_Original_Stiffness")
			input_params.solver_type = carl::CoupledSolverType::LATIN_ORIGINAL_STIFFNESS;
		if(input_params.solver_type_string == "CG")
		{
			input_params.solver_type = carl::CoupledSolverType::CG;
			if ( field_parser.search(1, "CGPreconditionerType") )
			{
				input_params.CG_precond_type_string = field_parser.next(input_params.CG_precond_type_string);
				if(input_params.CG_precond_type_string == "NONE")
					input_params.CG_precond_type = carl::BaseCGPrecondType::NO_PRECONDITIONER;
				if(input_params.CG_precond_type_string == "Coupling_operator")
					input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_OPERATOR;
				if(input_params.CG_precond_type_string == "Coupled_system_operator")
					input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLED_SYSTEM_OPERATOR;
				if(input_params.CG_precond_type_string == "Coupling_operator_jacobi")
					input_params.CG_precond_type = carl::BaseCGPrecondType::JACOBI;

			}
		}
	}

	if( field_parser.search(1,"ExportScalingData") )
	{
		input_params.b_ExportScalingData = true;
		input_params.scaling_data_file = field_parser.next(input_params.scaling_data_file);
	}
	else
	{
		input_params.b_ExportScalingData = false;
	}
};

};
#endif /* CARL_INPUT_PARSER_H_ */
