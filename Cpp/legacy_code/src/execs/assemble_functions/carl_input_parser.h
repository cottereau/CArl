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
/**	\brief Structure containing the parameters for the coupled solver test programs.
 *	
 *		Details on the parameters setup are found in the documentation of carl::get_input_params(GetPot& field_parser,
 *		coupling_generation_input_params& input_params).
 */
struct coupling_generation_input_params {
	std::string physical_params_file;	///< Physical parameters.

	std::string mesh_BIG_file;		///< Path to the macro (BIG) system mesh.
	std::string mesh_micro_file;	///< Path to the micro system mesh.

	std::string mesh_restrict_BIG_file;		///< Path to the restricted macro (BIG) system mesh.
	std::string mesh_restrict_micro_file;	///< Path to the restricted micro system mesh.

	std::string mesh_mediator_file;		///< Path to the mediator system mesh.
	std::string mesh_inter_file;		///< Common path to the intersection meshes.

	std::string mesh_weight_file;		///< Path to the mesh containing the weight region indices.

	/// Equivalence table between the macro (BIG) system mesh and its restriction mesh
	std::string equivalence_table_restrict_BIG_file;
	/// Equivalence table between the micro system mesh and its restriction mesh
	std::string equivalence_table_restrict_micro_file;
	/// Common path to the intersection tables
	std::string intersection_table_full;

	/// Path to the weight domain index files
	std::string weight_domain_idx_file;

	/// Path to the scaling data file (if the scaling is done)
	std::string scaling_data_file;

	bool b_UseMesh_BIG_AsMediator;		///< Use the macro (BIG) system's restricted mesh as the mediator? *Default*: true.
	bool b_UseMesh_micro_AsMediator;	///< Use the micro system's restricted mesh as the mediator? *Default*: false.
	bool b_UseMesh_extra_AsMediator;	///< Use an external mesh as the mediator? *Default*: false.
	bool b_Repartition_micro;			///< Repartition the micro system? *Default*: false. (DO NOT USE IT WITH CG!!!)
	bool b_UseRestartFiles;				///< Use the restart files? *Default*: false.
	bool b_PrintRestartFiles;			///< Export the restart files? *Default*: true.
	bool b_PrintOutput;					///< Save output? *Default*: true.
	bool b_ExportScalingData;			///< Export scaling data? *Default*: false.

	double mean_distance;	///< Mean element length (term `e` in the L2 coupling).

	// LATIN parameters 
	double k_dA;	///< [LATIN] Search parameter (decoupled, macro system).
	double k_dB;	///< [LATIN] Search parameter (decoupled, micro system).
	double k_cA;	///< [LATIN] Search parameter (coupled, macro system).
	double k_cB;	///< [LATIN] Search parameter (coupled, micro system).

	double LATIN_eps;		///< [LATIN] Precision.
	int LATIN_conv_max;		///< [LATIN] Maximum number of iterations.
	double LATIN_relax;		///< [LATIN] Relaxation.

	// More general / CG parameters
	double CG_coupled_conv_abs;		///< [CG] Absolute residual convergence.
	double CG_coupled_conv_rel;		///< [CG] Relative residual convergence.
	double CG_coupled_div;			///< [CG] Residual divergence.
	double CG_coupled_conv_corr;	///< [CG] Relative rigid body mode convergence.
	int CG_coupled_conv_max;		///< [CG] Maximum number of iterations.	
	
	carl::BaseCGPrecondType CG_precond_type;	///< [CG] Type of preconditionner.

	std::string coupled_convergence_output;		///< Convergence data output.
	std::string coupled_restart_file_base;		///< Common path to the solver restart files.

	std::string output_file_BIG;		///< Output file for the macro (BIG) system mesh.
	std::string output_file_micro;		///< Output file for the micro system mesh.

	carl::CoupledSolverType solver_type;	///< Coupled solver type.
};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *  + System and intersection meshes:
 *	  - `MeshA`, `-mA` or `--meshA` : path to the mesh A.
 *	  - `MeshB`, `-mB` or `--meshB` : path to the mesh B.
 *	  - `MeshI`, `-mI` or `--meshI` : common path to the intersection meshes and tables.
 *  + Arlequin weight parameters:
 *	  - `MeshWeight`, `-mW` or `--meshWeight` : path to the mesh defining the domains of the Arlequin weight parameters.
 *    - `WeightIndexes` or `--weightIdx` : path to the indices of the domains of the Arlequin weight parameters.
 *  + Restricted meshes and tables:
 *	  - `Mesh_A_Restriction`, `-mAR` or `--meshAR` : path to the restricted mesh A (formed by elements of the mesh A intersecting the coupling region).
 *	  - `Mesh_B_Restriction`, `-mBR` or `--meshBR` : path to the restricted mesh B (formed by elements of the mesh A intersecting the coupling region).
 *    - `Mesh_A_RestrictionEquivalenceTable` or `--tableRA` : path to the equivalence table between the mesh A and its restriction.
 *    - `Mesh_B_RestrictionEquivalenceTable` or `--tableRB` : path to the equivalence table between the mesh B and its restriction.
 *  + Physical parameters:
 *    - `PhysicalParameters`, `-p` or `--parameters` : physical parameters.
 *
 *  Optional parameters:
 *  + Mediator mesh boolean flags (only one can be used at each time)
 *    - `Use_A_AsMediator` : use the restricted mesh A for the mediator space (set by default).
 *    - `Use_B_AsMediator` : use the restricted mesh B for the mediator space.
 *    - `Use_extra_AsMediator` : use an extra mesh for the mediator (NOT IMPLEMENTED YET!).
 *  + Coupling parameters:
 *    - `SolverType` : coupled solver type. *Values*: "CG", "LATIN_MODIFIED_STIFFNESS" and "LATIN_ORIGINAL_STIFFNESS" (does not modifiy the stiffness matrices, **deprecated**). *Default*: "CG".
 *    - `CouplingMeshScale` or `--dist` : Mean element length (term `e` in the L2 coupling). *Default*: 0.2.
 *    - `CoupledConvergenceOutput` or `--CoupledConvOutput` : path for the convergence information file. *Default*: "coupled_convergence.dat".
 *  + Output and restart parameters:
 *    - `SkipOutput` (boolean flag): add this flag to skip writing the outputs.
 *    - `OutputEXOFileA`, `-oA` or `--outputA` : path to the output mesh for the system A. *Default*: "test_macro.exo".
 *    - `OutputEXOFileB`, `-oB` or `--outputB` : path to the output mesh for the system B. *Default*: "test_micro.exo".
 *    - `Use_restart_data` (boolean flag): add this flag to start the algorithm from the restart files.
 *    - `Save_restart_data` (boolean flag): add this flag to save the restart files.
 *    - `CoupledRestartDataFiles` : *Default*: "coupled_restart".
 *    - `ExportScalingData`: common path to the scaling files.
 *
 *  + LATIN optional parameters:
 *    - `LATINUseSameSearchParameters` : use the same values for all the search parameters. *Default*: true.
 *    - `LATINk` or `-k` : common value for the search parameters. *Default*: 2.5.
 *    - `LATINkDecoupledA` or `-kdA` : search parameters for the system A (decoupled step). *Default*: 2.5.
 *    - `LATINkDecoupledB` or `-kdB` : search parameters for the system B (decoupled step). *Default*: 2.5.
 *    - `LATINkCoupledA` or `-kcA` : search parameters for the system A (coupled step). *Default*: 2.5.
 *    - `LATINkCoupledA` or `-kcB` : search parameters for the system B (coupled step).  *Default*: 2.5.
 *    - `LATINEps` or `--LATINeps` : LATIN precision parameter. *Default*: 1e-2.
 *    - `LATINConvergneceLimit` or `--LATINconv` : LATIN maximum number of iterations. *Default*: 1e4.
 *    - `LATINRelax` or `--LATINrelax` : LATIN relaxation. *Default*: 0.8.
 *  + CG optional parameters:
 *    - `CGPreconditionerType` : CG preconditioner type. *Values*: "NONE", "Coupling_operator" or "Coupling_operator_jacobi". *Default*: "NONE".
 *    - `CoupledConvAbs` : CG absolute convergence on the residual. *Default*: 1e-20.
 *    - `CoupledConvRel` : CG relative convergence on the residual.  *Default*: 1e-5.
 *    - `CoupledCorrConvRel` : CG relative convergence on the rigid body corrections. *Default*: 1e-6.
 *    - `CoupledDiv` : CG residual divergence parameter. *Default*: 1e5.
 *    - `CoupledIterMax` : CG maximum number of iterations. *Default*: 1e4.
 *  + **Useless / possible candidate for deletion / for all that's sacred do not use**.
 *    - `Repartition_Micro` : repartition the micro system mesh, following the modifications due to the coupling (DO NOT USE WITH CG!!!) 
 */
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
		input_params.intersection_table_full = input_params.mesh_inter_file;
	} else {
		homemade_error_msg("Missing the intersection files path!");
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
	input_params.b_UseMesh_BIG_AsMediator = true;
	input_params.b_UseMesh_micro_AsMediator = false;
	input_params.b_UseMesh_extra_AsMediator = false;

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

	input_params.b_Repartition_micro = false;
	if (field_parser.search(1,"Repartition_Micro"))
	{
		input_params.b_Repartition_micro = true;
	}

	if(input_params.b_UseMesh_BIG_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_BIG_file;
	}
	if(input_params.b_UseMesh_micro_AsMediator)
	{
		input_params.mesh_mediator_file = input_params.mesh_restrict_micro_file;
	}
	if(input_params.b_UseMesh_extra_AsMediator)
	{
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

	// Set coupling parameters
	input_params.mean_distance = 0.2;

	if( field_parser.search(2, "--dist","CouplingMeshScale") )
	{
		input_params.mean_distance = field_parser.next(input_params.mean_distance);
	}

	input_params.coupled_convergence_output = "coupled_convergence.dat";
	if( field_parser.search(2, "--CoupledConvOutput","CoupledConvergenceOutput") )
	{
		input_params.coupled_convergence_output = field_parser.next(input_params.coupled_convergence_output);
	}

	// Set LATIN parameters
		bool bUseSameSearchCoeff = false;

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

	input_params.LATIN_eps = 1E-2;
	input_params.LATIN_conv_max = 10000;
	input_params.LATIN_relax = 0.8;

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

	// Set CG coupling solver convergence
	input_params.CG_coupled_conv_abs = 1e-20;
	input_params.CG_coupled_conv_rel = 1e-5;
	input_params.CG_coupled_div = 1e5;
	input_params.CG_coupled_conv_max = 1e4;
	input_params.CG_coupled_conv_corr =1e-6;

	if( field_parser.search(1,"CoupledConvAbs") )
	{
		input_params.CG_coupled_conv_abs = field_parser.next(input_params.CG_coupled_conv_abs);
	}
	if( field_parser.search(1,"CoupledConvRel") )
	{
		input_params.CG_coupled_conv_rel = field_parser.next(input_params.CG_coupled_conv_rel);
	}
	if( field_parser.search(1,"CoupledCorrConvRel") )
	{
		input_params.CG_coupled_conv_corr = field_parser.next(input_params.CG_coupled_conv_corr);
	}
	if( field_parser.search(1,"CoupledDiv") )
	{
		input_params.CG_coupled_div = field_parser.next(input_params.CG_coupled_div);
	}
	if( field_parser.search(1,"CoupledIterMax") )
	{
		input_params.CG_coupled_conv_max = field_parser.next(input_params.CG_coupled_conv_max);
	}

	// Restart files
	input_params.b_UseRestartFiles = false;
	input_params.b_PrintRestartFiles = false;
	input_params.coupled_restart_file_base = "coupled_restart";

	if( field_parser.search(1,"Use_restart_data") )
	{
		input_params.b_UseRestartFiles = true;
	}
	if( field_parser.search(1,"Save_restart_data") )
	{
		input_params.b_PrintRestartFiles = true;
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
	input_params.output_file_BIG = "test_macro.exo";
	input_params.output_file_micro = "test_micro.exo";
	if ( field_parser.search(3, "-oA","--outputA", "OutputEXOFileA") )
	{
		input_params.output_file_BIG = field_parser.next(input_params.output_file_BIG);
	}
	if ( field_parser.search(3, "-oB","--outputB", "OutputEXOFileB") )
	{
		input_params.output_file_micro = field_parser.next(input_params.output_file_micro);
	}

	input_params.solver_type = carl::CoupledSolverType::CG;
	input_params.CG_precond_type = carl::BaseCGPrecondType::NO_PRECONDITIONER;
	if ( field_parser.search(1, "SolverType") )
	{
		std::string solver_type_string = field_parser.next(solver_type_string);
		if(solver_type_string == "LATIN_Modified")
			input_params.solver_type = carl::CoupledSolverType::LATIN_MODIFIED_STIFFNESS;
		if(solver_type_string == "LATIN_Original_Stiffness")
			input_params.solver_type = carl::CoupledSolverType::LATIN_ORIGINAL_STIFFNESS;
		if(solver_type_string == "CG")
		{
			input_params.solver_type = carl::CoupledSolverType::CG;
			if ( field_parser.search(1, "CGPreconditionerType") )
			{
				std::string CG_precond_type_string = field_parser.next(CG_precond_type_string);
				if(CG_precond_type_string == "NONE")
					input_params.CG_precond_type = carl::BaseCGPrecondType::NO_PRECONDITIONER;
				if(CG_precond_type_string == "Coupling_operator")
					input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_OPERATOR;
				if(CG_precond_type_string == "Coupling_operator_jacobi")
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
