/*
 * carl_feti_iterate_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef CARL_FETI_ITERATE_INPUT_PARSER_H_
#define CARL_FETI_ITERATE_INPUT_PARSER_H_

#include "carl_headers.h"

namespace carl
{
/// Structure containing the parameters for the setup initialization of the FETI solver.
struct feti_iterate_params {
	// --- Parameters used directly by the CArl_FETI_iterate program (some are also used by the other CArl_FETI programs)

	// Cluster 
	ClusterSchedulerType scheduler; ///< Cluster scheduler software type. *Values*: PBS, SLURM (code not implemented for the later yet).

	// External solver commands
	std::string ext_solver_BIG;			///< Command used for the external solver for system A.
	std::string ext_solver_micro;		///< Command used for the external solver for system B.
	
	// Path to "scratch" folder
	std::string scratch_folder_path;	///< Path to the folder which will be used to save the temporary files during the solve operation

	// Rigid body mode options for the micro system
	bool bUseRigidBodyModes;			///< [RB] Use the rigid body modes for the micro system?
	std::string RB_vectors_base;		///< [RB] Common path base for the micro system's rigid body mode vectors.

	// Coupling matrices path
	std::string coupling_path_base;		///< Base of the coupling matrices path.

	// --- Parameters used to set up the other CArl_FETI input files

	// FETI / CG parameters
	double CG_coupled_conv_abs;		///< [CG] Absolute residual convergence.
	double CG_coupled_conv_rel;		///< [CG] Relative residual convergence.
	double CG_coupled_div;			///< [CG] Residual divergence.
	double CG_coupled_conv_corr;	///< [CG] Relative rigid body mode convergence.
	int CG_coupled_conv_max;		///< [CG] Maximum number of iterations.	
	
	carl::BaseCGPrecondType CG_precond_type;	///< [CG] Type of preconditionner.
};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *    - `ClusterSchedulerType` : scheduler type. *Values*: PBS or SLURM (code not implemented for the later yet).
 *	  - `ExtSolverA` : command line for the external solver for system A.
 *	  - `ExtSolverB` : command line for the external solver for system B.
 *	  - `ScratchFolderPath` : path to the folder where the temporary files used by the coupled solver will be saved.
 *    - `CouplingMatricesBase` : filename base of the coupling matrices files.
 *
 *  Boolean flags:
 *    - `UseRigidBodyModesB` : use the rigid body modes for system B.
 *
 *  Rigid body mode parameters (only read if `UseRigidBodyModesB` is used):
 *	  - `ExtForceSystemB` : path to the vector containing the external forces for the system B.
 *	  - `RBVectorBase` : filename base of the rigid body modes vectors.
 *
 *  Optional parameters:
 *  + FETI / CG optional parameters:
 *    - `CGPreconditionerType` : CG preconditioner type. *Values*: "NONE", "Coupling_operator" or "Coupling_operator_jacobi". *Default*: "Coupling_operator".
 *    - `CoupledConvAbs` : CG absolute convergence on the residual. *Default*: 1e-20.
 *    - `CoupledConvRel` : CG relative convergence on the residual.  *Default*: 1e-5.
 *    - `CoupledCorrConvRel` : CG relative convergence on the rigid body corrections. *Default*: 1e-6.
 *    - `CoupledDiv` : CG residual divergence parameter. *Default*: 100000.
 *    - `CoupledIterMax` : CG maximum number of iterations. *Default*: 1000.
 */
void get_input_params(GetPot& field_parser,
		feti_iterate_params& input_params) {

	if (field_parser.search(1, "ClusterSchedulerType")) {
		std::string cluster_scheduler_type;
		cluster_scheduler_type = field_parser.next(cluster_scheduler_type);
		if(cluster_scheduler_type == "PBS")
			input_params.scheduler = carl::ClusterSchedulerType::PBS;
		else if(cluster_scheduler_type == "SLURM")
			input_params.scheduler = carl::ClusterSchedulerType::SLURM;
		else
			homemade_error_msg("Invalid scheduler type!");
	} else {
		homemade_error_msg("Missing the scheduler type!");
	}

	if (field_parser.search(1, "ExtSolverA")) {
		input_params.ext_solver_BIG = field_parser.next(
				input_params.ext_solver_BIG);
	} else {
		homemade_error_msg("Missing the external solver A command line!");
	}

	if (field_parser.search(1, "ExtSolverB")) {
		input_params.ext_solver_micro = field_parser.next(
				input_params.ext_solver_micro);
	} else {
		homemade_error_msg("Missing the external solver B command line!");
	}

	if (field_parser.search(1, "ScratchFolderPath")) {
		input_params.scratch_folder_path = field_parser.next(
				input_params.scratch_folder_path);
	} else {
		homemade_error_msg("Missing the external scratch folder path!");
	}

	if (field_parser.search(1, "CouplingMatricesBase")) {
		input_params.coupling_path_base = field_parser.next(
				input_params.coupling_path_base);
	} else {
		homemade_error_msg("Missing the coupling matrices path!");
	}

	if (field_parser.search(1,"UseRigidBodyModesB"))
	{
		input_params.bUseRigidBodyModes = true;
		if (field_parser.search(1, "RBVectorBase")) {
			input_params.RB_vectors_base = field_parser.next(
					input_params.RB_vectors_base);
		} else {
			homemade_error_msg("Missing the system B's rigid body mode vectors!");
		}
	}
	else
	{
		input_params.bUseRigidBodyModes = false;
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
};

};
#endif /* CARL_FETI_ITERATE_INPUT_PARSER_H_ */
