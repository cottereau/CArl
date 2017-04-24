/*
 * carl_feti_setup_finish_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "carl_feti_setup_finish_input_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
		feti_setup_finish_params& input_params) {

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
	} else {
		input_params.bUseRigidBodyModes = false;
	}

	if ( field_parser.search(1, "CGPreconditionerType") )
	{
		std::string CG_precond_type_string = field_parser.next(CG_precond_type_string);
		if(CG_precond_type_string == "NONE")
			input_params.CG_precond_type = carl::BaseCGPrecondType::NO_PRECONDITIONER;
		else if(CG_precond_type_string == "Coupling_operator")
			input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_OPERATOR;
		else if(CG_precond_type_string == "Coupling_operator_jacobi")
			input_params.CG_precond_type = carl::BaseCGPrecondType::JACOBI;
		else
			homemade_error_msg("Invalid preconditionner type!");
	} else {
		homemade_error_msg("Missing preconditionner type!");
	}
};

};
