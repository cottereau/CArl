/*
 * carl_feti_iterate_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "carl_feti_iterate_input_parser.h"

namespace carl
{

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
				
		if (field_parser.search(1, "NbOfRBVectors")) {
			input_params.nb_of_rb_vectors = field_parser.next(
					input_params.nb_of_rb_vectors);
		} else {
			input_params.nb_of_rb_vectors = 6;
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

	if ( field_parser.search(1, "CGPreconditionerType") )
	{
		std::string CG_precond_type_string = field_parser.next(CG_precond_type_string);
		if(CG_precond_type_string == "NONE")
			input_params.CG_precond_type = carl::BaseCGPrecondType::NO_PRECONDITIONER;
		else if(CG_precond_type_string == "Coupling_operator")
			input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_OPERATOR;
		else if(CG_precond_type_string == "Coupling_operator_jacobi")
			input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_JACOBI;
		else
			homemade_error_msg("Invalid preconditionner type!");
	} else {
		homemade_error_msg("Missing preconditionner type!");
	}
};

};

