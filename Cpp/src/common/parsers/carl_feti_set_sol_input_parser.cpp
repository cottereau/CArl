/*
 * carl_feti_set_sol_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */


#include "carl_feti_set_sol_input_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
		feti_set_sol_params& input_params) {

	if (field_parser.search(1, "ClusterSchedulerType")) {
		std::string cluster_scheduler_type;
		cluster_scheduler_type = field_parser.next(cluster_scheduler_type);
		if(cluster_scheduler_type == "LOCAL")
		{
			std::cout << " !!! WARNING: " << std::endl;
			std::cout << "        Using the LOCAL job 'scheduler'. You will have to launch each script" << std::endl;
			std::cout << "     MANUALLY!!! Reason: MPI does not support recursive 'mpirun' calls" << std::endl;
			input_params.scheduler = carl::ClusterSchedulerType::LOCAL;
		}
		else if(cluster_scheduler_type == "PBS")
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

	if (field_parser.search(1, "ScratchFolderPathTimeMono")) {
		input_params.scratch_folder_path_time_mono = field_parser.next(
				input_params.scratch_folder_path_time_mono);
	} else {
		homemade_error_msg("Missing the external scratch folder path for Time_Mono!");
	}

	if (field_parser.search(1,"UseRigidBodyModesB"))
	{
		input_params.bUseRigidBodyModes = true;
	} else {
		input_params.bUseRigidBodyModes = false;
	}

	if (field_parser.search(1,"OutputFolder")) {
		input_params.output_folder = field_parser.next(
				input_params.output_folder);
	} else {
		homemade_error_msg("Missing the output filename base!");
	}
};

};
