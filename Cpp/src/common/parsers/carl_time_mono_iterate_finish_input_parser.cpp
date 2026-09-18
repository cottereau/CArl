/*
 * carl_time_mono_iterate_finish_input_parser.cpp
 *
 *  Created on: August 4, 2026
 *      Author: Romain Ruyssen
 */

#include "carl_time_mono_iterate_finish_input_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
		time_mono_iterate_finish_params& input_params) {

	//-------------------------------------------------------
	// Cluster scheduler type
	if (field_parser.search(1, "SchedulerType")) {
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

	//-------------------------------------------------------		
	// Path to "scratch" folder	
	if (field_parser.search(1, "ScratchFolderPath")) {
		input_params.scratch_folder_path = field_parser.next(
				input_params.scratch_folder_path);
		std::cout << input_params.scratch_folder_path << std::endl;
	} else {
		homemade_error_msg("Missing the scratch folder path!");
	}		

	//-------------------------------------------------------
	// Path for the copying function of the FETI outputs
	if (field_parser.search(1, "SolutionsFolderPath")) {
		input_params.feti_solution_path = field_parser.next(
				input_params.feti_solution_path);
		std::cout << input_params.feti_solution_path << std::endl;	
	}
	else {
		homemade_error_msg("Missing the FETI solution path!");
	}

	if (field_parser.search(1, "ResultsFolderPath")) {
		input_params.results_folder_path = field_parser.next(
				input_params.results_folder_path);
		std::cout << input_params.results_folder_path << std::endl;	
	}
	else {
		homemade_error_msg("Missing the results folder path!");
	}

	if (field_parser.search(1, "ResultsFileName")) {
		input_params.results_file_name = field_parser.next(
				input_params.results_file_name);
		std::cout << input_params.results_file_name << std::endl;	
	}
	else {
		homemade_error_msg("Missing the results file name!");
	}

	//-------------------------------------------------------
	// Path to the time mesh file for the time loop completion check	
	if (field_parser.search(1, "TimeMeshFile")) {
		input_params.time_mesh_file = field_parser.next(
				input_params.time_mesh_file);
		std::cout << input_params.time_mesh_file << std::endl;	
	}
	else {
		homemade_error_msg("Missing the time mesh file!");
	}

	};
};