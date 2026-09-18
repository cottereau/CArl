/*
 * carl_time_mono_setup_input_parser.cpp
 *
 *  Created on: July 22, 2026
 *      Author: Romain Ruyssen
 */

#include "carl_time_mono_setup_input_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
		time_mono_setup_params& input_params) {

	//-------------------------------------------------------
	// Cluster scheduler type
	if (field_parser.search(1, "ClusterSchedulerType")) {
		std::string cluster_scheduler_type;
		cluster_scheduler_type = field_parser.next(cluster_scheduler_type);
		if(cluster_scheduler_type == "LOCAL")
		{
			std::cout << " !!! WARNING: " << std::endl;
			std::cout << "        The LOCAL 'scheduler' type is only intended for small and fast test cases" << std::endl;
			std::cout << "     on computers without a job scheduler (PBS, SLURM). You will have to launch" << std::endl;
			std::cout << "     each script MANUALLY!!! Reason: MPI does not support recursive 'mpirun' calls" << std::endl;
			input_params.scheduler = carl::ClusterSchedulerType::LOCAL;
			input_params.script_filename = "";
		}
		else if(cluster_scheduler_type == "PBS") {
			input_params.scheduler = carl::ClusterSchedulerType::PBS;
			if (field_parser.search(1, "ScriptFile")) {
				input_params.script_filename = field_parser.next(input_params.script_filename);
			} else {
				homemade_error_msg("Missing the script file (needed for the PBS scheduler)!");
			}
		}
		else if(cluster_scheduler_type == "SLURM") {
			input_params.scheduler = carl::ClusterSchedulerType::SLURM;
			if (field_parser.search(1, "ScriptFile")) {
				input_params.script_filename = field_parser.next(input_params.script_filename);
			} else {
				homemade_error_msg("Missing the script file (needed for the SLURM scheduler)!");
			}
		}
		else
			homemade_error_msg("Invalid scheduler type!");
	} else {
		homemade_error_msg("Missing the scheduler type!");
	}

	//-------------------------------------------------------		
	// Dossier temporaire		
	if (field_parser.search(1, "ScratchFolderPath")) {
		input_params.scratch_folder_path = field_parser.next(
				input_params.scratch_folder_path);
		std::cout << input_params.scratch_folder_path << std::endl;
	} else {
		homemade_error_msg("Missing the scratch folder path!");
	}		

	//-------------------------------------------------------
	// Commandes pour l'assembleur du second membre des solveurs externes
	if (field_parser.search(1, "ExtSolverRhsAssemblyA")) {
		input_params.ext_solver_BIG_rhs_assembly = field_parser.next(
				input_params.ext_solver_BIG_rhs_assembly);
		std::cout << input_params.ext_solver_BIG_rhs_assembly << std::endl;
	} else {
		homemade_error_msg("Missing the BIG external solver command for rhs assembly!");
	}		
	if (field_parser.search(1, "ExtSolverRhsAssemblyB")) {
		input_params.ext_solver_micro_rhs_assembly = field_parser.next(
				input_params.ext_solver_micro_rhs_assembly);
		std::cout << input_params.ext_solver_micro_rhs_assembly << std::endl;
	} else {
		homemade_error_msg("Missing the micro external solver command for rhs assembly!");
	}	

	//-------------------------------------------------------
	// Commandes pour la mise à jour des conditions initiales des solveurs externes
	if (field_parser.search(1, "ExtSolverInitCondUpdateA")) {
		input_params.ext_solver_BIG_update_init_cond = field_parser.next(
				input_params.ext_solver_BIG_update_init_cond);
		std::cout << input_params.ext_solver_BIG_update_init_cond << std::endl;
	} else {
		homemade_error_msg("Missing the BIG external solver command for update of initial conditions!");
	}		
	if (field_parser.search(1, "ExtSolverInitCondUpdateB")) {
		input_params.ext_solver_micro_update_init_cond = field_parser.next(
				input_params.ext_solver_micro_update_init_cond);
		std::cout << input_params.ext_solver_micro_update_init_cond << std::endl;
	} else {
		homemade_error_msg("Missing the micro external solver command for update of initial conditions!");
	}

	//-------------------------------------------------------
	// Fichiers contenant le maillage temporel
	if (field_parser.search(1, "TimeMeshFile")) {
		input_params.time_mesh_file = field_parser.next(
				input_params.time_mesh_file);
		std::cout << input_params.time_mesh_file << std::endl;
	} else {
		homemade_error_msg("Missing the time mesh file !");
	}

	//-------------------------------------------------------
	// Fichiers contenant les paramètres de Newmark
	if (field_parser.search(1, "NewmarkParameters_Macro")) {
		input_params.newmark_parameters_A_file = field_parser.next(
				input_params.newmark_parameters_A_file);
		std::cout << input_params.newmark_parameters_A_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro Newmark parameters file !");
	}

	if (field_parser.search(1, "NewmarkParameters_Micro")) {
		input_params.newmark_parameters_B_file = field_parser.next(
				input_params.newmark_parameters_B_file);
		std::cout << input_params.newmark_parameters_B_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro Newmark parameters file !");
	}

	//-------------------------------------------------------
	// Fichiers contenant les maillages spatiaux
	if (field_parser.search(1, "Mesh_Macro")) {
		input_params.space_mesh_A_file = field_parser.next(
				input_params.space_mesh_A_file);
		std::cout << input_params.space_mesh_A_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro space mesh file !");
	}

	if (field_parser.search(1, "Mesh_Micro")) {
		input_params.space_mesh_B_file = field_parser.next(
				input_params.space_mesh_B_file);
		std::cout << input_params.space_mesh_B_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro space mesh file !");
	}

	//-------------------------------------------------------
	// Fichiers contenant les différentes matrices des domaines
	if (field_parser.search(1, "SysMatrix_Macro")) {
		input_params.domain_A_matrix_file = field_parser.next(
				input_params.domain_A_matrix_file);
		std::cout << input_params.domain_A_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_Macro_M")) {
		input_params.domain_A_M_matrix_file = field_parser.next(
				input_params.domain_A_M_matrix_file);
		std::cout << input_params.domain_A_M_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro M matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_Macro_K")) {
		input_params.domain_A_K_matrix_file = field_parser.next(
				input_params.domain_A_K_matrix_file);
		std::cout << input_params.domain_A_K_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro K matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_Macro_C")) {
		input_params.domain_A_C_matrix_file = field_parser.next(
				input_params.domain_A_C_matrix_file);
		std::cout << input_params.domain_A_C_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro C matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_Micro")) {
		input_params.domain_B_matrix_file = field_parser.next(
				input_params.domain_B_matrix_file);
		std::cout << input_params.domain_B_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_Micro_M")) {
		input_params.domain_B_M_matrix_file = field_parser.next(
				input_params.domain_B_M_matrix_file);
		std::cout << input_params.domain_B_M_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro M matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_Micro_K")) {
		input_params.domain_B_K_matrix_file = field_parser.next(
				input_params.domain_B_K_matrix_file);
		std::cout << input_params.domain_B_K_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro K matrix file !");
	}

	if (field_parser.search(1, "SysMatrix_Micro_C")) {
		input_params.domain_B_C_matrix_file = field_parser.next(
				input_params.domain_B_C_matrix_file);
		std::cout << input_params.domain_B_C_matrix_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro C matrix file !");
	}

	//-------------------------------------------------------
	// Fichiers contenant les conditions initiales des domaines
	if (field_parser.search(1, "InitCond_Macro_disp")) {
		input_params.domain_A_init_cond_disp_file = field_parser.next(
				input_params.domain_A_init_cond_disp_file);
		std::cout << input_params.domain_A_init_cond_disp_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro displacement initial condition file !");
	}

	if (field_parser.search(1, "InitCond_Macro_vel")) {
		input_params.domain_A_init_cond_vel_file = field_parser.next(
				input_params.domain_A_init_cond_vel_file);
		std::cout << input_params.domain_A_init_cond_vel_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro velocity initial condition file !");
	}

	if (field_parser.search(1, "InitCond_Macro_acc")) {
		input_params.domain_A_init_cond_acc_file = field_parser.next(
				input_params.domain_A_init_cond_acc_file);
		std::cout << input_params.domain_A_init_cond_acc_file << std::endl;
	} else {
		homemade_error_msg("Missing the Macro acceleration initial condition file !");
	}

	if (field_parser.search(1, "InitCond_Micro_disp")) {
		input_params.domain_B_init_cond_disp_file = field_parser.next(
				input_params.domain_B_init_cond_disp_file);
		std::cout << input_params.domain_B_init_cond_disp_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro displacement initial condition file !");
	}

	if (field_parser.search(1, "InitCond_Micro_vel")) {
		input_params.domain_B_init_cond_vel_file = field_parser.next(
				input_params.domain_B_init_cond_vel_file);
		std::cout << input_params.domain_B_init_cond_vel_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro velocity initial condition file !");
	}

	if (field_parser.search(1, "InitCond_Micro_acc")) {
		input_params.domain_B_init_cond_acc_file = field_parser.next(
				input_params.domain_B_init_cond_acc_file);
		std::cout << input_params.domain_B_init_cond_acc_file << std::endl;
	} else {
		homemade_error_msg("Missing the Micro acceleration initial condition file !");
	}

	//-------------------------------------------------------
	// Fichier contenant les paramètres de setup du solveur FETI 
	if (field_parser.search(1, "FETISetupParametersFile")) {
		input_params.feti_setup_params_file = field_parser.next(
				input_params.feti_setup_params_file);
		std::cout << input_params.feti_setup_params_file << std::endl;	
	}
	else {
		homemade_error_msg("Missing the FETI setup parameters file!");
	}

	//-------------------------------------------------------
	// Paramètres concernant les résultats
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
	} else {
		homemade_error_msg("Missing the results folder path!");
	}

	if (field_parser.search(1, "ResultsFileName")) {
		input_params.results_file_name = field_parser.next(
				input_params.results_file_name);
		std::cout << input_params.results_file_name << std::endl;
	} else {
		input_params.results_file_name = "out_";
	}

	}

};