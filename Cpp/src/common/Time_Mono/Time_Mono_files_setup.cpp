#include "Time_Mono_files_setup.h"

namespace carl
{
    
void Time_Mono_Files_Setup::print_time_mono_iterate_params(const std::string& output_filename)
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");

	std::ofstream output_file(output_filename);

	output_file << "SchedulerType " << carl::ClusterSchedulerType_to_string(m_input_params.scheduler) << std::endl;
	output_file << "ScratchFolderPath "  << m_input_params.scratch_folder_path 						  << std::endl;
	output_file << "SolutionsFolderPath " << m_input_params.feti_solution_path  						  << std::endl;
	output_file << "ResultsFolderPath "   << m_input_params.results_folder_path 						  << std::endl;
	output_file << "ResultsFileName "     << m_input_params.results_file_name   						  << std::endl;
	output_file << "TimeMeshFile "        << m_input_params.time_mesh_file      						  << std::endl;

	output_file.close();
	
}

void Time_Mono_Files_Setup::print_PBS_script(const std::string& output_filename, const std::string& job_name, const std::string& output_name, const std::string& error_name, const std::string& common_script, const std::string& command_to_run)
{
	std::ofstream output_script(output_filename);
	output_script << "#!/bin/bash" << std::endl;
	output_script << std::endl;
	output_script << "#PBS -S /bin/bash" << std::endl;
	output_script << "#PBS -N " << job_name << std::endl;
	output_script << "#PBS -o " << output_name << std::endl;
	output_script << "#PBS -e " << error_name << std::endl;
	output_script << common_script << std::endl;
	output_script << command_to_run << std::endl;
	output_script.close();
};

void Time_Mono_Files_Setup::set_Time_Mono_input_parameters(time_mono_setup_params& input_params)
{
	m_bInputParamsSet = true;
	m_input_params = input_params;
}

void Time_Mono_Files_Setup::set_scratch_folder()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	
	if(m_comm.rank() == 0)
	{
		std::string command_string;

		command_string = "rm -rf " + m_input_params.scratch_folder_path;
		carl::exec_command(command_string.c_str());

		command_string = "mkdir -p " + m_input_params.scratch_folder_path;
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;
	}

	m_bScratchFolderExists = true;

}

void Time_Mono_Files_Setup::generate_state_file()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

	if(m_comm.rank() == 0)
	{
		std::string state_file_path = m_input_params.scratch_folder_path + "/Time_Mono_state.dat";
		
		unsigned int it_init = 0;
		std::ofstream state_file;

		state_file.open(state_file_path);
		state_file.precision(15);
		state_file << it_init << std::endl;
		state_file.close();

	}

	m_bStateFileExists = true;

}

void Time_Mono_Files_Setup::generate_domain_state_folders()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

	if(m_comm.rank() == 0)
	{
		std::string command_string;

		command_string = "mkdir -p " + m_input_params.scratch_folder_path + "/domain_A_state";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;

		command_string = "mkdir -p " + m_input_params.scratch_folder_path + "/domain_B_state";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;
	}

	m_bDomainStateFoldersExists = true;

}

void Time_Mono_Files_Setup::set_results_folder()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");

	if(m_comm.rank() == 0)
	{
		std::string command_string;

		command_string = "rm -rf " + m_input_params.results_folder_path;
		carl::exec_command(command_string.c_str());

		command_string = "mkdir -p " + m_input_params.results_folder_path;
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;
	}

	m_bResultsFolderExists = true;

}

void Time_Mono_Files_Setup::copy_initial_conditions_files()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bResultsFolderExists,"Results folder not set yet!");

	if(m_comm.rank() == 0)
	{
		std::string command_string;

		// Copy the initial conditions files to the domain state folders in the scratch folder

		command_string = "cp " + m_input_params.domain_A_init_cond_disp_file + " " + m_input_params.scratch_folder_path + "/domain_A_state/";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;

		command_string = "cp " + m_input_params.domain_A_init_cond_vel_file + " " + m_input_params.scratch_folder_path + "/domain_A_state/";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;

		command_string = "cp " + m_input_params.domain_A_init_cond_acc_file + " " + m_input_params.scratch_folder_path + "/domain_A_state/";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;

		command_string = "cp " + m_input_params.domain_B_init_cond_disp_file + " " + m_input_params.scratch_folder_path + "/domain_B_state/";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;

		command_string = "cp " + m_input_params.domain_B_init_cond_vel_file + " " + m_input_params.scratch_folder_path + "/domain_B_state/";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;

		command_string = "cp " + m_input_params.domain_B_init_cond_acc_file + " " + m_input_params.scratch_folder_path + "/domain_B_state/";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;

		// Copy the initial displacements in the results folder with the good name

		command_string = "cp " + m_input_params.domain_A_init_cond_disp_file + " " + m_input_params.results_folder_path + "/" + m_input_params.results_file_name + "A_0.petscvec";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;

		command_string = "cp " + m_input_params.domain_B_init_cond_disp_file + " " + m_input_params.results_folder_path + "/" + m_input_params.results_file_name + "B_0.petscvec";
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;
		
	}

}

void Time_Mono_Files_Setup::generate_libmesh_external_solver_assemble_rhs_inputs()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

	m_ext_solver_assemble_rhs_A_input_filename = m_input_params.scratch_folder_path + "/assemble_dyn_rhs_A.txt";
	m_ext_solver_assemble_rhs_B_input_filename = m_input_params.scratch_folder_path + "/assemble_dyn_rhs_B.txt";
	
	if(m_comm.rank() == 0)
	{
		std::ofstream output_file_A(m_ext_solver_assemble_rhs_A_input_filename);
		output_file_A << "ScratchFolderPath " << m_input_params.scratch_folder_path       << std::endl;
		output_file_A << "TimeMeshFile "      << m_input_params.time_mesh_file            << std::endl;
		output_file_A << "Mesh "              << m_input_params.space_mesh_A_file         << std::endl;
		output_file_A << "NewmarkParameters " << m_input_params.newmark_parameters_A_file << std::endl;
		output_file_A << "SysMatrix "         << m_input_params.domain_A_matrix_file      << std::endl;
		output_file_A << "SysMatrix_M "       << m_input_params.domain_A_M_matrix_file    << std::endl;
		output_file_A << "SysMatrix_K "       << m_input_params.domain_A_K_matrix_file    << std::endl;
		output_file_A << "SysMatrix_C "       << m_input_params.domain_A_C_matrix_file    << std::endl;
		output_file_A.close();

		std::ofstream output_file_B(m_ext_solver_assemble_rhs_B_input_filename);
		output_file_B << "ScratchFolderPath " << m_input_params.scratch_folder_path       << std::endl;
		output_file_B << "TimeMeshFile "      << m_input_params.time_mesh_file            << std::endl;
		output_file_B << "Mesh "              << m_input_params.space_mesh_B_file         << std::endl;
		output_file_B << "NewmarkParameters " << m_input_params.newmark_parameters_B_file << std::endl;
		output_file_B << "SysMatrix "         << m_input_params.domain_B_matrix_file      << std::endl;
		output_file_B << "SysMatrix_M "       << m_input_params.domain_B_M_matrix_file    << std::endl;
		output_file_B << "SysMatrix_K "       << m_input_params.domain_B_K_matrix_file    << std::endl;
		output_file_B << "SysMatrix_C "       << m_input_params.domain_B_C_matrix_file    << std::endl;
		output_file_B.close();
	}

	m_bSetExternalSolversRhsInputFiles = true;
}

void Time_Mono_Files_Setup::generate_libmesh_external_solver_assemble_rhs_scripts()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversRhsInputFiles,"External solver rhs assemble input files not set yet!");

	m_ext_solver_assemble_rhs_A_script_filename = m_input_params.scratch_folder_path + "/ext_solver_rhs_assemble_A.sh";
	m_ext_solver_assemble_rhs_B_script_filename = m_input_params.scratch_folder_path + "/ext_solver_rhs_assemble_B.sh";

	switch (m_input_params.scheduler)
	{
		case ClusterSchedulerType::LOCAL :	this->generate_libmesh_external_solver_assemble_rhs_scripts_LOCAL();
						break;

		case ClusterSchedulerType::PBS :    this->generate_libmesh_external_solver_assemble_rhs_scripts_PBS();
						break;

		case ClusterSchedulerType::SLURM :	homemade_error_msg("Scheduler SLURM not implemented yet!");
						break;
		default : homemade_error_msg("Invalid scheduler name!");
	}

	m_bSetExternalSolversRhsScriptFiles = true;
}

void Time_Mono_Files_Setup::generate_libmesh_external_solver_assemble_rhs_scripts_LOCAL()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	
	if(m_comm.rank() == 0)
	{
		std::string command_to_run;

		command_to_run = m_input_params.ext_solver_BIG_rhs_assembly + " " + m_ext_solver_assemble_rhs_A_input_filename;
		std::ofstream output_script(m_ext_solver_assemble_rhs_A_script_filename);
		output_script << command_to_run << std::endl;
		output_script.close();

		command_to_run = m_input_params.ext_solver_micro_rhs_assembly + " " + m_ext_solver_assemble_rhs_B_input_filename;
		output_script.open(m_ext_solver_assemble_rhs_B_script_filename);
		output_script << command_to_run << std::endl;
		output_script.close();
	}
}

void Time_Mono_Files_Setup::generate_libmesh_external_solver_assemble_rhs_scripts_PBS()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	
	if(m_comm.rank() == 0)
	{
		// Get the full common script file into a string
		std::ifstream base_script(m_input_params.script_filename);
		std::string common_script((std::istreambuf_iterator<char>(base_script)),
									std::istreambuf_iterator<char>());
		base_script.close();

		std::string pbs_output;
		std::string pbs_error;
		std::string command_to_run;

		pbs_output = m_input_params.scratch_folder_path + "/output_ext_rhs_assemble_A.txt";
		pbs_error = m_input_params.scratch_folder_path + "/error_ext_rhs_assemble_A.txt";
		command_to_run = m_input_params.ext_solver_BIG_rhs_assembly + " " + m_ext_solver_assemble_rhs_A_input_filename;
		
		this->print_PBS_script(	m_ext_solver_assemble_rhs_A_script_filename, "ext_rhs_assemble_A",
							pbs_output, pbs_error, common_script,
							command_to_run);

		pbs_output = m_input_params.scratch_folder_path + "/output_ext_rhs_assemble_B.txt";
		pbs_error = m_input_params.scratch_folder_path + "/error_ext_rhs_assemble_B.txt";
		command_to_run = m_input_params.ext_solver_micro_rhs_assembly + " " + m_ext_solver_assemble_rhs_B_input_filename;
		
		this->print_PBS_script(	m_ext_solver_assemble_rhs_B_script_filename, "ext_rhs_assemble_B",
							pbs_output, pbs_error, common_script,
							command_to_run);

	}
}

void Time_Mono_Files_Setup::generate_libmesh_external_solver_update_init_cond_inputs()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

	m_ext_solver_update_init_cond_A_input_filename = m_input_params.scratch_folder_path + "/update_init_cond_A.txt";
	m_ext_solver_update_init_cond_B_input_filename = m_input_params.scratch_folder_path + "/update_init_cond_B.txt";
	
	if(m_comm.rank() == 0)
	{
		std::ofstream output_file_A(m_ext_solver_update_init_cond_A_input_filename);
		output_file_A << "TimeMeshFile "      << m_input_params.time_mesh_file                                   				  << std::endl;
		output_file_A << "SpaceMeshFile "     << m_input_params.space_mesh_A_file                                                 << std::endl;
		output_file_A << "NewmarkParameters " << m_input_params.newmark_parameters_A_file                                         << std::endl;
		output_file_A << "Solution "          << m_input_params.feti_solution_path + "/coupled_sol_A.petscvec"                    << std::endl;
		output_file_A << "InitCond_disp "     << m_input_params.scratch_folder_path + "/domain_A_state/init_cond_A_disp.petscvec" << std::endl;
		output_file_A << "InitCond_vel "      << m_input_params.scratch_folder_path + "/domain_A_state/init_cond_A_vel.petscvec"  << std::endl;
		output_file_A << "InitCond_acc "      << m_input_params.scratch_folder_path + "/domain_A_state/init_cond_A_acc.petscvec"  << std::endl;
		output_file_A.close();

		std::ofstream output_file_B(m_ext_solver_update_init_cond_B_input_filename);
		output_file_B << "TimeMeshFile "      << m_input_params.time_mesh_file                                   				  << std::endl;
		output_file_B << "SpaceMeshFile "     << m_input_params.space_mesh_B_file                                                 << std::endl;
		output_file_B << "NewmarkParameters " << m_input_params.newmark_parameters_B_file                                         << std::endl;
		output_file_B << "Solution "          << m_input_params.feti_solution_path + "/coupled_sol_B.petscvec"                    << std::endl;
		output_file_B << "InitCond_disp "     << m_input_params.scratch_folder_path + "/domain_B_state/init_cond_B_disp.petscvec" << std::endl;
		output_file_B << "InitCond_vel "      << m_input_params.scratch_folder_path + "/domain_B_state/init_cond_B_vel.petscvec"  << std::endl;
		output_file_B << "InitCond_acc "      << m_input_params.scratch_folder_path + "/domain_B_state/init_cond_B_acc.petscvec"  << std::endl;
		output_file_B.close();
	}

	m_bSetExternalSolversUpdateInitCondInputFiles = true;
}

void Time_Mono_Files_Setup::generate_libmesh_external_solver_update_init_cond_scripts()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversUpdateInitCondInputFiles,"External solver update init cond input files not set yet!");

	m_ext_solver_update_init_cond_A_script_filename = m_input_params.scratch_folder_path + "/ext_solver_update_init_cond_A.sh";
	m_ext_solver_update_init_cond_B_script_filename = m_input_params.scratch_folder_path + "/ext_solver_update_init_cond_B.sh";

	switch (m_input_params.scheduler)
	{
		case ClusterSchedulerType::LOCAL :	this->generate_libmesh_external_solver_update_init_cond_scripts_LOCAL();
						break;

		case ClusterSchedulerType::PBS :    this->generate_libmesh_external_solver_update_init_cond_scripts_PBS();
						break;

		case ClusterSchedulerType::SLURM :	homemade_error_msg("Scheduler SLURM not implemented yet!");
						break;
		default : homemade_error_msg("Invalid scheduler name!");
	}

	m_bSetExternalSolversUpdateInitCondScriptFiles = true;
}

void Time_Mono_Files_Setup::generate_libmesh_external_solver_update_init_cond_scripts_LOCAL()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	
	if(m_comm.rank() == 0)
	{
		std::string command_to_run;

		command_to_run = m_input_params.ext_solver_BIG_update_init_cond + " " + m_ext_solver_update_init_cond_A_input_filename;
		std::ofstream output_script(m_ext_solver_update_init_cond_A_script_filename);
		output_script << command_to_run << std::endl;
		output_script.close();

		command_to_run = m_input_params.ext_solver_micro_update_init_cond + " " + m_ext_solver_update_init_cond_B_input_filename;
		output_script.open(m_ext_solver_update_init_cond_B_script_filename);
		output_script << command_to_run << std::endl;
		output_script.close();
	}

}

void Time_Mono_Files_Setup::generate_libmesh_external_solver_update_init_cond_scripts_PBS()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	
	if(m_comm.rank() == 0)
	{
		// Get the full common script file into a string
		std::ifstream base_script(m_input_params.script_filename);
		std::string common_script((std::istreambuf_iterator<char>(base_script)),
									std::istreambuf_iterator<char>());
		base_script.close();

		std::string pbs_output;
		std::string pbs_error;
		std::string command_to_run;

		pbs_output = m_input_params.scratch_folder_path + "/output_ext_update_init_cond_A.txt";
		pbs_error = m_input_params.scratch_folder_path + "/error_ext_update_init_cond_A.txt";
		command_to_run = m_input_params.ext_solver_BIG_update_init_cond + " " + m_ext_solver_update_init_cond_A_input_filename;
		
		this->print_PBS_script(	m_ext_solver_update_init_cond_A_script_filename, "ext_update_init_cond_A",
							pbs_output, pbs_error, common_script,
							command_to_run);

		pbs_output = m_input_params.scratch_folder_path + "/output_ext_update_init_cond_B.txt";
		pbs_error = m_input_params.scratch_folder_path + "/error_ext_update_init_cond_B.txt";
		command_to_run = m_input_params.ext_solver_micro_update_init_cond + " " + m_ext_solver_update_init_cond_B_input_filename;
		
		this->print_PBS_script(	m_ext_solver_update_init_cond_B_script_filename, "ext_update_init_cond_B",
							pbs_output, pbs_error, common_script,
							command_to_run);
	}

}

void Time_Mono_Files_Setup::generate_FETI_setup_init_scripts()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

	m_FETI_setup_init_script_filename = m_input_params.scratch_folder_path + "/CArl_FETI_setup_init.sh";

	switch (m_input_params.scheduler)
	{
		case ClusterSchedulerType::LOCAL :	this->generate_FETI_setup_init_scripts_LOCAL();
						break;

		case ClusterSchedulerType::PBS :    this->generate_FETI_setup_init_scripts_PBS();
						break;

		case ClusterSchedulerType::SLURM :	homemade_error_msg("Scheduler SLURM not implemented yet!");
						break;
		default : homemade_error_msg("Invalid scheduler name!");
	}

	m_bSetFETISetupInitScriptFile = true;
}

void Time_Mono_Files_Setup::generate_FETI_setup_init_scripts_LOCAL()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	
	if(m_comm.rank() == 0)
	{
		std::string command_to_run;

		command_to_run = "./CArl_FETI_setup_init -i " + m_input_params.feti_setup_params_file;
		std::ofstream output_script(m_FETI_setup_init_script_filename);
		output_script << "mpirun -n " << m_comm.size() << " " << command_to_run << std::endl;
		output_script.close();

	}
	
}

void Time_Mono_Files_Setup::generate_FETI_setup_init_scripts_PBS()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	
	if(m_comm.rank() == 0)
	{
		// Get the full common script file into a string
		std::ifstream base_script(m_input_params.script_filename);
		std::string common_script((std::istreambuf_iterator<char>(base_script)),
									std::istreambuf_iterator<char>());
		base_script.close();

		std::string pbs_output;
		std::string pbs_error;
		std::string command_to_run;

		pbs_output = m_input_params.scratch_folder_path + "/output_CArl_FETI_setup_init.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_CArl_FETI_setup_init.txt";
		command_to_run = "mpirun -n " + std::to_string(m_comm.size()) + " ./CArl_FETI_setup_init -i " +
							m_input_params.feti_setup_params_file;

		this->print_PBS_script(	m_FETI_setup_init_script_filename, "CArl_setup_init",
							pbs_output, pbs_error, common_script,
							command_to_run);	
	}
	
}

void Time_Mono_Files_Setup::generate_Time_Mono_inputs()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

	// Pour l'exécutable CArl_Time_Mono_iterate_finish.
	m_CArl_Time_Mono_iterate_finish_input_filename = m_input_params.scratch_folder_path + "/CArl_Time_Mono_iterate_finish.txt";
	if(m_comm.rank() == 0)
	{
		this->print_time_mono_iterate_params(m_CArl_Time_Mono_iterate_finish_input_filename);
	}
	
	m_bSetCArlTimeMonoInputs = true;
}

void Time_Mono_Files_Setup::generate_Time_Mono_scripts()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetCArlTimeMonoInputs,"CArl_Time_Mono input files not set yet!");

	m_CArl_Time_Mono_iterate_finish_script_filename = m_input_params.scratch_folder_path + "/CArl_Time_Mono_iterate_finish.sh";

	switch (m_input_params.scheduler)
	{
		case ClusterSchedulerType::LOCAL :	this->generate_Time_Mono_scripts_LOCAL();
						break;

		case ClusterSchedulerType::PBS :    this->generate_Time_Mono_scripts_PBS();
						break;

		case ClusterSchedulerType::SLURM :	homemade_error_msg("Scheduler SLURM not implemented yet!");
						break;
		default : homemade_error_msg("Invalid scheduler name!");
	}

	m_bSetCArlTimeMonoScripts = true;
}

void Time_Mono_Files_Setup::generate_Time_Mono_scripts_LOCAL()
{
	
	if(m_comm.rank() == 0)
	{
		std::string command_to_run;

		command_to_run = "./CArl_Time_Mono_iterate_finish -i " + m_CArl_Time_Mono_iterate_finish_input_filename;
		std::ofstream output_script(m_CArl_Time_Mono_iterate_finish_script_filename);
		output_script << "mpirun -n " << m_comm.size() << " " << command_to_run << std::endl;
		output_script.close();
		
	}

}

void Time_Mono_Files_Setup::generate_Time_Mono_scripts_PBS()
{
	
	if(m_comm.rank() == 0)
	{
		// Get the full common script file into a string
		std::ifstream base_script(m_input_params.script_filename);
		std::string common_script((std::istreambuf_iterator<char>(base_script)),
									std::istreambuf_iterator<char>());
		base_script.close();

		std::string pbs_output;
		std::string pbs_error;
		std::string command_to_run;

		pbs_output = m_input_params.scratch_folder_path + "/output_CArl_Time_Mono_iterate_finish.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_CArl_Time_Mono_iterate_finish.txt";
		command_to_run = "mpirun -n " + std::to_string(m_comm.size()) + " ./CArl_Time_Mono_iterate_finish -i " +
							m_CArl_Time_Mono_iterate_finish_input_filename;

		this->print_PBS_script(	m_CArl_Time_Mono_iterate_finish_script_filename, "CArl_Time_Mono_iterate_finish",
							pbs_output, pbs_error, common_script,
							command_to_run);	
	}

}

void Time_Mono_Files_Setup::generate_Time_Mono_launch_scripts()
{

	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversRhsScriptFiles,"External solver rhs assembly scripts not set yet!");
	homemade_assert_msg(m_bSetCArlTimeMonoScripts,"CArl_Time_Mono script not set yet!");

	m_Time_Mono_iterate_init_script_filename   = m_input_params.scratch_folder_path + "/Time_Mono_iterate_init.sh";
	m_Time_Mono_iterate_finish_script_filename = m_input_params.scratch_folder_path + "/Time_Mono_iterate_finish.sh";
	
	switch (m_input_params.scheduler)
	{
		case ClusterSchedulerType::LOCAL :	this->generate_Time_Mono_launch_scripts_LOCAL();
						break;

		case ClusterSchedulerType::PBS :    this->generate_Time_Mono_launch_scripts_PBS();
						break;

		case ClusterSchedulerType::SLURM :	homemade_error_msg("Scheduler SLURM not implemented yet!");
						break;
		default : homemade_error_msg("Invalid scheduler name!");
	}

	m_bSetCArlTimeMonoLaunchScripts = true;
}

void Time_Mono_Files_Setup::generate_Time_Mono_launch_scripts_LOCAL()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversRhsScriptFiles,"External solver rhs assembly scripts not set yet!");
	homemade_assert_msg(m_bSetCArlTimeMonoScripts,"CArl_Time_Mono script not set yet!");

	// ONLY WORK ON PROCESSOR 0 !!!
	if(m_comm.rank() == 0)
	{
		// Création de Time_Mono_iterate_init.sh
		std::ofstream Time_Mono_iterate_init_script(m_Time_Mono_iterate_init_script_filename);
		Time_Mono_iterate_init_script << "#!/bin/bash" << std::endl;
		Time_Mono_iterate_init_script << std::endl;
		Time_Mono_iterate_init_script << ". " << m_ext_solver_assemble_rhs_A_script_filename << std::endl;
		Time_Mono_iterate_init_script << ". " << m_ext_solver_assemble_rhs_B_script_filename << std::endl;
		Time_Mono_iterate_init_script << ". " << m_FETI_setup_init_script_filename << std::endl;
		Time_Mono_iterate_init_script.close();

		// Création de Time_Mono_iterate_finish.sh
		std::ofstream Time_Mono_iterate_finish_script(m_Time_Mono_iterate_finish_script_filename);
		Time_Mono_iterate_finish_script << "#!/bin/bash" << std::endl;
		Time_Mono_iterate_finish_script << std::endl;
		Time_Mono_iterate_finish_script << ". " << m_ext_solver_update_init_cond_A_script_filename << std::endl;
		Time_Mono_iterate_finish_script << ". " << m_ext_solver_update_init_cond_B_script_filename << std::endl;
		Time_Mono_iterate_finish_script << ". " << m_CArl_Time_Mono_iterate_finish_script_filename << std::endl;
		Time_Mono_iterate_finish_script.close();
		
	}
}

void Time_Mono_Files_Setup::generate_Time_Mono_launch_scripts_PBS()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversRhsScriptFiles,"External solver rhs assembly scripts not set yet!");
	homemade_assert_msg(m_bSetCArlTimeMonoScripts,"CArl_Time_Mono script not set yet!");

	// ONLY WORK ON PROCESSOR 0 !!!
	if(m_comm.rank() == 0)
	{
		// Création de Time_Mono_iterate_init.sh
		std::ofstream Time_Mono_iterate_init_script(m_Time_Mono_iterate_init_script_filename);
		Time_Mono_iterate_init_script << "#!/bin/bash" << std::endl;
		Time_Mono_iterate_init_script << std::endl;
		Time_Mono_iterate_init_script << "job1_A=`qsub " << m_ext_solver_assemble_rhs_A_script_filename << "`" << std::endl;
		Time_Mono_iterate_init_script << "job1_B=`qsub " << m_ext_solver_assemble_rhs_B_script_filename << "`" << std::endl;
		Time_Mono_iterate_init_script << "job2=`qsub -W depend=afterok:$job1_A:$job1_B " 
										<< m_FETI_setup_init_script_filename << "`" << std::endl;
		Time_Mono_iterate_init_script.close();

		// Création de Time_Mono_iterate_finish.sh
		std::ofstream Time_Mono_iterate_finish_script(m_Time_Mono_iterate_finish_script_filename);
		Time_Mono_iterate_finish_script << "#!/bin/bash" << std::endl;
		Time_Mono_iterate_finish_script << std::endl;
		Time_Mono_iterate_finish_script << "job3_A=`qsub " << m_ext_solver_update_init_cond_A_script_filename << "`" << std::endl;
		Time_Mono_iterate_finish_script << "job3_B=`qsub " << m_ext_solver_update_init_cond_B_script_filename << "`" << std::endl;
		Time_Mono_iterate_finish_script << "job4=`qsub -W depend=afterok:$job3_A:$job3_B " 
										  << m_CArl_Time_Mono_iterate_finish_script_filename << "`" << std::endl;
		Time_Mono_iterate_finish_script.close();
	}
}

}