#include "solver_files_setup.h"

namespace carl
{
void Solver_Files_Setup::print_feti_setup_finish_params(const std::string& output_filename)
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");

	std::ofstream output_file(output_filename);
	output_file << "ClusterSchedulerType " << carl::ClusterSchedulerType_to_string(m_input_params.scheduler) << std::endl;
	output_file << "ScratchFolderPath " << m_input_params.scratch_folder_path << std::endl;
	output_file << "CouplingMatricesFolder " << m_input_params.coupling_folder_path << std::endl;
	output_file << "CGPreconditionerType " << carl::BaseCGPrecondType_to_string(m_input_params.CG_precond_type) << std::endl;

	if(m_input_params.bUseRigidBodyModes)
	{
		output_file << "UseRigidBodyModesB" << std::endl;
		output_file << "RBVectorBase " << m_input_params.RB_vectors_base << std::endl;
		output_file << "NbOfRBVectors " << m_input_params.nb_of_rb_vectors << std::endl;
	}
	output_file.close();
};

void Solver_Files_Setup::print_feti_iterate_params(const std::string& output_filename)
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");

	std::ofstream output_file(output_filename);
	output_file << "ClusterSchedulerType " << carl::ClusterSchedulerType_to_string(m_input_params.scheduler) << std::endl;
	output_file << "ScratchFolderPath " << m_input_params.scratch_folder_path << std::endl;
	output_file << "CouplingMatricesFolder " << m_input_params.coupling_folder_path << std::endl;

	if(m_input_params.bUseRigidBodyModes)
	{
		output_file << "UseRigidBodyModesB" << std::endl;
		output_file << "RBVectorBase " << m_input_params.RB_vectors_base << std::endl;
		output_file << "NbOfRBVectors " << m_input_params.nb_of_rb_vectors << std::endl;
	}

	output_file << "CGPreconditionerType " << carl::BaseCGPrecondType_to_string(m_input_params.CG_precond_type) << std::endl;
	output_file << "CoupledConvAbs " << m_input_params.CG_coupled_conv_abs << std::endl;
	output_file << "CoupledConvRel " << m_input_params.CG_coupled_conv_rel << std::endl;
	output_file << "CoupledCorrConvRel " << m_input_params.CG_coupled_conv_corr << std::endl;
	output_file << "CoupledDiv " << m_input_params.CG_coupled_div << std::endl;
	output_file << "CoupledIterMax " << m_input_params.CG_coupled_conv_max << std::endl;
	output_file.close();
};

void Solver_Files_Setup::print_feti_solution_params(const std::string& output_filename)
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");

	std::ofstream output_file(output_filename);
	output_file << "ScratchFolderPath " << m_input_params.scratch_folder_path << std::endl;
	output_file << "OutputFolder " << m_input_params.output_folder << std::endl;

	if(m_input_params.bUseRigidBodyModes)
	{
		output_file << "UseRigidBodyModesB" << std::endl;
	}
	output_file.close();
};

void Solver_Files_Setup::print_PBS_script(const std::string& output_filename, const std::string& job_name, const std::string& output_name, const std::string& error_name, const std::string& common_script, const std::string& command_to_run)
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

void Solver_Files_Setup::set_FETI_input_parameters(feti_setup_init_params& input_params)
{
	m_bInputParamsSet = true;
	m_input_params = input_params;
};

void Solver_Files_Setup::set_scratch_folder()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");

	m_bScratchFolderExists = true;
	
	if(m_comm.rank() == 0)
	{
		std::string command_string;

		command_string = "rm -rf " + m_input_params.scratch_folder_path;
		carl::exec_command(command_string.c_str());

		command_string = "mkdir -p " + m_input_params.scratch_folder_path;
		carl::exec_command(command_string.c_str());
		std::cout << command_string << std::endl;
	}
};

void Solver_Files_Setup::generate_libmesh_external_solver_inputs()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

	m_ext_solver_u0_A_input_filename = m_input_params.scratch_folder_path + "/ext_solver_u0_A.txt";
	m_ext_solver_u0_B_input_filename = m_input_params.scratch_folder_path + "/ext_solver_u0_B.txt";

	m_ext_solver_A_input_filename = m_input_params.scratch_folder_path + "/ext_solver_A.txt";
	m_ext_solver_B_input_filename = m_input_params.scratch_folder_path + "/ext_solver_B.txt";

	GetPot field_parser_A, field_parser_B;

	// Get the first input file
	field_parser_A.parse_input_file(m_input_params.ext_solver_BIG_input, "#", "\n", " \t\n");
	field_parser_B.parse_input_file(m_input_params.ext_solver_micro_input, "#", "\n", " \t\n");

	carl::libmesh_solve_linear_system_input_params solver_A_input_params, solver_B_input_params;

	carl::get_input_params(field_parser_A, solver_A_input_params);
	carl::get_input_params(field_parser_B, solver_B_input_params);

	// Set K_i * u_0,i  = F_i
	solver_A_input_params.output_base = m_input_params.scratch_folder_path + "/ext_solver_u0_A";
	solver_B_input_params.output_base = m_input_params.scratch_folder_path + "/ext_solver_u0_B";

	carl::print_input_params(m_ext_solver_u0_A_input_filename,solver_A_input_params);
	carl::print_input_params(m_ext_solver_u0_B_input_filename,solver_B_input_params);

	// Set K_i * y(k)_i = C_i^T * p(k) and K_i * x_i(kkk)  = C_i^T * phi(kkk)
	solver_A_input_params.output_base = m_input_params.scratch_folder_path + "/ext_solver_A";
	solver_B_input_params.output_base = m_input_params.scratch_folder_path + "/ext_solver_B";

	solver_A_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/ext_solver_A_rhs.petscvec";
	solver_B_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/ext_solver_B_rhs.petscvec";

	carl::print_input_params(m_ext_solver_A_input_filename,solver_A_input_params);
	carl::print_input_params(m_ext_solver_B_input_filename,solver_B_input_params);

	m_bSetExternalSolversInputFiles = true;
}

void Solver_Files_Setup::generate_libmesh_external_solver_scripts()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");

	m_ext_solver_u0_A_script_filename = m_input_params.scratch_folder_path + "/ext_solver_u0_A.sh";
	m_ext_solver_u0_B_script_filename = m_input_params.scratch_folder_path + "/ext_solver_u0_B.sh";

	m_ext_solver_A_script_filename = m_input_params.scratch_folder_path + "/ext_solver_A.sh";
	m_ext_solver_B_script_filename = m_input_params.scratch_folder_path + "/ext_solver_B.sh";

	switch (m_input_params.scheduler)
	{
		case ClusterSchedulerType::LOCAL :	this->generate_libmesh_external_solver_scripts_LOCAL();
						break;

		case ClusterSchedulerType::PBS :    this->generate_libmesh_external_solver_scripts_PBS();
						break;

		case ClusterSchedulerType::SLURM :	homemade_error_msg("Scheduler not implemented yet!");
						break;
		default : homemade_error_msg("Invalid scheduler name!");
	}

	m_bSetExternalSolversFiles = true;
}

void Solver_Files_Setup::generate_libmesh_external_solver_scripts_LOCAL()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");

	if(m_comm.rank() == 0)
	{
		std::string command_to_run;

		// Set the u0_i scripts
		command_to_run = m_input_params.ext_solver_BIG + " " + m_ext_solver_u0_A_input_filename;
		std::ofstream output_script(m_ext_solver_u0_A_script_filename);
		output_script << command_to_run << std::endl;
		output_script.close();

		command_to_run = m_input_params.ext_solver_micro + " " + m_ext_solver_u0_B_input_filename;
		output_script.open(m_ext_solver_u0_B_script_filename);
		output_script << command_to_run << std::endl;
		output_script.close();

		// Set the yk_i scripts
		command_to_run = m_input_params.ext_solver_BIG + " " + m_ext_solver_A_input_filename;
		output_script.open(m_ext_solver_A_script_filename);
		output_script << command_to_run << std::endl;
		output_script.close();

		command_to_run = m_input_params.ext_solver_micro + " " + m_ext_solver_B_input_filename;
		output_script.open(m_ext_solver_B_script_filename);
		output_script << command_to_run << std::endl;
		output_script.close();
	}
}

void Solver_Files_Setup::generate_libmesh_external_solver_scripts_PBS()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");

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

		// Set the u0_i scripts
		pbs_output = m_input_params.scratch_folder_path + "/output_ext_u0_A.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_ext_u0_A.txt";
		command_to_run = m_input_params.ext_solver_BIG + " " + m_ext_solver_u0_A_input_filename;

		this->print_PBS_script(	m_ext_solver_u0_A_script_filename, "ext_u0_A",
							pbs_output, pbs_error, common_script,
							command_to_run);

		pbs_output = m_input_params.scratch_folder_path + "/output_ext_u0_B.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_ext_u0_B.txt";
		command_to_run = m_input_params.ext_solver_micro + " " + m_ext_solver_u0_B_input_filename;

		this->print_PBS_script(	m_ext_solver_u0_B_script_filename, "ext_u0_B",
							pbs_output, pbs_error,common_script,
							command_to_run);

		// Set the yk_i scripts
		pbs_output = m_input_params.scratch_folder_path + "/output_ext_A.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_ext_A.txt";
		command_to_run = m_input_params.ext_solver_BIG + " " + m_ext_solver_A_input_filename;

		this->print_PBS_script(	m_ext_solver_A_script_filename, "ext_A",
							pbs_output, pbs_error,common_script,
							command_to_run);

		pbs_output = m_input_params.scratch_folder_path + "/output_ext_B.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_ext_B.txt";
		command_to_run = m_input_params.ext_solver_micro + " " + m_ext_solver_B_input_filename;

		this->print_PBS_script(	m_ext_solver_B_script_filename, "ext_B",
							pbs_output, pbs_error,common_script,
							command_to_run);
	}
}

void Solver_Files_Setup::generate_FETI_inputs()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

	m_CArl_FETI_setup_finish_input_filename = m_input_params.scratch_folder_path + "/CArl_FETI_setup_finish.txt";
	m_CArl_FETI_iterate_input_filename = m_input_params.scratch_folder_path + "/CArl_FETI_iterate.txt";
	m_CArl_FETI_solution_input_filename = m_input_params.scratch_folder_path + "/CArl_FETI_set_sol.txt";

	if(m_comm.rank() == 0)
	{
		this->print_feti_setup_finish_params(m_CArl_FETI_setup_finish_input_filename);
		this->print_feti_iterate_params(m_CArl_FETI_iterate_input_filename);
		this->print_feti_solution_params(m_CArl_FETI_solution_input_filename);
	}

	m_bSetCArlFETIInputs = true;
}

void Solver_Files_Setup::generate_FETI_scripts()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetCArlFETIInputs,"CArl_FETI input files not set yet!");

	if(m_comm.rank() == 0)
	{
		std::string command_string = "mkdir -p " + m_input_params.output_folder;
		carl::exec_command(command_string.c_str());
	}

	m_CArl_FETI_setup_finish_script_filename = m_input_params.scratch_folder_path + "/FETI_setup_finish.sh";
	m_CArl_FETI_iterate_script_filename = m_input_params.scratch_folder_path + "/FETI_iterate.sh";
	m_CArl_FETI_solution_script_filename = m_input_params.scratch_folder_path + "/FETI_solution.sh";

	switch (m_input_params.scheduler)
	{
		case ClusterSchedulerType::LOCAL :	this->generate_FETI_scripts_LOCAL();
						break;

		case ClusterSchedulerType::PBS :    this->generate_FETI_scripts_PBS();
						break;

		case ClusterSchedulerType::SLURM :	homemade_error_msg("Scheduler not implemented yet!");
						break;
		default : homemade_error_msg("Invalid scheduler name!");
	}

	m_bSetCArlFETIScripts = true;
}

void Solver_Files_Setup::generate_FETI_scripts_LOCAL()
{
	if(m_comm.rank() == 0)
	{
		std::string command_to_run;

		// Set the u0_i scripts
		command_to_run = "./CArl_FETI_setup_finish -i " + m_CArl_FETI_setup_finish_input_filename;
		std::ofstream output_script(m_CArl_FETI_setup_finish_script_filename);
		output_script << "mpirun -n " << m_comm.size() << " " << command_to_run << std::endl;
		output_script.close();
		
		
		command_to_run = "./CArl_FETI_iterate -i " + m_CArl_FETI_iterate_input_filename;
		output_script.open(m_CArl_FETI_iterate_script_filename);
		output_script << "mpirun -n " << m_comm.size() << " "  << command_to_run << std::endl;
		output_script.close();

		command_to_run = "./CArl_FETI_solution -i " + m_CArl_FETI_solution_input_filename;
		output_script.open(m_CArl_FETI_solution_script_filename);
		output_script << "mpirun -n " << m_comm.size() << " "  << command_to_run << std::endl;
		output_script.close();
	}
}

void Solver_Files_Setup::generate_FETI_scripts_PBS()
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

		// Set the u0_i scripts
		pbs_output = m_input_params.scratch_folder_path + "/output_CArl_FETI_setup_finish.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_CArl_FETI_setup_finish.txt";
		command_to_run = "mpirun -n " + std::to_string(m_comm.size()) + "./CArl_FETI_setup_finish -i " +
							m_CArl_FETI_setup_finish_input_filename;

		this->print_PBS_script(	m_CArl_FETI_setup_finish_script_filename, "CArl_setup_f",
							pbs_output, pbs_error, common_script,
							command_to_run);
		
		pbs_output = m_input_params.scratch_folder_path + "/output_CArl_FETI_iterate.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_CArl_FETI_iterate.txt";
		command_to_run = "mpirun -n " + std::to_string(m_comm.size()) + "./CArl_FETI_iterate -i " +
							m_CArl_FETI_iterate_input_filename;

		this->print_PBS_script(	m_CArl_FETI_iterate_script_filename, "CArl_iter",
							pbs_output, pbs_error, common_script,
							command_to_run);

		pbs_output = m_input_params.scratch_folder_path + "/output_CArl_FETI_solution.txt";
		pbs_error  = m_input_params.scratch_folder_path + "/error_CArl_FETI_solution.txt";
		command_to_run = "mpirun -n " + std::to_string(m_comm.size()) + "./CArl_FETI_solution -i " +
							m_CArl_FETI_solution_input_filename;

		this->print_PBS_script(	m_CArl_FETI_solution_script_filename, "CArl_sol",
							pbs_output, pbs_error, common_script,
							command_to_run);
	}
}

void Solver_Files_Setup::generate_FETI_launch_scripts()
{

	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversFiles,"External solver scripts not set yet!");
	homemade_assert_msg(m_bSetCArlFETIScripts,"CArl_FETI scripts not set yet!");

	m_FETI_init_launch_script_filename = m_input_params.scratch_folder_path + "/FETI_init_script.sh";
	m_FETI_iter_launch_script_filename = m_input_params.scratch_folder_path + "/FETI_iter_script.sh";
	m_FETI_sol_launch_script_filename  = m_input_params.scratch_folder_path + "/FETI_sol_script.sh";

	switch (m_input_params.scheduler)
	{
		case ClusterSchedulerType::LOCAL :	this->generate_FETI_launch_scripts_LOCAL();
						break;

		case ClusterSchedulerType::PBS :    this->generate_FETI_launch_scripts_PBS();
						break;

		case ClusterSchedulerType::SLURM :	homemade_error_msg("Scheduler not implemented yet!");
						break;
		default : homemade_error_msg("Invalid scheduler name!");
	}

	m_bSetFETILaunchScripts = true;
};

void Solver_Files_Setup::generate_FETI_launch_scripts_LOCAL()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversFiles,"External solver scripts not set yet!");
	homemade_assert_msg(m_bSetCArlFETIScripts,"CArl_FETI scripts not set yet!");

	// ONLY WORK ON PROCESSOR 0 !!!
	if(m_comm.rank() == 0)
	{
		std::ofstream FETI_init_script(m_FETI_init_launch_script_filename);
		FETI_init_script << "#!/bin/bash" << std::endl;
		FETI_init_script << std::endl;
		FETI_init_script << ". " << m_ext_solver_u0_A_script_filename << std::endl;
		FETI_init_script << ". " << m_ext_solver_u0_B_script_filename << std::endl;
		// --- We will only need job2_i if the rigid body modes are needed
		if(m_input_params.bUseRigidBodyModes)
		{
			FETI_init_script << ". " << m_ext_solver_A_script_filename << std::endl;
			FETI_init_script << ". " << m_ext_solver_B_script_filename << std::endl;
			FETI_init_script << ". " << m_CArl_FETI_setup_finish_script_filename << std::endl;
		} else {
			FETI_init_script << ". " << m_CArl_FETI_setup_finish_script_filename << std::endl;
		}
		FETI_init_script.close();

		std::ofstream FETI_iter_script(m_FETI_iter_launch_script_filename);
		FETI_iter_script << "#!/bin/bash" << std::endl;
		FETI_iter_script << std::endl;
		FETI_iter_script << ". " << m_ext_solver_A_script_filename << std::endl;
		FETI_iter_script << ". " << m_ext_solver_B_script_filename << std::endl;
		FETI_iter_script << ". " << m_CArl_FETI_iterate_script_filename << std::endl;
		FETI_iter_script.close();

		std::ofstream FETI_set_sol_script(m_FETI_sol_launch_script_filename);
		FETI_set_sol_script << "#!/bin/bash" << std::endl;
		FETI_set_sol_script << std::endl;
		FETI_set_sol_script << ". " << m_ext_solver_A_script_filename << std::endl;
		FETI_set_sol_script << ". " << m_ext_solver_B_script_filename << std::endl;
		FETI_set_sol_script << ". " << m_CArl_FETI_solution_script_filename << std::endl;
		FETI_set_sol_script.close();
	}
};

void Solver_Files_Setup::generate_FETI_launch_scripts_PBS()
{
	homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
	homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSetExternalSolversFiles,"External solver scripts not set yet!");
	homemade_assert_msg(m_bSetCArlFETIScripts,"CArl_FETI scripts not set yet!");

	// ONLY WORK ON PROCESSOR 0 !!!
	if(m_comm.rank() == 0)
	{
		std::ofstream FETI_init_script(m_FETI_init_launch_script_filename);
		FETI_init_script << "#!/bin/bash" << std::endl;
		FETI_init_script << std::endl;
		FETI_init_script << "job1_A = `qsub " << m_ext_solver_u0_A_script_filename << "`" << std::endl;
		FETI_init_script << "job1_B = `qsub " << m_ext_solver_u0_B_script_filename << "`" << std::endl;
		// --- We will only need job2_i if the rigid body modes are needed
		if(m_input_params.bUseRigidBodyModes)
		{
			FETI_init_script << "job2_A = `qsub " << m_ext_solver_A_script_filename << "`" << std::endl;
			FETI_init_script << "job2_B = `qsub " << m_ext_solver_B_script_filename << "`" << std::endl;
			FETI_init_script << "job3 = `qsub -W depend=afterok:$job1_A:$job1_B:$job2_A:$job2_B "
							<< m_CArl_FETI_setup_finish_script_filename << "`" << std::endl;
		} else {
			FETI_init_script << "job3 = `qsub -W depend=afterok:$job1_A:$job1_B "
							<< m_CArl_FETI_setup_finish_script_filename << "`" << std::endl;
		}
		FETI_init_script.close();

		std::ofstream FETI_iter_script(m_FETI_iter_launch_script_filename);
		FETI_iter_script << "#!/bin/bash" << std::endl;
		FETI_iter_script << std::endl;
		FETI_iter_script << "job4_A = `qsub " << m_ext_solver_A_script_filename << "`" << std::endl;
		FETI_iter_script << "job4_B = `qsub " << m_ext_solver_B_script_filename << "`" << std::endl;
		FETI_iter_script << "job5 = `qsub -W depend=afterok:$job4_A:$job4_B "
							<< m_CArl_FETI_iterate_script_filename << "`" << std::endl;
		FETI_iter_script.close();

		std::ofstream FETI_set_sol_script(m_FETI_sol_launch_script_filename);
		FETI_set_sol_script << "#!/bin/bash" << std::endl;
		FETI_set_sol_script << std::endl;
		FETI_set_sol_script << "job6_A = `qsub " << m_ext_solver_A_script_filename << "`" << std::endl;
		FETI_set_sol_script << "job6_B = `qsub " << m_ext_solver_B_script_filename << "`" << std::endl;
		FETI_set_sol_script << "job7 = `qsub -W depend=afterok:$job4_A:$job4_B "
							<< m_CArl_FETI_solution_script_filename << "`" << std::endl;
		FETI_set_sol_script.close();
	}
};

}