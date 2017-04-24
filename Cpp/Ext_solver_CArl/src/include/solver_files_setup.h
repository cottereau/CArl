/*
 * solver_files_setup.h
 *
 *  Created on: Apr 23, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef SOLVER_FILES_SETUP_H_
#define SOLVER_FILES_SETUP_H_

#include "carl_headers.h"
#include "carl_feti_setup_init_input_parser.h"
#include "libmesh_solve_linear_system_input_parser.h"

namespace carl
{
class	Solver_Files_Setup
{
protected:

	libMesh::Parallel::Communicator& m_comm;

	feti_setup_init_params m_input_params;
	bool m_bInputParamsSet;
	bool m_bScratchFolderExists;

	bool m_bSetExternalSolversInputFiles;
	std::string m_ext_solver_u0_A_input_filename;
	std::string m_ext_solver_u0_B_input_filename;

	std::string m_ext_solver_x0_A_input_filename;
	std::string m_ext_solver_x0_B_input_filename;

	std::string m_ext_solver_yk_A_input_filename;
	std::string m_ext_solver_yk_B_input_filename;

	std::string m_ext_solver_xf_A_input_filename;
	std::string m_ext_solver_xf_B_input_filename;

	bool m_bSetExternalSolversFiles;
	std::string m_ext_solver_u0_A_script_filename;
	std::string m_ext_solver_u0_B_script_filename;

	std::string m_ext_solver_x0_A_script_filename;
	std::string m_ext_solver_x0_B_script_filename;

	std::string m_ext_solver_yk_A_script_filename;
	std::string m_ext_solver_yk_B_script_filename;

	std::string m_ext_solver_xf_A_script_filename;
	std::string m_ext_solver_xf_B_script_filename;

	bool m_bSetCArlFETIInputs;
	std::string m_CArl_FETI_init_finish_input_filename;
	std::string m_CArl_FETI_iterate_input_filename;
	std::string m_CArl_FETI_set_sol_input_filename;

	bool m_bSetCArlFETIScripts;
	std::string m_CArl_FETI_init_finish_script_filename;
	std::string m_CArl_FETI_iterate_script_filename;
	std::string m_CArl_FETI_set_sol_script_filename;

	bool m_bSetFETILaunchScripts;
	std::string m_FETI_init_launch_script_filename;
	std::string m_FETI_iter_launch_script_filename;
	std::string m_FETI_sol_launch_script_filename;

	void print_feti_setup_finish_params(const std::string& output_filename);
	void print_feti_iterate_params(const std::string& output_filename);
	void print_feti_solution_params(const std::string& output_filename);
	void print_PBS_script(const std::string& output_filename, const std::string& job_name, const std::string& output_name, const std::string& error_name, const std::string& common_script, const std::string& command_to_run);
	
	Solver_Files_Setup();
public:
	Solver_Files_Setup(libMesh::Parallel::Communicator& comm) :
		m_comm { comm },
		m_bInputParamsSet { false },
		m_bScratchFolderExists { false },
		m_bSetExternalSolversInputFiles { false },
		m_bSetExternalSolversFiles { false },
		m_bSetCArlFETIInputs { false },
		m_bSetCArlFETIScripts { false },
		m_bSetFETILaunchScripts { false }
	{
	};

	Solver_Files_Setup(libMesh::Parallel::Communicator& comm, feti_setup_init_params& input_params) :
		m_comm { comm },
		m_bInputParamsSet { false },
		m_bScratchFolderExists { false },
		m_bSetExternalSolversInputFiles { false },
		m_bSetExternalSolversFiles { false },
		m_bSetCArlFETIInputs { false },
		m_bSetCArlFETIScripts { false },
		m_bSetFETILaunchScripts { false }
	{
		this->set_FETI_input_parameters(input_params);
	}

	void set_FETI_input_parameters(feti_setup_init_params& input_params);

	void set_scratch_folder();

	void generate_libmesh_external_solver_scripts();

	void generate_libmesh_external_solver_inputs();

	void generate_FETI_scripts();

	void generate_FETI_inputs();

	void generate_FETI_launch_scripts();
};
}

#endif /* SOLVER_FILES_SETUP_H_ */