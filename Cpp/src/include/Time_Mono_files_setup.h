/*
 * Time_Mono_files_setup.h
 *
 *  Created on: July 22, 2026
 *      Author: Romain Ruyssen
 */

#ifndef TIME_MONO_FILES_SETUP_H_
#define TIME_MONO_FILES_SETUP_H_

#include "carl_headers.h"
#include "carl_time_mono_setup_input_parser.h"

namespace carl
{
class	Time_Mono_Files_Setup
{
protected:

    libMesh::Parallel::Communicator& m_comm;

	time_mono_setup_params m_input_params;        // Structure contenant tous les paramètres d'entrée pour le setup du solveur Time_Mono.
	bool m_bInputParamsSet;                      // Booléen indiquant si les paramètres d'entrée ont été définis.
	bool m_bScratchFolderExists;                 // Booléen indiquant si le dossier temporaire a été créé.
	bool m_bStateFileExists;                     // Booléen indiquant si le fichier d'état a été créé.
	bool m_bDomainStateFoldersExists;            // Booléen indiquant si les dossiers d'état des domaines ont été créés.
	bool m_bResultsFolderExists;                 // Booléen indiquant si le dossier de résultats a été créé.
	
	
	bool m_bSetExternalSolversRhsInputFiles;				 // Booléen indiquant si les paramètres d'entrée pour le solveur Time_Mono ont été définis.
	bool m_bSetExternalSolversRhsScriptFiles;
	std::string m_ext_solver_assemble_rhs_A_input_filename;
	std::string m_ext_solver_assemble_rhs_B_input_filename;
	std::string m_ext_solver_assemble_rhs_A_script_filename;
	std::string m_ext_solver_assemble_rhs_B_script_filename;

	bool m_bSetExternalSolversUpdateInitCondInputFiles;				 // Booléen indiquant si les paramètres d'entrée pour le solveur Time_Mono ont été définis.
	bool m_bSetExternalSolversUpdateInitCondScriptFiles;
	std::string m_ext_solver_update_init_cond_A_input_filename;
	std::string m_ext_solver_update_init_cond_B_input_filename;
	std::string m_ext_solver_update_init_cond_A_script_filename;
	std::string m_ext_solver_update_init_cond_B_script_filename;

	bool m_bSetFETISetupInitScriptFile;
	std::string m_FETI_setup_init_script_filename;

	bool m_bSetCArlTimeMonoInputs;				 // Booléen indiquant si les paramètres d'entrée pour le solveur Time_Mono ont été définis.
	std::string m_CArl_Time_Mono_iterate_finish_input_filename;
	bool m_bSetCArlTimeMonoScripts;				 // Booléen indiquant si les scripts pour le solveur Time_Mono ont été définis.
	std::string m_CArl_Time_Mono_iterate_finish_script_filename;
	bool m_bSetCArlTimeMonoLaunchScripts;        // Booléen indiquant si les scripts de lancement des deux état de Time_Mono ont été générés.
	std::string m_Time_Mono_iterate_init_script_filename;
	std::string m_Time_Mono_iterate_finish_script_filename;
	
	void print_time_mono_iterate_params(const std::string& output_filename); // Méthode pour générer le fichier d'inputs du script CArl_Time_Mono_setup_finish.
	void print_PBS_script(const std::string& output_filename, const std::string& job_name, const std::string& output_name, const std::string& error_name, const std::string& common_script, const std::string& command_to_run);

	Time_Mono_Files_Setup();

public:
	
	Time_Mono_Files_Setup(libMesh::Parallel::Communicator& comm, time_mono_setup_params& input_params) :
        m_comm { comm },
		m_bInputParamsSet { false },
		m_bScratchFolderExists { false },
		m_bStateFileExists { false },
		m_bDomainStateFoldersExists { false },
		m_bResultsFolderExists { false },
		m_bSetExternalSolversRhsInputFiles { false },
		m_bSetExternalSolversRhsScriptFiles { false },
		m_bSetExternalSolversUpdateInitCondInputFiles { false },
		m_bSetExternalSolversUpdateInitCondScriptFiles { false },
		m_bSetFETISetupInitScriptFile { false },
		m_bSetCArlTimeMonoInputs { false },
		m_bSetCArlTimeMonoScripts { false },
		m_bSetCArlTimeMonoLaunchScripts { false }
	{
		this->set_Time_Mono_input_parameters(input_params);
	}

	void set_Time_Mono_input_parameters(time_mono_setup_params& input_params);
	void set_scratch_folder();
	void generate_state_file();
	void generate_domain_state_folders();
	void set_results_folder();
	void copy_initial_conditions_files();

	void generate_libmesh_external_solver_assemble_rhs_inputs();
	void generate_libmesh_external_solver_assemble_rhs_scripts();
	void generate_libmesh_external_solver_assemble_rhs_scripts_LOCAL();
	void generate_libmesh_external_solver_assemble_rhs_scripts_PBS();

	void generate_libmesh_external_solver_update_init_cond_inputs();
	void generate_libmesh_external_solver_update_init_cond_scripts();
	void generate_libmesh_external_solver_update_init_cond_scripts_LOCAL();
	void generate_libmesh_external_solver_update_init_cond_scripts_PBS();

	// Pour le moment, l'input est directement prérempli dans le dossier scratch 
	void generate_FETI_setup_init_scripts(); 
	void generate_FETI_setup_init_scripts_LOCAL();
	void generate_FETI_setup_init_scripts_PBS();

	void generate_Time_Mono_inputs();
	void generate_Time_Mono_scripts();
	void generate_Time_Mono_scripts_LOCAL();
	void generate_Time_Mono_scripts_PBS();
	void generate_Time_Mono_launch_scripts();
	void generate_Time_Mono_launch_scripts_LOCAL();
	void generate_Time_Mono_launch_scripts_PBS();

};
}

#endif /* SOLVER_FILES_SETUP_H_ */