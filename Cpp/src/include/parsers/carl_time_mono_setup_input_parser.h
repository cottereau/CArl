/*
 * carl_time_mono_setup_input_parser.h
 *
 *  Created on: July 22, 2026
 *      Author: Romain Ruyssen
 */

#ifndef CARL_TIME_MONO_SETUP_INPUT_PARSER_H_
#define CARL_TIME_MONO_SETUP_INPUT_PARSER_H_

#include "carl_headers.h"

namespace carl
{
/// Structure containing the parameters for the setup initialization of the Time Mono solver.
struct time_mono_setup_params {

	// Cluster 
	ClusterSchedulerType scheduler; ///< Cluster scheduler software type. *Values*: PBS, SLURM (code not implemented for the later yet).
	
	// Path to "scratch" folder
	std::string scratch_folder_path;	///< Path to the folder which will be used to save the temporary files during the solve operation
	std::string script_filename;        ///< Path to the file used to generate the scripts.
	
	// External solver commands
	std::string ext_solver_BIG_rhs_assembly;	    ///< Command used for the external solver rhs assembly for system A.
	std::string ext_solver_micro_rhs_assembly;		///< Command used for the external solver rhs assembly for system B.
	std::string ext_solver_BIG_update_init_cond;	///< Command used for the external solver update initial conditions for system A.
	std::string ext_solver_micro_update_init_cond;	///< Command used for the external solver update initial conditions for system B.

	// Path to the time mesh file
	std::string time_mesh_file;	        ///< Path to the file containing the time mesh. 

	// Path to the newmark parameters file
	std::string newmark_parameters_A_file;      ///< Path to the file containing the newmark parameters for domain A.
	std::string newmark_parameters_B_file;      ///< Path to the file containing the newmark parameters for domain B.

	// Pathes to the space mesh files
	std::string space_mesh_A_file;		///< Path to the file containing the space mesh for domain A.
	std::string space_mesh_B_file;		///< Path to the file containing the space mesh for domain B.

	// Domains matrices files
	std::string domain_A_matrix_file;	///< Path to the file containing the matrix for domain A.
	std::string domain_A_M_matrix_file;	///< Path to the file containing the matrix for domain A.
	std::string domain_A_K_matrix_file;	///< Path to the file containing the matrix for domain A.
	std::string domain_A_C_matrix_file;	///< Path to the file containing the matrix for domain A.

	std::string domain_B_matrix_file;	///< Path to the file containing the matrix for domain B.
	std::string domain_B_M_matrix_file;	///< Path to the file containing the matrix for domain B.
	std::string domain_B_K_matrix_file;	///< Path to the file containing the matrix for domain B.
	std::string domain_B_C_matrix_file;	///< Path to the file containing the matrix for domain B.

	// Domains initial conditions files
	std::string domain_A_init_cond_disp_file;	///< Path to the file containing the initial displacement for domain A.
	std::string domain_A_init_cond_vel_file;	///< Path to the file containing the initial velocity for domain A.
	std::string domain_A_init_cond_acc_file;	///< Path to the file containing the initial acceleration for domain A.
	
	std::string domain_B_init_cond_disp_file;	///< Path to the file containing the initial displacement for domain B.
	std::string domain_B_init_cond_vel_file;	///< Path to the file containing the initial velocity for domain B.
	std::string domain_B_init_cond_acc_file;	///< Path to the file containing the initial acceleration for domain B.

	// Path to the FETI setup parameters file
	std::string feti_setup_params_file;		///< Path to the file containing the F

	// Path to the FETI solution folder
	std::string feti_solution_path;			///< FETI solution path

	// Path to the "results" folder
	std::string results_folder_path;	///< Path to the folder where the results will be saved.

	// Name of the results files
	std::string results_file_name;		///< Name of the results files. The files will be saved in the `results_folder_path` folder.

};

/**	\brief Parser function for the Time_Mono solver.

	Required parameters:
    - `ScratchFolderPath` : path to the folder where the temporary files used by the Time_Mono solver will be saved.
    - `TimeMeshFile` : path to the file containing the time mesh.
	- 'ResultsFolderPath' : path to the folder where the results will be saved.
	- `ResultsFileName` : name of the results files. The files will be saved in the `results_folder_path` folder.
	- `ResultsFileType` : type of the results files.
 */
void get_input_params(GetPot& field_parser,
		time_mono_setup_params& input_params);

}
#endif /* CARL_TIME_MONO_SETUP_INPUT_PARSER_H_ */
