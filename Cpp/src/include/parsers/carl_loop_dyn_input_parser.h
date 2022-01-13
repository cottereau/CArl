/*
 * carl_loop_dyn_input_parser.h 
 *
 *  Created on: Nov 17, 2021
 *      Author: Chensheng Luo
 */

#ifndef CARL_LOOP_DYN_INPUT_PARSER_H_
#define CARL_LOOP_DYN_INPUT_PARSER_H_

#include "carl_headers.h"

namespace carl
{
/// Structure containing the parameters for begining the dynamic loop
struct feti_loop_dyn_params {
	// --- Parameters used directly by the CArl_loop_dyn program

	// Cluster 
	ClusterSchedulerType scheduler;     ///< Cluster scheduler software type. *Values*: PBS, SLURM (code not implemented for the later yet).
    std::string script_filename;         ///< Cluster script file name
	
	// Path to "scratch" folder
	std::string scratch_folder_path;	///< Path to the folder which will be used to save the temporary files during the solve operation
    std::string result_folder_path;     ///< Path to the folder which will be used to save the result of all calculation.
	
    // Path of calculated matrix and vectors
    std::string tilde_matrix_A;          ///< mass tilde matrix for A
    std::string coupling_matrix_A;       ///< coupling matrix for A
    std::string stiffness_matrix_A;      ///< stiffness matrix for A
    std::string rhs_vector_A;            ///< rhs vector for A
    std::string force_folder_A;

    std::string tilde_matrix_B;          ///< mass tilde matrix for B
    std::string coupling_matrix_B;       ///< coupling matrix for B
    std::string stiffness_matrix_B;      ///< stiffness matrix for B
    std::string rhs_vector_B;            ///< rhs vector for B
    std::string force_folder_B;

    std::string interpolation_matrix;    ///< interpolation matrix

    std::string ext_solver_A_input;       ///Input parameters of the external solver for system A.
    std::string ext_solver_B_input;       ///Input parameters of the external solver for system B.
    std::string ext_solver_general_input; ///Input parameters of the external solver for general system(to solve interpolation).

    std::string ext_solver_launch_script;       ///Input parameters of the external solve for system A.
    std::string general_entry_file_path;         /// Path of entry file

    // Newmark Parameter
    double betaA; 
    double betaB; 
    double gammaA;
    double gammaB;
    double deltatA;  // in s
    double deltatB;  // in s
    double simulation_duration; // in s

    int inner_loop_times;
    int outer_loop_times;

    carl::BaseCGPrecondType CG_precond_type;	///< [CG] Type of preconditionner.
};

void get_input_params(GetPot& field_parser,
		feti_loop_dyn_params& input_params);
};
#endif /* CARL_LOOP_DYN_INPUT_PARSER_H_ */