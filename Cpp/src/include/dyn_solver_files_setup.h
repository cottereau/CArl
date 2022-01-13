/*
 * dyn_solver_files_setup.h
 *
 *  Created on: Dec 2, 2021
 *      Author: Severin Meo, Chensheng Luo
 */

#ifndef DYN_SOLVER_FILES_SETUP_H_
#define DYN_SOLVER_FILES_SETUP_H_

#include "carl_headers.h"
#include "carl_loop_dyn_input_parser.h"
#include "libmesh_solve_linear_system_input_parser.h"
#include "carl_loop_dyn_iteration_progression_parser.h"

namespace carl
{
class Dyn_Solver_Files_Setup
{
protected:
    libMesh::Parallel::Communicator& m_comm;

    feti_loop_dyn_params m_input_params;
    bool m_bInputParamsSet;
    bool m_bScratchFolderExists;


    bool m_bSetExternalSolversInputFiles ;
    std::string m_ext_solver_Afree_input_filename;
    std::string m_ext_solver_Bfree_input_filename;
    std::string m_ext_solver_couplingIV_input_filename;
    std::string m_ext_solver_Blink_input_filename;
    std::string m_ext_solver_Alink_input_filename;

    bool m_bSetExternalSolversScripts;
    std::string m_ext_solver_Afree_script_filename;
    std::string m_ext_solver_Bfree_script_filename;
    std::string m_ext_solver_couplingIV_script_filename;
    std::string m_ext_solver_Blink_script_filename;
    std::string m_ext_solver_Alink_script_filename;

    bool m_bSetInnerOperationScripts;
    std::string m_inner_operation_Afree_script_filename;
    std::string m_inner_operation_Bfree_script_filename;
    std::string m_inner_operation_setup_interpolation_script_filename;
    std::string m_inner_operation_prepare_Blink_script_filename;
    std::string m_inner_operation_Blink_script_filename;
    std::string m_inner_operation_prepare_Alink_script_filename;
    std::string m_inner_operation_Alink_script_filename;

    bool m_bSetCombinedScripts;
    std::string m_Afree_combined_filename;
    std::string m_Bfree_combined_filename;
    std::string m_coupling_combined_filename;
    std::string m_Blink_combined_filename;
    std::string m_Alink_combined_filename;

    bool m_bSetProgressionInput;
    std::string progression_input_filename;

    void generate_libmesh_external_solver_scripts_SLURM();
    void generate_combined_scripts_SLURM();
    void generate_inner_operation_scripts_SLURM();
    void print_SLURM_script(const std::string& output_filename, 
        const std::string& job_name, 
        const std::string& output_name, 
        const std::string& error_name, 
        const std::string& common_script, 
        const std::string& command_to_run);

    
  
public:
    Dyn_Solver_Files_Setup(libMesh::Parallel::Communicator& comm) :
    m_comm { comm },
    m_bInputParamsSet { false },
    m_bScratchFolderExists { false },
    m_bSetExternalSolversInputFiles { false },
    m_bSetExternalSolversScripts { false },
    m_bSetInnerOperationScripts{ false },
    m_bSetCombinedScripts{ false },
    m_bSetProgressionInput{ false }
  {
  };

    Dyn_Solver_Files_Setup(libMesh::Parallel::Communicator& comm, feti_loop_dyn_params& input_params) :
    m_comm { comm },
    m_bInputParamsSet { false },
    m_bScratchFolderExists { false },
    m_bSetExternalSolversInputFiles { false },
    m_bSetExternalSolversScripts { false },
    m_bSetInnerOperationScripts{ false },
    m_bSetCombinedScripts{ false },
    m_bSetProgressionInput{ false }
  {
    this->set_FETI_input_parameters(input_params);
  };

  void set_FETI_input_parameters(feti_loop_dyn_params& input_params);

  void set_scratch_folder();

  void generate_libmesh_external_solver_inputs();
  void generate_libmesh_external_solver_script();
  void generate_inner_operation_script();
  void generate_combined_scripts();
  void generate_progression_inputs();
  
};
}

#endif /* DYN_SOLVER_FILES_SETUP_H_ */
