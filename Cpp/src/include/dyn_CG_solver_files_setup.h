/*
 * dyn_CG_solver_files_setup.h
 *
 *  Created on: April 26, 2022
 *      Author: Chensheng Luo
 */

#ifndef DYN_CG_SOLVER_FILES_SETUP_H_
#define DYN_CG_SOLVER_FILES_SETUP_H_

#include "carl_headers.h"
#include "carl_loop_dyn_input_parser.h"
#include "libmesh_solve_linear_system_input_parser.h"
#include "carl_loop_dyn_iteration_progression_parser.h"

namespace carl
{
/** \brief **DYN-CG** Class generating all necessary `.sh` script and input for CArl dynamic time loop.
 * 
 *  This class is used in CArl_loop_dyn_setup.cpp to generate all the external solver lancing script and 
 *  input files, thay is to say:
 *      - Scratch folder, result folder
 *      - External solver launching script, *i.e.* `ext_solver_***.sh`
 *      - External solver input files, *i.e.* `ext_solver_***.txt`
 *      - `CArl_loop_dyn_***` launching script, *i.e.* `inner_ope_***.sh`
 *      - 6 general steps launching script, *i.e.* `Afree.sh` and so on
 *      - Progression recording file, *i.e.* `iteration_progression.txt`
 *
 *
 */
class Dyn_CG_Solver_Files_Setup
{
protected:
    libMesh::Parallel::Communicator& m_comm;

    feti_loop_dyn_params m_input_params;
    bool m_bInputParamsSet;
    bool m_bScratchFolderExists;

    bool m_bSetExternalSolversInputFiles ;
    std::string m_ext_solver_Afree_input_filename;
    std::string m_ext_solver_Bfree_input_filename;
    std::string m_ext_solver_CG_A_input_filename;
    std::string m_ext_solver_CG_B_input_filename;

    bool m_bSetExternalSolversScripts;
    std::string m_ext_solver_Afree_script_filename;
    std::string m_ext_solver_Bfree_script_filename;
    std::string m_ext_solver_CG_A_script_filename;
    std::string m_ext_solver_CG_B_script_filename;

    bool m_bSetInnerOperationScripts;
    std::string m_inner_operation_Afree_script_filename;
    std::string m_inner_operation_Bfree_script_filename;
    std::string m_inner_operation_coupling_init_script_filename;
    std::string m_inner_operation_coupling_finish_script_filename;
    std::string m_inner_operation_coupling_iterate_script_filename;
    std::string m_inner_operation_coupling_solution_script_filename;
    std::string m_inner_operation_Blink_script_filename;
    std::string m_inner_operation_Alink_script_filename;

    bool m_bSetCombinedScripts;
    std::string m_Afree_combined_filename;
    std::string m_Bfree_combined_filename;
    std::string m_coupling_combined_filename;
    std::string m_coupling_init_combined_filename;
    std::string m_coupling_iterate_combined_filename;
    std::string m_coupling_solution_combined_filename;

    bool m_bSetProgressionInput;
    std::string progression_input_filename;

    void generate_libmesh_external_solver_scripts_SLURM();   ///< Generate SLURM script to launch external solver
    void generate_libmesh_external_solver_scripts_LOCAL();   ///< Generate SLURM script to launch external solver
    void generate_libmesh_external_solver_scripts_PBS();   ///< Generate SLURM script to launch external solver
    void generate_combined_scripts_SLURM();    ///< Generate `Afree.sh` and so on for SLRUM
    void generate_inner_operation_scripts_SLURM();   ///< Generate `inner_ope_Afree.sh` and so on for SLRUM
    void print_SLURM_script(const std::string& output_filename, 
        const std::string& job_name, 
        const std::string& output_name, 
        const std::string& error_name, 
        const std::string& common_script, 
        const std::string& command_to_run);   ///< Print SLURM script into a file to launch a job.
    void print_LOCAL_script(const std::string& output_filename,
 const std::string& job_name,
 const std::string& output_name,
 const std::string& error_name,
 const std::string& common_script,
 const std::string& command_to_run);
    void print_PBS_script(const std::string& output_filename,
 const std::string& job_name,
 const std::string& output_name,
 const std::string& error_name,
 const std::string& common_script,
 const std::string& command_to_run);
    
  
public:
    Dyn_CG_Solver_Files_Setup(libMesh::Parallel::Communicator& comm) :
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

    Dyn_CG_Solver_Files_Setup(libMesh::Parallel::Communicator& comm, feti_loop_dyn_params& input_params) :
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

  void set_scratch_folder();   ///< Make scratch folder, result folder and force folder

  void generate_libmesh_external_solver_inputs();   ///< Generate external solver input file `ext_solver_***.txt`
  void generate_libmesh_external_solver_script();   ///< Generate external solver script `ext_solver_***.sh`
  void generate_inner_operation_script();   ///< Generate `CArl_loop_dyn_***` launching script `inner_ope_***.sh`
  void generate_combined_scripts();    ///< Generate 6 general steps launching script, *i.e.* `Afree.sh` and so on. for SLRUM
  void generate_progression_inputs();   ///< Generate `iteration_progression.txt` used to record progression
};
}

#endif /* DYN_CG_SOLVER_FILES_SETUP_H_ */
