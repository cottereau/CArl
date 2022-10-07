/*
 * \file dyn_CG_solver_files_setup.cpp
 *
 *  Created on: April 26, 2022
 *      Author: Chensheng Luo
 * 
 * \brief **DYN-CG**   functions responsible for generating all scripts of Solving step of CG solver
 */

#include "dyn_CG_solver_files_setup.h"

namespace carl
{

void Dyn_CG_Solver_Files_Setup::print_LOCAL_script(const std::string& output_filename,
 const std::string& job_name,
 const std::string& output_name,
 const std::string& error_name,
 const std::string& common_script,
 const std::string& command_to_run)
{
  //In this file, it's not necessary to 
  std::ofstream output_script(output_filename);
  output_script << command_to_run << std::endl;
  output_script.close();
};

void Dyn_CG_Solver_Files_Setup::print_PBS_script(const std::string& output_filename,
 const std::string& job_name,
 const std::string& output_name,
 const std::string& error_name,
 const std::string& common_script,
 const std::string& command_to_run)
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

void Dyn_CG_Solver_Files_Setup::print_SLURM_script(const std::string& output_filename,
 const std::string& job_name,
 const std::string& output_name,
 const std::string& error_name,
 const std::string& common_script,
 const std::string& command_to_run)
{
  std::ofstream output_script(output_filename);
  output_script << "#!/bin/bash" << std::endl;
  output_script << std::endl;
  output_script << "#SBATCH --job-name=" << job_name << std::endl;
  output_script << "#SBATCH --output=" << output_name << std::endl;
  output_script << "#SBATCH --error=" << error_name << std::endl;
  output_script << common_script << std::endl;
  output_script << command_to_run << std::endl;
  output_script.close();
};

void Dyn_CG_Solver_Files_Setup::set_FETI_input_parameters(feti_loop_dyn_params& input_params)
{
  m_bInputParamsSet = true;
  m_input_params = input_params;
};

void Dyn_CG_Solver_Files_Setup::set_scratch_folder()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");

  m_bScratchFolderExists = true;
  
  if(m_comm.rank() == 0)
  {
    std::string command_string;

    command_string = "rm -rf " + m_input_params.result_folder_path;
    carl::exec_command(command_string.c_str());

    command_string = "mkdir -p " + m_input_params.result_folder_path;
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "rm -rf " + m_input_params.scratch_folder_path;
    carl::exec_command(command_string.c_str());


    command_string = "mkdir -p " + m_input_params.scratch_folder_path;
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "mkdir -p " + m_input_params.scratch_folder_path + "/CG_solver";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;
  }
};

void Dyn_CG_Solver_Files_Setup::generate_libmesh_external_solver_inputs()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

  m_bSetExternalSolversInputFiles = true;
  
  m_ext_solver_Afree_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Afree.txt";
  m_ext_solver_Bfree_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Bfree.txt"; 
  m_ext_solver_CG_A_input_filename = m_input_params.scratch_folder_path + "/CG_solver/ext_solver_CG_A.txt"; 
  m_ext_solver_CG_B_input_filename = m_input_params.scratch_folder_path + "/CG_solver/ext_solver_CG_B.txt"; 


// Get general input parameters
  GetPot field_parser_A;
  field_parser_A.parse_input_file(m_input_params.ext_solver_A_input, "#", "\n", " \t\n");

  GetPot field_parser_B;                                                                   
  field_parser_B.parse_input_file(m_input_params.ext_solver_B_input, "#", "\n", " \t\n");  

  GetPot field_parser_CG_A;                                                                        
  field_parser_CG_A.parse_input_file(m_input_params.ext_solver_A_input, "#", "\n", " \t\n"); 

  GetPot field_parser_CG_B;                                                                        
  field_parser_CG_B.parse_input_file(m_input_params.ext_solver_B_input, "#", "\n", " \t\n"); 

// Get general input parameters
  carl::libmesh_solve_linear_system_input_params solver_Afree_input_params;
  carl::get_input_params(field_parser_A, solver_Afree_input_params);

  carl::libmesh_solve_linear_system_input_params solver_Bfree_input_params;   
  carl::get_input_params(field_parser_B, solver_Bfree_input_params);          

  carl::libmesh_solve_linear_system_input_params solver_CG_A_input_params;   
  carl::get_input_params(field_parser_CG_A, solver_CG_A_input_params); 

  carl::libmesh_solve_linear_system_input_params solver_CG_B_input_params;   
  carl::get_input_params(field_parser_CG_B, solver_CG_B_input_params);            

  // Set input parameters
  solver_Afree_input_params.sys_matrix_file = m_input_params.matrix_A.mass_tilde;
  solver_Afree_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_A_free.petscvec";
  solver_Afree_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_A_free";

  solver_Bfree_input_params.sys_matrix_file = m_input_params.matrix_B.mass_tilde;                                      
  solver_Bfree_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_B_free.petscvec"; 
  solver_Bfree_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_B_free";                 

  solver_CG_A_input_params.sys_matrix_file = m_input_params.scratch_folder_path + "/CG_solver/Mtilde_A.petscmat"; 
  solver_CG_A_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/CG_solver/ext_solver_A_rhs.petscvec"; 
  solver_CG_A_input_params.output_base = m_input_params.scratch_folder_path + "/CG_solver/ext_solver_A";          
                          
  solver_CG_B_input_params.sys_matrix_file = m_input_params.scratch_folder_path + "/CG_solver/Mtilde_B.petscmat"; 
  solver_CG_B_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/CG_solver/ext_solver_B_rhs.petscvec"; 
  solver_CG_B_input_params.output_base = m_input_params.scratch_folder_path + "/CG_solver/ext_solver_B";          

  carl::print_input_params(m_ext_solver_Afree_input_filename,solver_Afree_input_params);
  carl::print_input_params(m_ext_solver_Bfree_input_filename, solver_Bfree_input_params); 
  carl::print_input_params(m_ext_solver_CG_A_input_filename, solver_CG_A_input_params); 
  carl::print_input_params(m_ext_solver_CG_B_input_filename, solver_CG_B_input_params); 

  
}

void Dyn_CG_Solver_Files_Setup::generate_libmesh_external_solver_script()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");

  m_ext_solver_Afree_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Afree_acc.sh";
  m_ext_solver_Bfree_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Bfree_acc.sh"; 
  m_ext_solver_CG_A_script_filename = m_input_params.scratch_folder_path + "/CG_solver/ext_solver_CG_A.sh"; 
  m_ext_solver_CG_B_script_filename = m_input_params.scratch_folder_path + "/CG_solver/ext_solver_CG_B.sh"; 


  switch (m_input_params.scheduler)
  {
    case ClusterSchedulerType::LOCAL :  this->generate_libmesh_external_solver_scripts_LOCAL();
            break;

    case ClusterSchedulerType::PBS :    this->generate_libmesh_external_solver_scripts_PBS();
            break;

    case ClusterSchedulerType::SLURM :  this->generate_libmesh_external_solver_scripts_SLURM();
            break;
    default : homemade_error_msg("Invalid scheduler name!");
  }

  m_bSetExternalSolversScripts = true;
}

void Dyn_CG_Solver_Files_Setup::generate_libmesh_external_solver_scripts_LOCAL()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");

  if(m_comm.rank() == 0)
  {
    std::string command_to_run;

    // Set the Afree_acc scripts
    command_to_run = m_input_params.ext_solver_launch_script_A + " " + m_ext_solver_Afree_input_filename;
    std::ofstream output_script(m_ext_solver_Afree_script_filename);
    output_script << command_to_run << std::endl;
    output_script.close();

    // Set the Bfree_acc scripts
    command_to_run = m_input_params.ext_solver_launch_script_B + " " + m_ext_solver_Bfree_input_filename;     
    std::ofstream output_scriptB(m_ext_solver_Bfree_script_filename);
    output_scriptB << command_to_run << std::endl;
    output_scriptB.close();


    // Set the CG_A scripts                               
    command_to_run = m_input_params.ext_solver_launch_script_A + " " + m_ext_solver_CG_A_input_filename; 
    std::ofstream output_scriptC(m_ext_solver_CG_A_script_filename);
    output_scriptC << command_to_run << std::endl;
    output_scriptC.close();

    // Set the CG_B scripts
    command_to_run = m_input_params.ext_solver_launch_script_B + " " + m_ext_solver_CG_B_input_filename;     
    std::ofstream output_scriptD(m_ext_solver_CG_B_script_filename);
    output_scriptD << command_to_run << std::endl;
    output_scriptD.close();
  }
}

void Dyn_CG_Solver_Files_Setup::generate_libmesh_external_solver_scripts_PBS()
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

    // Set the Afree_acc scripts
    pbs_output = m_input_params.scratch_folder_path + "/output_Afree_acc.txt";
    pbs_error  = m_input_params.scratch_folder_path + "/error_Afree_acc.txt";
    command_to_run = m_input_params.ext_solver_launch_script_A + " " + m_ext_solver_Afree_input_filename;

    this->print_PBS_script(m_ext_solver_Afree_script_filename, "Afree_acc",
              pbs_output, pbs_error, common_script,
              command_to_run);

    // Set the Bfree_acc scripts
    pbs_output = m_input_params.scratch_folder_path + "/output_Bfree_acc.txt";                               
    pbs_error = m_input_params.scratch_folder_path + "/error_Bfree_acc.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_B + " " + m_ext_solver_Bfree_input_filename;     

    this->print_PBS_script(m_ext_solver_Bfree_script_filename, "Bfree_acc", 
              pbs_output, pbs_error, common_script,
              command_to_run);

    // Set the CG_A scripts
    pbs_output = m_input_params.scratch_folder_path + "/CG_solver/output_CG_A.txt";                            
    pbs_error = m_input_params.scratch_folder_path + "/CG_solver/error_CG_A.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_A + " " + m_ext_solver_CG_A_input_filename; 


    this->print_PBS_script(m_ext_solver_CG_A_script_filename, "CG_A",
             pbs_output, pbs_error, common_script,
             command_to_run);   

    // Set the CG_B scripts
    pbs_output = m_input_params.scratch_folder_path + "/CG_solver/output_CG_B.txt";                               
    pbs_error = m_input_params.scratch_folder_path + "/CG_solver/error_CG_B.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_B + " " + m_ext_solver_CG_B_input_filename;     

    this->print_PBS_script(m_ext_solver_CG_B_script_filename, "CG_B", 
            pbs_output, pbs_error, common_script,
            command_to_run); 
    }
  }

void Dyn_CG_Solver_Files_Setup::generate_libmesh_external_solver_scripts_SLURM(){
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

    std::string slurm_output;
    std::string slurm_error;
    std::string command_to_run;

    // Set the Afree_acc scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Afree_acc.txt";
    slurm_error  = m_input_params.scratch_folder_path + "/error_Afree_acc.txt";
    command_to_run = m_input_params.ext_solver_launch_script_A + " " + m_ext_solver_Afree_input_filename;

    this->print_SLURM_script(m_ext_solver_Afree_script_filename, "Afree_acc",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the Bfree_acc scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Bfree_acc.txt";                               
    slurm_error = m_input_params.scratch_folder_path + "/error_Bfree_acc.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_B + " " + m_ext_solver_Bfree_input_filename;     

    this->print_SLURM_script(m_ext_solver_Bfree_script_filename, "Bfree_acc", 
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the CG_A scripts
    slurm_output = m_input_params.scratch_folder_path + "/CG_solver/output_CG_A.txt";                            
    slurm_error = m_input_params.scratch_folder_path + "/CG_solver/error_CG_A.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_A + " " + m_ext_solver_CG_A_input_filename; 


    this->print_SLURM_script(m_ext_solver_CG_A_script_filename, "CG_A",
             slurm_output, slurm_error, common_script,
             command_to_run);   

    // Set the CG_B scripts
    slurm_output = m_input_params.scratch_folder_path + "/CG_solver/output_CG_B.txt";                               
    slurm_error = m_input_params.scratch_folder_path + "/CG_solver/error_CG_B.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_B + " " + m_ext_solver_CG_B_input_filename;     

    this->print_SLURM_script(m_ext_solver_CG_B_script_filename, "CG_B", 
            slurm_output, slurm_error, common_script,
            command_to_run); 
  }
}

void Dyn_CG_Solver_Files_Setup::generate_inner_operation_script()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");


  m_inner_operation_Afree_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Afree.slurm";
  m_inner_operation_Bfree_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Bfree.slurm";           
  m_inner_operation_coupling_init_script_filename = m_input_params.scratch_folder_path + "/CG_solver/inner_coupling_init.slurm"; 
  m_inner_operation_coupling_finish_script_filename = m_input_params.scratch_folder_path + "/CG_solver/inner_coupling_finish.slurm"; 
  m_inner_operation_coupling_iterate_script_filename = m_input_params.scratch_folder_path + "/CG_solver/inner_coupling_iterate.slurm";
  m_inner_operation_coupling_solution_script_filename = m_input_params.scratch_folder_path + "/CG_solver/inner_coupling_solution.slurm";
  m_inner_operation_Blink_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Blink.slurm";           
  m_inner_operation_Alink_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Alink.slurm";       



  switch (m_input_params.scheduler)
  {
    case ClusterSchedulerType::LOCAL :  homemade_error_msg("Scheduler not implemented yet!");
            break;

    case ClusterSchedulerType::PBS :    homemade_error_msg("Scheduler not implemented yet!");
            break;

    case ClusterSchedulerType::SLURM :  this->generate_inner_operation_scripts_SLURM();
            break;
    default : homemade_error_msg("Invalid scheduler name!");
  }

  m_bSetInnerOperationScripts = true;
}

void Dyn_CG_Solver_Files_Setup::generate_inner_operation_scripts_SLURM(){
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");

  if(m_comm.rank() == 0)
  {
    // Get the full common script file into a string
    std::ifstream base_script(m_input_params.script_filename);
    std::string common_script((std::istreambuf_iterator<char>(base_script)),
                  std::istreambuf_iterator<char>());
    base_script.close();

    std::string slurm_output;
    std::string slurm_error;
    std::string command_to_run;

    // Set the Afree scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Afree_speed.txt";
    slurm_error  = m_input_params.scratch_folder_path + "/error_Afree_speed.txt";
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_Afree -i "+ m_input_params.general_entry_file_path;

    this->print_SLURM_script(m_inner_operation_Afree_script_filename, "Afree_speed",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the Bfree scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Bfree_speed.txt";                            
    slurm_error  = m_input_params.scratch_folder_path + "/error_Bfree_speed.txt";                             
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_Bfree -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_Bfree_script_filename, "Bfree_speed",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                

    // Set the coupling init scripts
    slurm_output = m_input_params.scratch_folder_path + "/CG_solver/output_coupling_init.txt";                             
    slurm_error = m_input_params.scratch_folder_path + "/CG_solver/error_coupling_init.txt";                               
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_coupling_setup_init -i "+ m_input_params.general_entry_file_path+" -sub_pc_factor_shift_type positive_definite";  

    this->print_SLURM_script(m_inner_operation_coupling_init_script_filename, "CG_init",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the coupling finish scripts
    slurm_output = m_input_params.scratch_folder_path + "/CG_solver/output_coupling_finish.txt";                             
    slurm_error = m_input_params.scratch_folder_path + "/CG_solver/error_coupling_finish.txt";                               
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_coupling_setup_finish -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_coupling_finish_script_filename, "CG_finish",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the coupling iterate scripts
    slurm_output = m_input_params.scratch_folder_path + "/CG_solver/output_coupling_iterate.txt";                             
    slurm_error = m_input_params.scratch_folder_path + "/CG_solver/error_coupling_iterate.txt";                               
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_coupling_iterate -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_coupling_iterate_script_filename, "CG_iterate",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the coupling solution scripts
    slurm_output = m_input_params.scratch_folder_path + "/CG_solver/output_coupling_solution.txt";                             
    slurm_error = m_input_params.scratch_folder_path + "/CG_solver/error_coupling_solution.txt";                               
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_coupling_solution -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_coupling_solution_script_filename, "CG_solution",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the Blink scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Blink_speed.txt";                            
    slurm_error  = m_input_params.scratch_folder_path + "/error_Blink_speed.txt";                             
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_Blink -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_Blink_script_filename, "Blink_speed",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                
                                                 

    // Set the Alink scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Alink_speed.txt";                            
    slurm_error  = m_input_params.scratch_folder_path + "/error_Alink_speed.txt";                             
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_Alink -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_Alink_script_filename, "Alink_speed",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                

  }
}

void Dyn_CG_Solver_Files_Setup::generate_combined_scripts(){
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");
  homemade_assert_msg(m_bSetInnerOperationScripts,"Inner operation script files not set yet!");

  m_Afree_combined_filename = m_input_params.scratch_folder_path + "/Afree.sh";
  m_Bfree_combined_filename = m_input_params.scratch_folder_path + "/Bfree.sh";
  m_coupling_combined_filename = m_input_params.scratch_folder_path + "/coupling.sh";                     
  m_coupling_init_combined_filename = m_input_params.scratch_folder_path + "/CG_solver/coupling_init.sh";
  m_coupling_iterate_combined_filename = m_input_params.scratch_folder_path + "/CG_solver/coupling_iterate.sh";                  
  m_coupling_solution_combined_filename = m_input_params.scratch_folder_path + "/CG_solver/coupling_solution.sh";

  switch (m_input_params.scheduler)
  {
    case ClusterSchedulerType::LOCAL :  homemade_error_msg("Scheduler not implemented yet!");
            break;

    case ClusterSchedulerType::PBS :    homemade_error_msg("Scheduler not implemented yet!");
            break;

    case ClusterSchedulerType::SLURM :  this->generate_combined_scripts_SLURM();
            break;
    default : homemade_error_msg("Invalid scheduler name!");
    }

  m_bSetCombinedScripts = true;

}

void Dyn_CG_Solver_Files_Setup::generate_combined_scripts_SLURM(){
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");
  homemade_assert_msg(m_bSetInnerOperationScripts,"Inner operation script files not set yet!");

  if(m_comm.rank() == 0)
  {
    // Afree_script
    std::ofstream Afree_script(m_Afree_combined_filename);
    Afree_script << "#!/bin/bash" << std::endl;
    Afree_script << std::endl;
    Afree_script << "job_acc=$(sbatch " << m_ext_solver_Afree_script_filename << ")" << std::endl;
    Afree_script << "job_reste=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name Afree_acc) "<< m_inner_operation_Afree_script_filename << ")" << std::endl;
    Afree_script.close();

    // Bfree_script
    std::ofstream Bfree_script(m_Bfree_combined_filename);
    Bfree_script << "#!/bin/bash" << std::endl;
    Bfree_script << std::endl;
    Bfree_script << "job_acc=$(sbatch " << m_ext_solver_Bfree_script_filename << ")" << std::endl;
    Bfree_script << "job_reste=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name Bfree_acc) "<< m_inner_operation_Bfree_script_filename << ")" << std::endl;
    Bfree_script.close();

    //coupling_script
    std::ofstream coupling_script(m_coupling_combined_filename);
    coupling_script << "#!/bin/bash" << std::endl;
    coupling_script << std::endl;
    coupling_script << "job=$(sbatch " << m_inner_operation_coupling_init_script_filename << ")" << std::endl;
    coupling_script.close();

    //coupling_init_script
    std::ofstream coupling_init_script(m_coupling_init_combined_filename);
    coupling_init_script << "#!/bin/bash" << std::endl;
    coupling_init_script << std::endl;
    coupling_init_script << "job_CGA=$(sbatch " << m_ext_solver_CG_A_script_filename << ")" << std::endl;
    coupling_init_script << "job_CGB=$(sbatch " << m_ext_solver_CG_B_script_filename << ")" << std::endl;
    coupling_init_script << "job_inter=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name CG_A):$(squeue --noheader --format %i --name CG_B) " << m_inner_operation_coupling_finish_script_filename << ")" << std::endl;
    coupling_init_script.close();  

    //coupling_iterate_script
    std::ofstream coupling_iterate_script(m_coupling_iterate_combined_filename);
    coupling_iterate_script << "#!/bin/bash" << std::endl;
    coupling_iterate_script << std::endl;
    coupling_iterate_script << "job_CGA=$(sbatch " << m_ext_solver_CG_A_script_filename << ")" << std::endl;
    coupling_iterate_script << "job_CGB=$(sbatch " << m_ext_solver_CG_B_script_filename << ")" << std::endl;
    coupling_iterate_script << "job_inter=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name CG_A):$(squeue --noheader --format %i --name CG_B) " << m_inner_operation_coupling_iterate_script_filename << ")" << std::endl;
    coupling_iterate_script.close();  

    //coupling_solution_script
    std::ofstream coupling_solution_script(m_coupling_solution_combined_filename);
    coupling_solution_script << "#!/bin/bash" << std::endl;
    coupling_solution_script << std::endl;
    coupling_solution_script << "job_CGA=$(sbatch " << m_ext_solver_CG_A_script_filename << ")" << std::endl;
    coupling_solution_script << "job_CGB=$(sbatch " << m_ext_solver_CG_B_script_filename << ")" << std::endl;
    coupling_solution_script << "job_inter=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name CG_A):$(squeue --noheader --format %i --name CG_B) " << m_inner_operation_coupling_solution_script_filename << ")" << std::endl;
    coupling_solution_script.close();  
  }
}


void Dyn_CG_Solver_Files_Setup::generate_progression_inputs(){
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");
  homemade_assert_msg(m_bSetInnerOperationScripts,"Inner operation script files not set yet!");
  homemade_assert_msg(m_bSetCombinedScripts,"Combined scrpit files not set yet!");

  if(m_comm.rank() == 0)
  {
  progression_input_filename = m_input_params.scratch_folder_path + "/iteration_progression.txt";

  carl::feti_loop_dyn_iteration_progression_params progression_params;
  progression_params.inner_loop_progression = 0;
  progression_params.outer_loop_progression = 0;

  print_input_params(progression_input_filename,progression_params);
  }
  m_bSetProgressionInput=true;

}





}
