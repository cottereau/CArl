/*
 * \file dyn_DI_solver_files_setup.cpp
 *
 *  Created on: Dec 2, 2021
 *      Author: Chensheng Luo, Severin Meo
 * 
 * \brief **DYN-DI**   functions responsible for generating all scripts of Solving step of CG solver
 */
#include "dyn_DI_solver_files_setup.h"



namespace carl
{

void Dyn_DI_Solver_Files_Setup::print_SLURM_script(const std::string& output_filename, const std::string& job_name, const std::string& output_name, const std::string& error_name, const std::string& common_script, const std::string& command_to_run)
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

void Dyn_DI_Solver_Files_Setup::set_FETI_input_parameters(feti_loop_dyn_params& input_params)
{
  m_bInputParamsSet = true;
  m_input_params = input_params;
};

void Dyn_DI_Solver_Files_Setup::set_scratch_folder()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");

  m_bScratchFolderExists = true;
  
  std::string command_string;

  if(m_comm.rank() == 0)
  {

    command_string = "rm -rf " + m_input_params.result_folder_path;
    std::cout << command_string << std::endl;
    carl::exec_command(command_string.c_str());

    command_string = "mkdir -p " + m_input_params.result_folder_path;
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "rm -rf " + m_input_params.scratch_folder_path;
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "mkdir -p " + m_input_params.scratch_folder_path;
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;
  }
};

void Dyn_DI_Solver_Files_Setup::generate_libmesh_external_solver_inputs()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

  m_bSetExternalSolversInputFiles = true;

  m_ext_solver_Afree_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Afree.txt";
  m_ext_solver_Bfree_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Bfree.txt"; 
  m_ext_solver_coupling_input_filename = m_input_params.scratch_folder_path + "/ext_solver_coupling.txt"; 
  m_ext_solver_Blink_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Blink.txt"; 
  m_ext_solver_Alink_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Alink.txt"; 

// Get general input parameters
  GetPot field_parser_A;
  field_parser_A.parse_input_file(m_input_params.ext_solver_A_input, "#", "\n", " \t\n");

  GetPot field_parser_B;                                                                   
  field_parser_B.parse_input_file(m_input_params.ext_solver_B_input, "#", "\n", " \t\n");  

  GetPot field_parser_coupling;                                                                        
  field_parser_coupling.parse_input_file(m_input_params.ext_solver_general_input, "#", "\n", " \t\n"); 

// Get general input parameters
  carl::libmesh_solve_linear_system_input_params solver_Afree_input_params;
  carl::get_input_params(field_parser_A, solver_Afree_input_params);

  carl::libmesh_solve_linear_system_input_params solver_Bfree_input_params;   
  carl::get_input_params(field_parser_B, solver_Bfree_input_params);          

  carl::libmesh_solve_linear_system_input_params solver_coupling_input_params;   
  carl::get_input_params(field_parser_coupling, solver_coupling_input_params); 

  carl::libmesh_solve_linear_system_input_params solver_Blink_input_params;   
  carl::get_input_params(field_parser_B, solver_Blink_input_params); 

  carl::libmesh_solve_linear_system_input_params solver_Alink_input_params; 
  carl::get_input_params(field_parser_A, solver_Alink_input_params);    

  // Set input parameters
  solver_Afree_input_params.sys_matrix_file = m_input_params.matrix_A.mass_tilde;
  solver_Afree_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_A_free.petscvec";
  solver_Afree_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_A_free";

  solver_Bfree_input_params.sys_matrix_file = m_input_params.matrix_B.mass_tilde;                                      
  solver_Bfree_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_B_free.petscvec"; 
  solver_Bfree_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_B_free";                 


  solver_coupling_input_params.sys_matrix_file = m_input_params.interpolation_matrix; 
  solver_coupling_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_interpolation_vec.petscvec"; 
  solver_coupling_input_params.output_base = m_input_params.scratch_folder_path + "/coupling_interpolation_vec";          
                          
  solver_Blink_input_params.sys_matrix_file = m_input_params.matrix_B.mass_tilde; 
  solver_Blink_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_B_link.petscvec"; 
  solver_Blink_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_B_link";          

  solver_Alink_input_params.sys_matrix_file = m_input_params.matrix_A.mass_tilde; 
  solver_Alink_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_A_link.petscvec"; 
  solver_Alink_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_A_link";          

  carl::print_input_params(m_ext_solver_Afree_input_filename,solver_Afree_input_params);

  carl::print_input_params(m_ext_solver_Bfree_input_filename, solver_Bfree_input_params); 

  carl::print_input_params(m_ext_solver_coupling_input_filename, solver_coupling_input_params); 

  carl::print_input_params(m_ext_solver_Blink_input_filename, solver_Blink_input_params); 

  carl::print_input_params(m_ext_solver_Alink_input_filename, solver_Alink_input_params); 

  
}

void Dyn_DI_Solver_Files_Setup::generate_libmesh_external_solver_script()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");

  m_ext_solver_Afree_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Afree_acc.slurm";
  m_ext_solver_Bfree_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Bfree_acc.slurm"; 
  m_ext_solver_coupling_script_filename = m_input_params.scratch_folder_path + "/ext_solver_coupling.slurm"; 
  m_ext_solver_Blink_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Blink_acc.slurm"; 
  m_ext_solver_Alink_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Alink_acc.slurm"; 

  switch (m_input_params.scheduler)
  {
    case ClusterSchedulerType::LOCAL :  homemade_error_msg("Scheduler not implemented yet!");
            break;

    case ClusterSchedulerType::PBS :    homemade_error_msg("Scheduler not implemented yet!");
            break;

    case ClusterSchedulerType::SLURM :  this->generate_libmesh_external_solver_scripts_SLURM();
            break;
    default : homemade_error_msg("Invalid scheduler name!");
  }

  m_bSetExternalSolversScripts = true;
}

void Dyn_DI_Solver_Files_Setup::generate_libmesh_external_solver_scripts_SLURM(){
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

    // Set the coupling scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_coupling.txt";                            
    slurm_error = m_input_params.scratch_folder_path + "/error_coupling.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_general + " " + m_ext_solver_coupling_input_filename; 


    this->print_SLURM_script(m_ext_solver_coupling_script_filename, "coupling",
             slurm_output, slurm_error, common_script,
             command_to_run);   

    // Set the Blink scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Blink_acc.txt";                               
    slurm_error = m_input_params.scratch_folder_path + "/error_Blink_acc.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_B + " " + m_ext_solver_Blink_input_filename;     

    this->print_SLURM_script(m_ext_solver_Blink_script_filename, "Blink_acc", 
            slurm_output, slurm_error, common_script,
            command_to_run); 

    // Set the Alink scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Alink_acc.txt";                               
    slurm_error = m_input_params.scratch_folder_path + "/error_Alink_acc.txt";                                
    command_to_run = m_input_params.ext_solver_launch_script_A + " " + m_ext_solver_Alink_input_filename;     

    this->print_SLURM_script(m_ext_solver_Alink_script_filename, "Alink_acc", 
            slurm_output, slurm_error, common_script,
            command_to_run); 
  }
}

void Dyn_DI_Solver_Files_Setup::generate_inner_operation_script()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");


  m_inner_operation_Afree_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Afree.slurm";
  m_inner_operation_Bfree_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Bfree.slurm";           
  m_inner_operation_coupling_script_filename = m_input_params.scratch_folder_path + "/setup_interpolation.slurm"; 
  m_inner_operation_pre_Blink_script_filename = m_input_params.scratch_folder_path + "/prepare_Blink.slurm"; 
  m_inner_operation_pre_Alink_script_filename = m_input_params.scratch_folder_path + "/prepare_Alink.slurm";
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

void Dyn_DI_Solver_Files_Setup::generate_inner_operation_scripts_SLURM(){
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
    slurm_output = m_input_params.scratch_folder_path + "/output_coupling.txt";                             
    slurm_error = m_input_params.scratch_folder_path + "/error_coupling.txt";                               
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_coupling -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_coupling_script_filename, "set_coupling",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the coupling finish scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_preBlink.txt";                             
    slurm_error = m_input_params.scratch_folder_path + "/error_preBlink.txt";                               
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_pre_Blink -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_pre_Blink_script_filename, "pre_Blink",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the coupling iterate scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_preAlink.txt";                             
    slurm_error = m_input_params.scratch_folder_path + "/error_preAlink.txt";                               
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_pre_Alink -i "+ m_input_params.general_entry_file_path;  

    this->print_SLURM_script(m_inner_operation_pre_Alink_script_filename, "pre_Alink",
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

void Dyn_DI_Solver_Files_Setup::generate_combined_scripts(){
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");
  homemade_assert_msg(m_bSetInnerOperationScripts,"Inner operation script files not set yet!");

  m_Afree_combined_filename = m_input_params.scratch_folder_path + "/Afree.sh";
  m_Bfree_combined_filename = m_input_params.scratch_folder_path + "/Bfree.sh";                     
  m_coupling_combined_filename = m_input_params.scratch_folder_path + "/coupling.sh";
  m_Blink_combined_filename = m_input_params.scratch_folder_path + "/Blink.sh";                  
  m_Alink_combined_filename = m_input_params.scratch_folder_path + "/Alink.sh";                     

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

void Dyn_DI_Solver_Files_Setup::generate_combined_scripts_SLURM(){
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

    // coupling_script
    std::ofstream coupling_script(m_coupling_combined_filename);
    coupling_script << "#!/bin/bash" << std::endl;
    coupling_script << std::endl;
    coupling_script << "job_setup=$(sbatch " << m_inner_operation_coupling_script_filename << ")" << std::endl;
    coupling_script << "job_interpolation=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name set_coupling) "<< m_ext_solver_coupling_script_filename << ")" << std::endl;
    coupling_script << "job_prepare=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name coupling) "<< m_inner_operation_pre_Blink_script_filename << ")" << std::endl;
    coupling_script.close();

    // Blink_script

    std::ofstream Blink_script(m_Blink_combined_filename);
    Blink_script << "#!/bin/bash" << std::endl;
    Blink_script << std::endl;
    Blink_script << "job_acc=$(sbatch " << m_ext_solver_Blink_script_filename << ")" << std::endl;
    Blink_script << "job_reste=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name Blink_acc) "<< m_inner_operation_Blink_script_filename << ")" << std::endl;
    Blink_script.close();


    std::ofstream Alink_script(m_Alink_combined_filename);
    Alink_script << "#!/bin/bash" << std::endl;
    Alink_script << std::endl;
    Alink_script << "job_prepare=$(sbatch "<< m_inner_operation_pre_Alink_script_filename << ")" << std::endl;
    Alink_script << "job_acc=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name pre_Alink) " << m_ext_solver_Alink_script_filename << ")" << std::endl;
    Alink_script << "job_reste=$(sbatch --dependency=afterok:$(squeue --noheader --format %i --name Alink_acc) "<< m_inner_operation_Alink_script_filename << ")" << std::endl;
    Alink_script.close();


  }
}


void Dyn_DI_Solver_Files_Setup::generate_progression_inputs(){
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
