#include "dyn_solver_files_setup.h"

namespace carl
{

void Dyn_Solver_Files_Setup::print_SLURM_script(const std::string& output_filename, const std::string& job_name, const std::string& output_name, const std::string& error_name, const std::string& common_script, const std::string& command_to_run)
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

void Dyn_Solver_Files_Setup::set_FETI_input_parameters(feti_loop_dyn_params& input_params)
{
  m_bInputParamsSet = true;
  m_input_params = input_params;
};

void Dyn_Solver_Files_Setup::set_scratch_folder()
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

void Dyn_Solver_Files_Setup::generate_libmesh_external_solver_inputs()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");

  m_ext_solver_Afree_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Afree.txt";
  m_ext_solver_Bfree_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Bfree.txt"; // sev
  m_ext_solver_couplingIV_input_filename = m_input_params.scratch_folder_path + "/ext_solver_interpolation.txt"; //sev
  m_ext_solver_Blink_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Blink.txt"; // sev
  m_ext_solver_Alink_input_filename = m_input_params.scratch_folder_path + "/ext_solver_Blink.txt"; //sev


// Get general input parameters
  GetPot field_parser_A;
  field_parser_A.parse_input_file(m_input_params.ext_solver_A_input, "#", "\n", " \t\n");

  GetPot field_parser_B;                                                                   //sev
  field_parser_B.parse_input_file(m_input_params.ext_solver_B_input, "#", "\n", " \t\n");  //sev

  GetPot field_parser_couplingIV;                                                                        //sev
  field_parser_couplingIV.parse_input_file(m_input_params.ext_solver_general_input, "#", "\n", " \t\n"); //sev

// Get general input parameters
  carl::libmesh_solve_linear_system_input_params solver_Afree_input_params;
  carl::get_input_params(field_parser_A, solver_Afree_input_params);

  carl::libmesh_solve_linear_system_input_params solver_Bfree_input_params;   //sev
  carl::get_input_params(field_parser_B, solver_Bfree_input_params);          //sev

  carl::libmesh_solve_linear_system_input_params solver_couplingIV_input_params;   //sev
  carl::get_input_params(field_parser_couplingIV, solver_couplingIV_input_params); //sev

  carl::libmesh_solve_linear_system_input_params solver_Blink_input_params;   // sev
  carl::get_input_params(field_parser_B, solver_Blink_input_params);          // sev

  carl::libmesh_solve_linear_system_input_params solver_Alink_input_params;   // sev
  carl::get_input_params(field_parser_B, solver_Alink_input_params);          // sev

  // Set input parameters
  solver_Afree_input_params.sys_matrix_file = m_input_params.tilde_matrix_A;
  solver_Afree_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_A_free.petscvec";
  solver_Afree_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_A_free";

  solver_Bfree_input_params.sys_matrix_file = m_input_params.tilde_matrix_B;                                      //sev
  solver_Bfree_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_B_free.petscvec"; //sev
  solver_Bfree_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_B_free";                 //sev

  solver_couplingIV_input_params.sys_matrix_file = m_input_params.interpolation_matrix;                                     //sev
  solver_couplingIV_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_interpolation_vec.petscvec"; //sev
  solver_couplingIV_input_params.output_base = m_input_params.scratch_folder_path + "/coupling_interpolation_vec";          //sev
 
  solver_Blink_input_params.sys_matrix_file = m_input_params.tilde_matrix_B;                                      //sev
  solver_Blink_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_B_link.petscvec"; //sev
  solver_Blink_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_B_link";                 //sev

  solver_Blink_input_params.sys_matrix_file = m_input_params.tilde_matrix_A;                                      //sev
  solver_Blink_input_params.sys_rhs_vec_file = m_input_params.scratch_folder_path + "/rhs_vec_A_link.petscvec"; //sev
  solver_Blink_input_params.output_base = m_input_params.scratch_folder_path + "/this_acc_A_link";                 //sev


  carl::print_input_params(m_ext_solver_Afree_input_filename,solver_Afree_input_params);

  carl::print_input_params(m_ext_solver_Bfree_input_filename, solver_Bfree_input_params); //sev

  carl::print_input_params(m_ext_solver_couplingIV_input_filename, solver_couplingIV_input_params); //sev

  carl::print_input_params(m_ext_solver_Blink_input_filename, solver_Blink_input_params); //sev

  carl::print_input_params(m_ext_solver_Alink_input_filename, solver_Alink_input_params); //sev


  m_bSetExternalSolversInputFiles = true;
}

void Dyn_Solver_Files_Setup::generate_libmesh_external_solver_script()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");

  m_ext_solver_Afree_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Afree_acc.sh";
  m_ext_solver_Bfree_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Bfree_acc.sh"; //sev
  m_ext_solver_couplingIV_script_filename = m_input_params.scratch_folder_path + "/ext_solver_coupling.sh"; //sev
  m_ext_solver_Blink_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Blink_acc.sh"; //sev
  m_ext_solver_Alink_script_filename = m_input_params.scratch_folder_path + "/ext_solver_Alink_acc.sh"; //sev


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

void Dyn_Solver_Files_Setup::generate_libmesh_external_solver_scripts_SLURM(){
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
    command_to_run = m_input_params.ext_solver_launch_script + " " + m_ext_solver_Afree_input_filename;

    this->print_SLURM_script(m_ext_solver_Afree_script_filename, "Afree_acc",
              slurm_output, slurm_error, common_script,
              command_to_run);

    // Set the Bfree_acc scripts
    slurm_output = m_input_params.scratch_folder_path + "/out_Bfree_acc.txt";                               //sev
    slurm_error = m_input_params.scratch_folder_path + "/error_Bfree_acc.txt";                                //sev
    command_to_run = m_input_params.ext_solver_launch_script + " " + m_ext_solver_Bfree_input_filename;     //sev

    this->print_SLURM_script(m_ext_solver_Bfree_input_filename, "Bfree_acc", 
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                //sev

    // Set the couplingIV scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_coupling.txt";                            //sev
    slurm_error = m_input_params.scratch_folder_path + "/error_coupling.txt";                                //sev
    command_to_run = m_input_params.ext_solver_launch_script + " " + m_ext_solver_couplingIV_input_filename; //sev


    this->print_SLURM_script(m_ext_solver_couplingIV_input_filename, "coupling",
             slurm_output, slurm_error, common_script,
             command_to_run);                                                                                //sev

    // Set the Blink scripts
    slurm_output = m_input_params.scratch_folder_path + "/out_Blink_acc.txt";                               //sev
    slurm_error = m_input_params.scratch_folder_path + "/error_Blink_acc.txt";                                //sev
    command_to_run = m_input_params.ext_solver_launch_script + " " + m_ext_solver_Blink_input_filename;     //sev

    this->print_SLURM_script(m_ext_solver_Blink_input_filename, "Blink_acc", 
            slurm_output, slurm_error, common_script,
            command_to_run);                                                                                // sev

    // Set the Alink scripts
    slurm_output = m_input_params.scratch_folder_path + "/out_Alink_acc.txt";                               //sev
    slurm_error = m_input_params.scratch_folder_path + "/error_Alink_acc.txt";                                //sev
    command_to_run = m_input_params.ext_solver_launch_script + " " + m_ext_solver_Alink_input_filename;     //sev

    this->print_SLURM_script(m_ext_solver_Alink_input_filename, "Alink_acc", 
            slurm_output, slurm_error, common_script,
            command_to_run);

  }
}

void Dyn_Solver_Files_Setup::generate_inner_operation_script()
{
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");


  m_inner_operation_Afree_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Afree.sh";
  m_inner_operation_Bfree_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Bfree.sh";           //sev
  m_inner_operation_setup_interpolation_script_filename = m_input_params.scratch_folder_path + "/setup_interpolation.sh"; //sev
  m_inner_operation_prepare_Blink_script_filename = m_input_params.scratch_folder_path + "/prepare_Blink.sh"; //sev
  m_inner_operation_Blink_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Blink.sh";       //sev
  m_inner_operation_prepare_Alink_script_filename = m_input_params.scratch_folder_path + "/prepare_Alink";    //sev
  m_inner_operation_Alink_script_filename = m_input_params.scratch_folder_path + "/inner_ope_Alink.sh";       //sev



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

void Dyn_Solver_Files_Setup::generate_inner_operation_scripts_SLURM(){
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
    slurm_output = m_input_params.scratch_folder_path + "/output_Bfree_speed.txt";                            //sev
    slurm_error  = m_input_params.scratch_folder_path + "/error_Bfree_speed.txt";                             //sev
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_Bfree -i "+ m_input_params.general_entry_file_path;  //sev

    this->print_SLURM_script(m_inner_operation_Bfree_script_filename, "Bfree_speed",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                //sev

    // Set the couplingIV scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_setup_interpolation.txt";                             //sev
    slurm_error = m_input_params.scratch_folder_path + "/error_setup_interpolation.txt";                               //sev
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_inter -i "+ m_input_params.general_entry_file_path;  //sev

    this->print_SLURM_script(m_inner_operation_setup_interpolation_script_filename, "",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                    //sev

    slurm_output = m_input_params.scratch_folder_path + "output_pre_Blink.txt";                                   //sev
    slurm_error = m_input_params.scratch_folder_path + "/error_pre_Blink.txt";                                    //sev
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_pre_Blink -i "+ m_input_params.general_entry_file_path;  //sev

    this->print_SLURM_script(m_inner_operation_prepare_Blink_script_filename, "",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                    //sev

    // Set the Blink scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Blink_speed.txt";                            //sev
    slurm_error  = m_input_params.scratch_folder_path + "/error_Blink_speed.txt";                             //sev
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_Blink -i "+ m_input_params.general_entry_file_path;  //sev

    this->print_SLURM_script(m_inner_operation_Blink_script_filename, "Blink_speed",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                //sev

    slurm_output = m_input_params.scratch_folder_path + "output_pre_Alink.txt";                                   //sev
    slurm_error = m_input_params.scratch_folder_path + "/error_pre_Alink.txt";                                    //sev
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_pre_Alink -i "+ m_input_params.general_entry_file_path;  //sev

    this->print_SLURM_script(m_inner_operation_prepare_Alink_script_filename, "",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                    //sev

    // Set the Alink scripts
    slurm_output = m_input_params.scratch_folder_path + "/output_Alink_speed.txt";                            //sev
    slurm_error  = m_input_params.scratch_folder_path + "/error_Alink_speed.txt";                             //sev
    command_to_run = "srun -n 4 $CARLBUILD/CArl_loop_dyn_Alink -i "+ m_input_params.general_entry_file_path;  //sev

    this->print_SLURM_script(m_inner_operation_Alink_script_filename, "Alink_speed",
              slurm_output, slurm_error, common_script,
              command_to_run);                                                                                //sev

  }
}

void Dyn_Solver_Files_Setup::generate_combined_scripts(){
  homemade_assert_msg(m_bInputParamsSet,"Input parameters not set yet!");
  homemade_assert_msg(m_bScratchFolderExists,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSetExternalSolversInputFiles,"External solver input files not set yet!");
  homemade_assert_msg(m_bSetExternalSolversScripts,"External solver script files not set yet!");
  homemade_assert_msg(m_bSetInnerOperationScripts,"Inner operation script files not set yet!");

  m_Afree_combined_filename = m_input_params.scratch_folder_path + "/Afree.sh";
  m_Bfree_combined_filename = m_input_params.scratch_folder_path + "/Bfree.sh";                     //sev
  m_coupling_combined_filename = m_input_params.scratch_folder_path + "/coupling.sh";
  m_Blink_combined_filename = m_input_params.scratch_folder_path + "/coupling.sh";                  //sev
  m_Alink_combined_filename = m_input_params.scratch_folder_path + "/Alink.sh";                     //sev

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

void Dyn_Solver_Files_Setup::generate_combined_scripts_SLURM(){
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
    Afree_script << "job_acc=`sbatch " << m_ext_solver_Afree_script_filename << "`" << std::endl;
    Afree_script << "job_reste=`sbatch --dependency=afterok:$job_acc"<< m_inner_operation_Afree_script_filename << "`" << std::endl;
    Afree_script.close();

    // Bfree_script
    std::ofstream Bfree_script(m_Bfree_combined_filename);
    Bfree_script << "#!/bin/bash" << std::endl;
    Bfree_script << std::endl;
    Bfree_script << "job_acc=`sbatch " << m_ext_solver_Bfree_script_filename << "`" << std::endl;
    Bfree_script << "job_reste=`sbatch --dependency=afterok:$job_acc"<< m_inner_operation_Bfree_script_filename << "`" << std::endl;
    Bfree_script.close();

    //coupling_script
    std::ofstream coupling_script(m_coupling_combined_filename);
    coupling_script << "#!/bin/bash" << std::endl;
    coupling_script << std::endl;
    coupling_script << "job_setup=`sbatch " << m_inner_operation_setup_interpolation_script_filename << "`" << std::endl;
    coupling_script << "job_interpolation=`sbatch --dependency=afterok:$job_setup"<< m_ext_solver_couplingIV_script_filename << "`" << std::endl;
    coupling_script << "job_prepare=`sbatch --dependency=afterok:$job_interpolation" << m_inner_operation_prepare_Blink_script_filename << "`" << std::endl;
    coupling_script.close();

    // Blink_script
    std::ofstream Blink_script(m_Blink_combined_filename);                               //Sev
    Blink_script << "#!/bin/bash" << std::endl;                                          //sev
    Blink_script << std::endl;                                                            //sev
    Blink_script << "job_acc=`sbatch " << m_ext_solver_Blink_script_filename << "`" << std::endl; //Sev
    Blink_script << "job_reste=`sbatch --dependency=afterok:$job_acc"<< m_inner_operation_Blink_script_filename << "`" << std::endl; //Sev
    Blink_script.close();                                                                  // sev

    // Alink_script
    std::ofstream Alink_script(m_Alink_combined_filename);                               // Sev
    Alink_script << "#!/bin/bash" << std::endl;                                          // sev
    Alink_script << std::endl;                                                               //sev
    Alink_script << "job_prepare=`sbatch --dependency=afterok:$job_interpolation" << m_inner_operation_prepare_Alink_script_filename << "`" << std::endl;    //sev
    Alink_script << "job_acc=`sbatch " << m_ext_solver_Alink_script_filename << "`" << std::endl;                                    //sev
    Alink_script << "job_reste=`sbatch --dependency=afterok:$job_acc"<< m_inner_operation_Alink_script_filename << "`" << std::endl; //sev
    Alink_script.close();                                                                //sev


  }
}


void Dyn_Solver_Files_Setup::generate_progression_inputs(){
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
  m_bSetProgressionInput=true;
  }

}





}