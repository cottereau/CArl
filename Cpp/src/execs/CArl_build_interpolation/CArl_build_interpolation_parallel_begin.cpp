/*
 *
 *  Created on: Oct 12, 2022
 *      Author: Chensheng Luo
 */

#include "CArl_build_interpolation.h"
/** \file CArl_build_interpolation_parallel_begin.cpp
 * \brief **DYN-DI** Program used to calculate interpolation matrix by parallel calculation.

 Usage: `./CArl_build_interpolation_matrix_parallel_begin -i [input file]`

The program CArl_build_interpolation_parallel_begin and CArl_build_interpolation_parallel_end will calculate the interpolation matrix \f$H=\beta^A (\Delta t^A)^2 C^A (M^A)^{-1} (C^A)^T+\beta^B (\Delta t^B)^2 C^B (M^B)^{-1} (C^B)^T \f$. 
This calculation is done in following three steps:
  1. Separating \f$ (C^i)^T \f$ into several columns: C^i_{ttt} (begin)
  2. Inverting by solving \f$ M^i X^i_{ttt}=C^i_{ttt} \f$ (solving)
  3. By putting \f$X^i_{ttt}\f$ together to form a matrix \f$G^i\f$ and then calculate \f$H=\beta^A (\Delta t^A)^2 C^A G^A+\beta^B (\Delta t^B)^2 C^B G^B \f (end)$

The input file is parsed by the get_input_params(GetPot& field_parser, carl_build_interpolation_params& input_params) function. 

 */

void generate_external_solver_inputs(std::string matrix_name,
          std::string vector_name,
          std::string solution_output_name,
          std::string template_name,
          std::string file_name);

void generate_external_solver_script(libMesh::Parallel::Communicator& m_comm,
          carl::ClusterSchedulerType scheduler,
          std::string script_save_filename,
          std::string template_name,
          std::string output_log_file,
          std::string error_file,
          std::string input_filename,
          std::string job_name,
          std::string combined_file_name,
          std::string external_solver);

void generate_external_solver_scripts_SLURM(libMesh::Parallel::Communicator& m_comm,
          std::string script_save_filename,
          std::string template_name,
          std::string slurm_output,
          std::string slurm_error,
          std::string input_filename,
          std::string job_name,
          std::string combined_file_name,
          std::string external_solver);

void print_SLURM_script(const std::string& output_filename, 
          const std::string& job_name, 
          const std::string& output_name, 
          const std::string& error_name, 
          const std::string& common_script, 
          const std::string& command_to_run);

int main(int argc, char *argv[]) 
{
  // --- Setup part
  // - Initialize libMesh
  libMesh::LibMeshInit init(argc, argv);

  // libMesh's C++ / MPI communicator wrapper
  libMesh::Parallel::Communicator& WorldComm = init.comm();

  // Number of processors and processor rank.
  int rank = WorldComm.rank();
  int nodes = WorldComm.size();

  libMesh::PerfLog perf_log("Main program");

  // -  Set up input_params
  // Command line parser
  GetPot command_line(argc, argv);

  // File parser
  GetPot field_parser;

  // If there is an input file, parse it to get the parameters. Else, parse the command line
  std::string input_filename;
  if (command_line.search(2, "--inputfile", "-i")) 
  {
    input_filename = command_line.next(input_filename);
    field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
  } 
  else
  {
    field_parser = command_line;
  }

  // Parse file content to the input_params
  carl_build_interpolation_params input_params;
  get_input_params(field_parser, input_params);

  if(input_params.m_bParallelPossible==false){
    homemade_error_msg("[CArl Error] Parallel Not Possible, check your parameter!");
  }

  // - Definition of matrix

  Mat coupling_matrix_A;
  Mat coupling_matrix_B;
  Mat transpose_coupling_matrix_A;
  Mat transpose_coupling_matrix_B;
  Vec output_vec_A;
  Vec output_vec_B;
  PetscInt CmA,CnA,CmB,CnB;
  // - Load matrix

  perf_log.push("Load PETSc Matrices");

  MatCreate(WorldComm.get(),&coupling_matrix_A);
  MatCreate(WorldComm.get(),&coupling_matrix_B);
  MatCreate(WorldComm.get(),&transpose_coupling_matrix_A);
  MatCreate(WorldComm.get(),&transpose_coupling_matrix_B);

  carl::read_PETSC_matrix(coupling_matrix_A, 
                input_params.path_macro_coupling_matrix, WorldComm.get());

  carl::read_PETSC_matrix(coupling_matrix_B, 
                input_params.path_micro_coupling_matrix, WorldComm.get());

  perf_log.pop("Load PETSc Matrices");

  // - Transpose matrix

  perf_log.push("Transpose matrix");

  MatGetSize(coupling_matrix_A, &CnA, &CmA);
  MatGetSize(coupling_matrix_B, &CnB, &CmB);
  MatSetSizes(transpose_coupling_matrix_A, PETSC_DECIDE, PETSC_DECIDE, CmA, CnA);
  MatSetSizes(transpose_coupling_matrix_B, PETSC_DECIDE, PETSC_DECIDE, CmB, CnB);
  MatTranspose(coupling_matrix_A, MAT_INITIAL_MATRIX, &transpose_coupling_matrix_A);
  MatTranspose(coupling_matrix_B, MAT_INITIAL_MATRIX, &transpose_coupling_matrix_B);
  
  perf_log.pop("Transpose matrix");

  // - Make folders

  perf_log.push("Make folders");

  std::string command_string;

  if(WorldComm.rank() == 0)
  {

    command_string = "rm -rf " + input_params.output_base + "coupling_transpose_A/";
    std::cout << command_string << std::endl;
    carl::exec_command(command_string.c_str());

    command_string = "mkdir -p " + input_params.output_base + "coupling_transpose_A/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "rm -rf " + input_params.output_base + "coupling_transpose_B/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "mkdir -p " + input_params.output_base + "coupling_transpose_B/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "rm -rf " + input_params.output_base + "inverse_lancer/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "mkdir -p " + input_params.output_base + "inverse_lancer/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "rm -rf " + input_params.output_base + "mass_inverted_coupling_A/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "mkdir -p " + input_params.output_base + "mass_inverted_coupling_A/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "rm -rf " + input_params.output_base + "mass_inverted_coupling_B/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "mkdir -p " + input_params.output_base + "mass_inverted_coupling_B/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "rm -rf " + input_params.output_base + "outputs/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;

    command_string = "mkdir -p " + input_params.output_base + "outputs/";
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;
  }
  
  perf_log.pop("Make folders");

  // - Output all columns vectors 

  perf_log.push("Output vectors");


  VecCreate(WorldComm.get(), &output_vec_A);
  VecSetSizes(output_vec_A, PETSC_DECIDE, CmA);
  VecSetType(output_vec_A, VECMPI);
  PetscInt i;

  for (i = 0; i < CnA; i++){
    MatGetColumnVector(transpose_coupling_matrix_A, output_vec_A, i);
    carl::write_PETSC_vector(output_vec_A, input_params.output_base + "coupling_transpose_A/CTA_"+std::to_string(i)+".petscvec",0,WorldComm.get(),1);
  }

  VecCreate(WorldComm.get(), &output_vec_B);
  VecSetSizes(output_vec_B, PETSC_DECIDE, CmB);
  VecSetType(output_vec_B, VECMPI);

  for (i = 0; i < CnB; i++){
    MatGetColumnVector(transpose_coupling_matrix_B, output_vec_B, i);
    carl::write_PETSC_vector(output_vec_B, input_params.output_base + "coupling_transpose_B/CTB_"+std::to_string(i)+".petscvec",0,WorldComm.get(),1);
  }
  
  perf_log.pop("Output vectors");

  // - Prepare all potential lancers and inputs

  perf_log.push("Prepare lancers and inputs");

  if(WorldComm.rank() == 0)
  {
    std::ofstream combined_script(input_params.output_base + "inverse_lancer/mass_inverse_all.sh");
    combined_script << "#!/bin/bash" << std::endl;
    combined_script << std::endl;
    combined_script.close();
  }

  for (i = 0; i < CnA; i++){
    generate_external_solver_inputs(input_params.path_tilde_matrix_A,
          input_params.output_base + "coupling_transpose_A/CTA_"+std::to_string(i)+".petscvec",
          input_params.output_base + "mass_inverted_coupling_A/MCTA_"+std::to_string(i)+".petscvec",
          input_params.ext_solver_A_input,
          input_params.output_base + "inverse_lancer/mass_inverse_input_A_"+std::to_string(i)+".txt");
    
    generate_external_solver_script(WorldComm,
          input_params.scheduler,
          input_params.output_base + "inverse_lancer/inverse_mass_A_"+std::to_string(i)+".slurm",
          input_params.script_filename,
          input_params.output_base + "outputs/output_A_"+std::to_string(i)+".out",
          input_params.output_base + "outputs/error_A_"+std::to_string(i)+".out",
          input_params.output_base + "inverse_lancer/mass_inverse_input_A_"+std::to_string(i)+".txt",
          "invA"+std::to_string(i),
          input_params.output_base + "inverse_lancer/mass_inverse_all.sh",
          input_params.ext_solver_launch_script_A);
  }

  for (i = 0; i < CnB; i++){
    generate_external_solver_inputs(input_params.path_tilde_matrix_B,
          input_params.output_base + "coupling_transpose_B/CTB_"+std::to_string(i)+".petscvec",
          input_params.output_base + "mass_inverted_coupling_B/MCTB_"+std::to_string(i)+".petscvec",
          input_params.ext_solver_B_input,
          input_params.output_base + "inverse_lancer/mass_inverse_input_B_"+std::to_string(i)+".txt");
    
    generate_external_solver_script(WorldComm,
          input_params.scheduler,
          input_params.output_base + "inverse_lancer/inverse_mass_B_"+std::to_string(i)+".slurm",
          input_params.script_filename,
          input_params.output_base + "outputs/output_B_"+std::to_string(i)+".out",
          input_params.output_base + "outputs/error_B_"+std::to_string(i)+".out",
          input_params.output_base + "inverse_lancer/mass_inverse_input_B_"+std::to_string(i)+".txt",
          "invB"+std::to_string(i),
          input_params.output_base + "inverse_lancer/mass_inverse_all.sh",
          input_params.ext_solver_launch_script_B);
  }

  // TODO: add final parallel end!
  //if(WorldComm.rank() == 0)
  //{
    //std::ofstream combined_script(input_params.output_base + "inverse_lancer/mass_inverse_all.sh",std::ios_base:: app);
    //combined_script.close();
  //}
  
  perf_log.pop("Prepare lancers and inputs");

  // - Execute files

  perf_log.push("Execute files");
  if(WorldComm.rank() == 0){
     std::string script_command = ". " + input_params.output_base + "inverse_lancer/mass_inverse_all.sh";
     if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
     {
        std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
        std::cout << script_command << std::endl ;
     } else {
        carl::exec_command(script_command);
     }
  }

  
  perf_log.pop("Execute files");

  // - Clean up

  perf_log.push("Clean up");
  MatDestroy(&coupling_matrix_A);
  MatDestroy(&coupling_matrix_B);
  VecDestroy(&output_vec_A);
  VecDestroy(&output_vec_B);

  perf_log.pop("Clean up");
  
  return 0;
}

void generate_external_solver_inputs(std::string matrix_name,
          std::string vector_name,
          std::string solution_output_name,
          std::string template_name,
          std::string file_name)
{ 
  // Get general input parameters
  GetPot field_parser;
  field_parser.parse_input_file(template_name, "#", "\n", " \t\n");

  // Get general input parameters
  carl::libmesh_solve_linear_system_input_params solver_input_params;
  carl::get_input_params(field_parser, solver_input_params);

  // Set input parameters
  solver_input_params.sys_matrix_file = matrix_name;
  solver_input_params.sys_rhs_vec_file = vector_name;
  solver_input_params.output_base = solution_output_name;

  carl::print_input_params(file_name,solver_input_params);  
}

void generate_external_solver_script(libMesh::Parallel::Communicator& m_comm,
          carl::ClusterSchedulerType scheduler,
          std::string script_save_filename,
          std::string template_name,
          std::string output_log_file,
          std::string error_file,
          std::string input_filename,
          std::string job_name,
          std::string combined_file_name,
          std::string external_solver)
{
  switch (scheduler)
  {
    case carl::ClusterSchedulerType::LOCAL :  homemade_error_msg("Scheduler not implemented yet!");
            break;

    case carl::ClusterSchedulerType::PBS :    homemade_error_msg("Scheduler not implemented yet!");
            break;

    case carl::ClusterSchedulerType::SLURM :  generate_external_solver_scripts_SLURM(m_comm,
          script_save_filename,
          template_name,
          output_log_file,
          error_file,
          input_filename,
          job_name,
          combined_file_name,
          external_solver);
          break;
    default : homemade_error_msg("Invalid scheduler name!");
  }
};

void generate_external_solver_scripts_SLURM(libMesh::Parallel::Communicator& m_comm,
          std::string script_save_filename,
          std::string template_name,
          std::string slurm_output,
          std::string slurm_error,
          std::string input_filename,
          std::string job_name,
          std::string combined_file_name,
          std::string external_solver)
{
  if(m_comm.rank() == 0)
  {
    // Get the full common script file into a string
    std::ifstream base_script(template_name);
    std::string common_script((std::istreambuf_iterator<char>(base_script)),
                  std::istreambuf_iterator<char>());
    base_script.close();

    std::string command_to_run;

    command_to_run = external_solver + " " + input_filename;

    print_SLURM_script(script_save_filename, job_name, slurm_output, slurm_error, common_script,
              command_to_run);

    std::ofstream combined_script(combined_file_name,std::ios_base:: app);
    combined_script << "job=$(sbatch " << script_save_filename << ")" << std::endl;
    combined_script.close();
  }
};

void print_SLURM_script(const std::string& output_filename, 
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