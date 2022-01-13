/*
 * carl_loop_dyn_input_parser.cpp
 *
 *  Created on: Nov 17, 2021
 *      Author: Chensheng Luo
 */

#include "carl_loop_dyn_input_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
    feti_loop_dyn_params& input_params) {

  if (field_parser.search(1, "ClusterSchedulerType")) {
    std::string cluster_scheduler_type;
    cluster_scheduler_type = field_parser.next(cluster_scheduler_type);
    if(cluster_scheduler_type == "LOCAL")
    {
      std::cout << " !!! WARNING: " << std::endl;
      std::cout << "        The LOCAL 'scheduler' type is only intended for small and fast test cases" << std::endl;
      std::cout << "     on computers without a job scheduler (PBS, SLURM). You will have to launch" << std::endl;
      std::cout << "     each script MANUALLY!!! Reason: MPI does not support recursive 'mpirun' calls" << std::endl;
      input_params.scheduler = carl::ClusterSchedulerType::LOCAL;
      input_params.script_filename = "";
    }
    else if(cluster_scheduler_type == "PBS") {
      input_params.scheduler = carl::ClusterSchedulerType::PBS;
      if (field_parser.search(1, "ScriptFile")) {
        input_params.script_filename = field_parser.next(input_params.script_filename);
      } else {
        homemade_error_msg("Missing the script file (needed for the PBS scheduler)!");
      }
    }
    else if(cluster_scheduler_type == "SLURM") {
      input_params.scheduler = carl::ClusterSchedulerType::SLURM;
      if (field_parser.search(1, "ScriptFile")) {
        input_params.script_filename = field_parser.next(input_params.script_filename);
      } else {
        homemade_error_msg("Missing the script file (needed for the SLURM scheduler)!");
      }
    }
    else
      homemade_error_msg("Invalid scheduler type!");
  } else {
    homemade_error_msg("Missing the scheduler type!");
  }

  if (field_parser.search(1, "ScratchFolderPath")) {
    input_params.scratch_folder_path = field_parser.next(
        input_params.scratch_folder_path);
    std::cout << input_params.scratch_folder_path << std::endl;
  } else {
    homemade_error_msg("Missing scratch folder path command line!");
  }

  if (field_parser.search(1, "ResultFolderPath")) {
    input_params.result_folder_path = field_parser.next(
        input_params.result_folder_path);
    std::cout << input_params.result_folder_path << std::endl;
  } else {
    homemade_error_msg("Missing result folder path command line!");
  }

  if (field_parser.search(1, "TildeMatrixA")) {
    input_params.tilde_matrix_A = field_parser.next(
        input_params.tilde_matrix_A);
    std::cout << input_params.tilde_matrix_A << std::endl;
  } else {
    homemade_error_msg("Missing mass tilde matrix A command line!");
  }

  if (field_parser.search(1, "CouplingMatrixA")) {
    input_params.coupling_matrix_A = field_parser.next(
        input_params.coupling_matrix_A);
    std::cout << input_params.coupling_matrix_A << std::endl;
  } else {
    homemade_error_msg("Missing coupling matrix A command line!");
  }

  if (field_parser.search(1, "StiffnessMatrixA")) {
    input_params.stiffness_matrix_A = field_parser.next(
        input_params.stiffness_matrix_A);
    std::cout << input_params.stiffness_matrix_A << std::endl;
  } else {
    homemade_error_msg("Missing stiffness matrix A command line!");
  }

  if (field_parser.search(1, "RHSVectorA")) {
    input_params.rhs_vector_A = field_parser.next(
        input_params.rhs_vector_A);
    std::cout << input_params.rhs_vector_A << std::endl;
  } else {
    homemade_error_msg("Missing RHS vector A command line!");
  }

  if (field_parser.search(1, "ForceFolderA")) {
    input_params.force_folder_A = field_parser.next(
        input_params.force_folder_A);
    std::cout << input_params.force_folder_A << std::endl;
  } else {
    homemade_error_msg("Missing force folder A command line!");
  }

    if (field_parser.search(1, "TildeMatrixB")) {
    input_params.tilde_matrix_B = field_parser.next(
        input_params.tilde_matrix_B);
    std::cout << input_params.tilde_matrix_B << std::endl;
  } else {
    homemade_error_msg("Missing mass tilde matrix B command line!");
  }

  if (field_parser.search(1, "CouplingMatrixB")) {
    input_params.coupling_matrix_B = field_parser.next(
        input_params.coupling_matrix_B);
    std::cout << input_params.coupling_matrix_B << std::endl;
  } else {
    homemade_error_msg("Missing coupling matrix B command line!");
  }

  if (field_parser.search(1, "StiffnessMatrixB")) {
    input_params.stiffness_matrix_B = field_parser.next(
        input_params.stiffness_matrix_B);
    std::cout << input_params.stiffness_matrix_B << std::endl;
  } else {
    homemade_error_msg("Missing stiffness matrix B command line!");
  }

  if (field_parser.search(1, "RHSVectorB")) {
    input_params.rhs_vector_B = field_parser.next(
        input_params.rhs_vector_B);
    std::cout << input_params.rhs_vector_B << std::endl;
  } else {
    homemade_error_msg("Missing RHS vector B command line!");
  }

  if (field_parser.search(1, "StiffnessMatrixB")) {
    input_params.stiffness_matrix_B = field_parser.next(
        input_params.stiffness_matrix_B);
    std::cout << input_params.stiffness_matrix_B << std::endl;
  } else {
    homemade_error_msg("Missing stiffness matrix B command line!");
  }

  if (field_parser.search(1, "ForceFolderB")) {
    input_params.force_folder_B = field_parser.next(
        input_params.force_folder_B);
    std::cout << input_params.force_folder_B << std::endl;
  } else {
    homemade_error_msg("Missing force folder B command line!");
  }

  if (field_parser.search(1, "ExtSolverInputA")) {
    input_params.ext_solver_A_input = field_parser.next(
        input_params.ext_solver_A_input);
    std::cout << input_params.ext_solver_A_input << std::endl;
  } else {
    homemade_error_msg("Missing ext solver A input command line!");
  }

  if (field_parser.search(1, "ExtSolverInputB")) {
    input_params.ext_solver_B_input = field_parser.next(
        input_params.ext_solver_B_input);
    std::cout << input_params.ext_solver_B_input << std::endl;
  } else {
    homemade_error_msg("Missing ext solver B input command line!");
  }

  if (field_parser.search(1, "ExtSolverInputGeneral")) {
    input_params.ext_solver_general_input = field_parser.next(
        input_params.ext_solver_general_input);
    std::cout << input_params.ext_solver_general_input << std::endl;
  } else {
    homemade_error_msg("Missing ext solver general input command line!");
  }

  if (field_parser.search(1, "ExtSolverLaunchScript")) {
    input_params.ext_solver_launch_script = field_parser.next(
        input_params.ext_solver_launch_script);
    std::cout << input_params.ext_solver_launch_script << std::endl;
  } else {
    homemade_error_msg("Missing ext solver launch script command line!");
  }

  if (field_parser.search(1, "GeneralEntryFilePath")) {
    input_params.general_entry_file_path = field_parser.next(
        input_params.general_entry_file_path);
    std::cout << input_params.general_entry_file_path << std::endl;
  } else {
    homemade_error_msg("Missing scratch folder path command line!");
  }


  if (field_parser.search(1, "InterpolationMatrix")) {
    input_params.interpolation_matrix = field_parser.next(
        input_params.interpolation_matrix);
    std::cout << input_params.interpolation_matrix << std::endl;
  } else {
    homemade_error_msg("Missing interpolation matrix command line!");
  }


  if (field_parser.search(1, "deltatA")) {
    input_params.deltatA = field_parser.next(input_params.deltatA);
  } else {
    input_params.deltatA = 0.25;
  }

  if (field_parser.search(1, "betaA")) {
    input_params.betaA = field_parser.next(input_params.betaA);
  } else {
    input_params.betaA = 0.25;
  }

  if (field_parser.search(1, "gammaA")) {
    input_params.gammaA = field_parser.next(input_params.gammaA);
  } else {
    input_params.gammaA = 0.5;
  }

  if (field_parser.search(1, "deltatB")) {
    input_params.deltatB = field_parser.next(input_params.deltatB);
  } else {
    input_params.deltatB = 0.05;
  }

  if (field_parser.search(1, "betaB")) {
    input_params.betaB = field_parser.next(input_params.betaB);
  } else {
    input_params.betaB = 0.25;
  }

  if (field_parser.search(1, "gammaB")) {
    input_params.gammaB = field_parser.next(input_params.gammaB);
  } else {
    input_params.gammaB = 0.5;
  }

  if (field_parser.search(1, "SimulationDuration")) {
    input_params.simulation_duration = field_parser.next(input_params.simulation_duration);
  } else {
    input_params.simulation_duration = 10;
  }

  input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_OPERATOR;
  if ( field_parser.search(1, "CGPreconditionerType") )
  {
    std::string CG_precond_type_string = field_parser.next(CG_precond_type_string);
    if(CG_precond_type_string == "NONE")
      input_params.CG_precond_type = carl::BaseCGPrecondType::NO_PRECONDITIONER;
    else if(CG_precond_type_string == "Coupling_operator")
      input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_OPERATOR;
    else if(CG_precond_type_string == "Coupling_operator_jacobi")
      input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_JACOBI;
  }

  input_params.inner_loop_times=(int)(input_params.deltatA/input_params.deltatB);
  input_params.outer_loop_times=(int)(input_params.simulation_duration/input_params.deltatA);
  };
}