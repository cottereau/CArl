/*
 * carl_build_interpolation_input_parser.cpp
 *
 *
 *  Created on: Sep 26, 2021
 *      Author: Chensheng Luo
 */

#include "carl_build_interpolation_input_parser.h"

void get_input_params(GetPot& field_parser,carl_build_interpolation_params& input_params) 
{
  input_params.m_bParallelPossible = true;

  //Set mass tilde matrix for A
  if (field_parser.search(1, "mTildeA")) {
    input_params.path_tilde_matrix_A = field_parser.next(input_params.path_tilde_matrix_A);
  } else {
    homemade_error_msg("[CArl Parameters]Missing mass tilde file!");
  }

  //Set mass tilde matrix for B
  if (field_parser.search(1, "mTildeB")) {
    input_params.path_tilde_matrix_B = field_parser.next(input_params.path_tilde_matrix_B);
  } else {
    homemade_error_msg("[CArl Parameters]Missing mass tilde file!");
  }

  //Set coupling matrix for A
  if (field_parser.search(1, "couplingMatrixA")) {
    input_params.path_macro_coupling_matrix = field_parser.next(input_params.path_macro_coupling_matrix);
  } else {
    homemade_error_msg("[CArl Parameters]Missing coupling matrix for A file!");
  }

  //Set coupling matrix for B
  if (field_parser.search(1, "couplingMatrixB")) {
    input_params.path_micro_coupling_matrix = field_parser.next(input_params.path_micro_coupling_matrix);
  } else {
    homemade_error_msg("[CArl Parameters]Missing coupling matrix for B file!");
  }

  //Set output folder
  if (field_parser.search(1, "OutputFolder")) {
    input_params.output_base = field_parser.next(input_params.output_base);
  } else {
    homemade_error_msg("[CArl Parameters]Missing output base!");
  }

    //Set output folder
  if (field_parser.search(1, "InterInputFile")) {
    input_params.interpolation_input_file = field_parser.next(input_params.interpolation_input_file);
  } else {
    homemade_error_msg("[CArl Parameters]Missing interpolation input file!");
  }


  if (field_parser.search(2, "NewmarkParametersA","NewmarkParameters")){
      std::string filename;
      filename = field_parser.next(filename);
      GetPot newmark_parser;
      newmark_parser.parse_input_file(filename, "#", "\n", " \t\n");
      carl::get_newmark_params(newmark_parser, input_params.newmark_A);
    } else{
      homemade_error_msg("[CArl Parameters]Missing A Newmark parameters file!");
    }

  if (field_parser.search(2, "NewmarkParametersB","NewmarkParameters")){
      std::string filename;
      filename = field_parser.next(filename);
      GetPot newmark_parser;
      newmark_parser.parse_input_file(filename, "#", "\n", " \t\n");
      carl::get_newmark_params(newmark_parser, input_params.newmark_B);
    } else{
      homemade_error_msg("[CArl Parameters]Missing B Newmark parameters file!");
    }

  if (field_parser.search(1, "InverseNameA")) {
    input_params.path_inv_matrix_A = field_parser.next(input_params.path_inv_matrix_A);
  } else {
    input_params.path_inv_matrix_A = "inv_dense_A";
    std::cout << " [CArl Parameters]WARNING! No invert mass matrix A output base!, default number is chosen as "<< input_params.path_inv_matrix_A << std::endl;
  }

  if (field_parser.search(1, "InverseNameB")) {
    input_params.path_inv_matrix_B = field_parser.next(input_params.path_inv_matrix_B);
  } else {
    input_params.path_inv_matrix_B = "inv_dense_B";
    std::cout << " [CArl Parameters]WARNING! No invert mass matrix B output base!, default number is chosen as "<< input_params.path_inv_matrix_B << std::endl;
  }

  if (field_parser.search(1, "InvertThreshold")){
    input_params.invert_threshold = field_parser.next(input_params.invert_threshold);
  } else {
    input_params.invert_threshold = + std::numeric_limits<double>::infinity(); // No threshold is set!
    std::cout << " [CArl Parameters]WARNING! No inversion norm residual threshold is chosen! " << std::endl;
    std::cout << " Please check manually result_test_inversion.txt in your output base." << std::endl;
  }

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
    input_params.scheduler = carl::ClusterSchedulerType::NONE;
    input_params.script_filename = "";
    input_params.m_bParallelPossible = false;
    std::cout << " [CArl Parameters]WARNING! No cluster scheduler type, parallel inversion not possible!" << std::endl;
  }

  if (field_parser.search(2, "ExtSolverInputA", "ExtSolverInput")) {
    input_params.ext_solver_A_input = field_parser.next(
        input_params.ext_solver_A_input);
    std::cout << input_params.ext_solver_A_input << std::endl;
  } else {
    input_params.ext_solver_A_input = "";
    input_params.m_bParallelPossible = false;
    std::cout << " [CArl Parameters]WARNING! No external solver A input, parallel inversion not possible!" << std::endl;
  }

  if (field_parser.search(2, "ExtSolverInputB", "ExtSolverInput")) {
    input_params.ext_solver_B_input = field_parser.next(
        input_params.ext_solver_B_input);
    std::cout << input_params.ext_solver_B_input << std::endl;
  } else {
    input_params.ext_solver_B_input = "";
    input_params.m_bParallelPossible = false;
    std::cout << " [CArl Parameters]WARNING! No external solver B input, parallel inversion not possible!" << std::endl;
  }

  if (field_parser.search(2, "ExtSolverLaunchScriptA", "ExtSolverLaunchScript")) {
    input_params.ext_solver_launch_script_A = field_parser.next(
        input_params.ext_solver_launch_script_A);
    std::cout << input_params.ext_solver_launch_script_A << std::endl;
  } else {
    input_params.ext_solver_launch_script_A = "";
    input_params.m_bParallelPossible = false;
    std::cout << " [CArl Parameters]WARNING! No external solver A launch script, parallel inversion not possible!" << std::endl;
  }
  

  if (field_parser.search(2, "ExtSolverLaunchScriptB", "ExtSolverLaunchScript")) {
    input_params.ext_solver_launch_script_B = field_parser.next(
        input_params.ext_solver_launch_script_B);
    std::cout << input_params.ext_solver_launch_script_B << std::endl;
  } else {
    input_params.ext_solver_launch_script_B = "";
    input_params.m_bParallelPossible = false;
    std::cout << " [CArl Parameters]WARNING! No external solver B launch script, parallel inversion not possible!" << std::endl;
  }
  



}

/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */