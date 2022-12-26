/*
 * carl_loop_dyn_input_parser.cpp
 *
 *  Created on: Nov 17, 2021
 *      Author: Chensheng Luo
 */

#include "carl_loop_dyn_input_parser.h"

namespace carl
{

bool is_multiple(double big,double small,double tolerance=1e-3){
  int divide=(int)(big/small);
  if(abs(big-small*divide)<(tolerance*small)){
    return true;
  }else{
    std::cout<<"WARNING: "<<big<<" is not multiple of "<<small<<"with a rest"<<abs(big-small*divide)<<std::endl;
    return false;
  }
}

void get_input_params(GetPot& field_parser,
    feti_loop_dyn_params& input_params){

  std::cout << " - Reading input files Begin ......" << std::endl;

  if (field_parser.search(1, "ClusterSchedulerType")) {
    std::string cluster_scheduler_type;
    cluster_scheduler_type = field_parser.next(cluster_scheduler_type);
    if(cluster_scheduler_type == "LOCAL")
    {
      std::cout << " !!! WARNING: " << std::endl;
      std::cout << "        The LOCAL 'scheduler' type is only intended for small and fast test cases" << std::endl;
      std::cout << "     on computers without a job scheduler (PBS, SLURM). You will have to launch  " << std::endl;
      std::cout << "     each script MANUALLY!!! " << std::endl;
      input_params.scheduler = carl::ClusterSchedulerType::LOCAL;
      input_params.script_filename = "";
    }
    else if(cluster_scheduler_type == "PBS") {
      input_params.scheduler = carl::ClusterSchedulerType::PBS;
      std::cout << "- Scheduler: PBS" << std::endl;
      if (field_parser.search(1, "ScriptFile")) {
        input_params.script_filename = field_parser.next(input_params.script_filename);
        std::cout << "- ScriptFile: "<< input_params.script_filename << std::endl;
      } else {
        homemade_error_msg("[CArl Parameters]ERROR! Missing the script file (needed for the PBS scheduler)!");
      }
    }
    else if(cluster_scheduler_type == "SLURM") {
      input_params.scheduler = carl::ClusterSchedulerType::SLURM;
      std::cout << "- Scheduler: SLURM" << std::endl;
      if (field_parser.search(1, "ScriptFile")) {
        input_params.script_filename = field_parser.next(input_params.script_filename);
        std::cout << "- ScriptFile: "<< input_params.script_filename << std::endl;
      } else {
        homemade_error_msg("[CArl Parameters]ERROR! Missing the script file (needed for the SLURM scheduler)!");
      }
    }
    else
      homemade_error_msg("[CArl Parameters]ERROR! Invalid scheduler type!");
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing the scheduler type!");
  }

  if (field_parser.search(1, "ScratchFolderPath")) {
    input_params.scratch_folder_path = field_parser.next(
        input_params.scratch_folder_path);
    std::cout << "- ScratchFolderPath: "<< input_params.scratch_folder_path << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing scratch folder path command line!");
  }

  if (field_parser.search(1, "ResultFolderPath")) {
    input_params.result_folder_path = field_parser.next(
        input_params.result_folder_path);
    std::cout << "- ResultFolderPath: "<< input_params.result_folder_path << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing result folder path command line!");
  }

  if (field_parser.search(1, "TildeMatrixA")) {
    input_params.matrix_A.mass_tilde = field_parser.next(
        input_params.matrix_A.mass_tilde);
    std::cout << "- MassTildeMatrixA: "<< input_params.matrix_A.mass_tilde << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing mass tilde matrix A command line!");
  }

  if (field_parser.search(1, "StiffnessMatrixA")) {
    input_params.matrix_A.stiffness = field_parser.next(
        input_params.matrix_A.stiffness);
    std::cout << "- StiffnessMatrixA: "<< input_params.matrix_A.stiffness << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing stiffness matrix A command line!");
  }

  if (field_parser.search(1, "DampingMatrixA")) {
    input_params.matrix_A.damping = field_parser.next(
        input_params.matrix_A.damping);
    std::cout << "- DampingMatrixA: "<< input_params.matrix_A.damping << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing damping matrix A command line!");
  }

    if (field_parser.search(1, "TildeMatrixB")) {
    input_params.matrix_B.mass_tilde = field_parser.next(
        input_params.matrix_B.mass_tilde);
    std::cout << "- MassTildeMatrixB: "<< input_params.matrix_B.mass_tilde << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing mass tilde matrix B command line!");
  }

  if (field_parser.search(1, "StiffnessMatrixB")) {
    input_params.matrix_B.stiffness = field_parser.next(
        input_params.matrix_B.stiffness);
    std::cout << "- StiffnessMatrixB: " << input_params.matrix_B.stiffness << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing stiffness matrix B command line!");
  }

  if (field_parser.search(1, "DampingMatrixB")) {
    input_params.matrix_B.damping = field_parser.next(
        input_params.matrix_B.damping);
    std::cout << "- DampingMatrixB: " << input_params.matrix_B.damping << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing damping matrix B command line!");
  }

  if (field_parser.search(1, "CouplingMatricesFolder")) {
    input_params.coupling_folder_path = field_parser.next(
        input_params.coupling_folder_path);
    std::cout << "- CouplingMatricesFolder: " << input_params.coupling_folder_path << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing the coupling matrices path!");
  }

  input_params.matrix_A.coupling = input_params.coupling_folder_path + "/coupling_matrix_macro.petscmat";
  input_params.matrix_B.coupling = input_params.coupling_folder_path + "/coupling_matrix_micro.petscmat";

  if (field_parser.search(1, "InterpolationMatrix")) {
    input_params.interpolation_matrix = field_parser.next(
        input_params.interpolation_matrix);
    std::cout << "- InterpolationMatrix: " << input_params.interpolation_matrix << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing interpolation matrix command line!");
  }



  if (field_parser.search(1, "InitialDispA")) {
    input_params.initial_A.disp = field_parser.next(
        input_params.initial_A.disp);
    std::cout << "- InitialDispA: "<< input_params.initial_A.disp << std::endl;
  } else {
    std::cout << "- InitialDispA: void [DEFAULT]" << std::endl;
    input_params.initial_A.disp = "";
  }
  

  if (field_parser.search(1, "InitialDispB")) {
    input_params.initial_B.disp = field_parser.next(
        input_params.initial_B.disp);
    std::cout << "- InitialDispB: "<< input_params.initial_B.disp << std::endl;
  } else {
    std::cout << "- InitialDispB: void [DEFAULT]" << std::endl;
    input_params.initial_B.disp = "";
  }
  

  if (field_parser.search(1, "InitialSpeedA")) {
    input_params.initial_A.speed = field_parser.next(
        input_params.initial_A.speed);
    std::cout << "- InitialSpeedA: "<< input_params.initial_A.speed << std::endl;
  } else {
    std::cout << "- InitialSpeedA: void [DEFAULT]" << std::endl;
    input_params.initial_A.speed = "";
  }
  

  if (field_parser.search(1, "InitialSpeedB")) {
    input_params.initial_B.speed = field_parser.next(
        input_params.initial_B.speed);
    std::cout << "- InitialSpeedB: "<< input_params.initial_B.speed << std::endl;
  } else {
    std::cout << "- InitialSpeedB: void [DEFAULT]" << std::endl;
    input_params.initial_B.speed = "";
  }
  



  if (field_parser.search(1, "ForcePrepareMethod")) {
    std::string prepare_method;
    prepare_method = field_parser.next(prepare_method);
    if(prepare_method == "ModalSinus"){
      input_params.force_prepare_method = carl::ForcePrepareMethod::MODAL_SINUS;
      std::cout << "- ForcePrepareMethod: ModalSinus"<< std::endl;
    }else if(prepare_method == "ModalConstant"){
      input_params.force_prepare_method = carl::ForcePrepareMethod::MODAL_CONSTANT;
      std::cout << "- ForcePrepareMethod: ModalConstant"<< std::endl;
    }else if(prepare_method == "ModalLinear"){
      input_params.force_prepare_method = carl::ForcePrepareMethod::MODAL_LINEAR;
      std::cout << "- ForcePrepareMethod: ModalLinear"<< std::endl;
    }else if(prepare_method == "ModalProduct"){
      input_params.force_prepare_method = carl::ForcePrepareMethod::MODAL_PRODUCT;
      std::cout << "- ForcePrepareMethod: ModalProduct"<< std::endl;
    }else{
      homemade_error_msg("[CArl Parameters]ERROR! Invalid force prepare method!");
    }
  }else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing force prepare method!");
  }

  if (field_parser.search(1, "ForcePrepareParams")) {
    input_params.force_prepare_params = field_parser.next(
        input_params.force_prepare_params);
    std::cout << "- ForcePrepareParams: "<< input_params.force_prepare_params << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing force prepare params input command line!");
  }

  if (field_parser.search(2, "ExtSolverInputA", "ExtSolverInput")) {
    input_params.ext_solver_A_input = field_parser.next(
        input_params.ext_solver_A_input);
    std::cout << "- ExtSolverInputA: "<< input_params.ext_solver_A_input << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing ext solver A input command line!");
  }

  if (field_parser.search(2, "ExtSolverInputB", "ExtSolverInput")) {
    input_params.ext_solver_B_input = field_parser.next(
        input_params.ext_solver_B_input);
    std::cout << "- ExtSolverInputB: "<< input_params.ext_solver_B_input << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing ext solver B input command line!");
  }

  if (field_parser.search(2, "ExtSolverInputInterpolation", "ExtSolverInput")) {
    input_params.ext_solver_general_input = field_parser.next(
        input_params.ext_solver_general_input);
    std::cout << "- ExtSolverInputInterpolation: "<< input_params.ext_solver_general_input << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing ext solver interpolation input command line!");
  }

  if (field_parser.search(2, "ExtSolverLaunchScriptA", "ExtSolverLaunchScript")) {
    input_params.ext_solver_launch_script_A = field_parser.next(
        input_params.ext_solver_launch_script_A);
    std::cout << "- ExtSolverLaunchScriptA: "<< input_params.ext_solver_launch_script_A << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing ext solver launch script command line for A!");
  }
  

  if (field_parser.search(2, "ExtSolverLaunchScriptB", "ExtSolverLaunchScript")) {
    input_params.ext_solver_launch_script_B = field_parser.next(
        input_params.ext_solver_launch_script_B);
    std::cout << "- ExtSolverLaunchScriptB: "<< input_params.ext_solver_launch_script_B << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing ext solver launch script command line for B!");
  }
  
  //if (input_params.dyn_solver==carl::DynamicSolver::DI){
    if (field_parser.search(2, "ExtSolverLaunchScriptInterpolation", "ExtSolverLaunchScript")) {
      input_params.ext_solver_launch_script_general = field_parser.next(
          input_params.ext_solver_launch_script_general);
      std::cout << "- ExtSolverLaunchScriptInterpolation: "<< input_params.ext_solver_launch_script_general << std::endl;
    } else{
      homemade_error_msg("[CArl Parameters]ERROR! Missing ext solver launch script command line!");
    }
  //}



  if (field_parser.search(1, "GeneralEntryFilePath")) {
    input_params.general_entry_file_path = field_parser.next(
        input_params.general_entry_file_path);
    std::cout << "- GeneralEntryFilePath: "<< input_params.general_entry_file_path << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing scratch folder path command line!");
  }


  if (field_parser.search(2, "NewmarkParametersA","NewmarkParameters")){
        std::string filename;
        filename = field_parser.next(filename);
        std::cout << "- NewmarkParametersFilePathA: "<< filename << std::endl;
        GetPot newmark_parser;
        newmark_parser.parse_input_file(filename, "#", "\n", " \t\n");
        carl::get_newmark_params(newmark_parser, input_params.newmark_A);
      } else{
        homemade_error_msg("[CArl Parameters]ERROR! Missing A Newmark parameters file!");
      }

    if (field_parser.search(2, "NewmarkParametersB","NewmarkParameters")){
        std::string filename;
        filename = field_parser.next(filename);
        std::cout << "- NewmarkParametersFilePathB: "<< filename << std::endl;
        GetPot newmark_parser;
        newmark_parser.parse_input_file(filename, "#", "\n", " \t\n");
        carl::get_newmark_params(newmark_parser, input_params.newmark_B);
      } else{
        homemade_error_msg("[CArl Parameters]ERROR! Missing B Newmark parameters file!");
      }

  if (field_parser.search(1, "SimulationDuration")) {
    input_params.simulation_duration = field_parser.next(input_params.simulation_duration);
    std::cout << "- SimulationDuration: "<< input_params.simulation_duration << std::endl;
  } else {
    input_params.simulation_duration = 1;
    std::cout << "- SimulationDuration: "<< input_params.simulation_duration  << "[No entry, default value is taken]"<< std::endl;
  }

  input_params.inner_loop_times=std::max((int)(input_params.newmark_A.deltat/input_params.newmark_B.deltat),1);

  if (is_multiple(input_params.newmark_A.deltat,input_params.newmark_B.deltat)==false){
    homemade_error_msg("[Newmark Parameter]ERROR! Delta t for A is not the integer multiple of Delta t for B, please reset and restart the calculation from beginning!");
  }
  
  input_params.outer_loop_times=std::max((int)(input_params.simulation_duration/input_params.newmark_A.deltat),1);

  if (is_multiple(input_params.simulation_duration,input_params.newmark_A.deltat)==false){
    input_params.simulation_duration = input_params.outer_loop_times*input_params.newmark_A.deltat;
    std::cout << "[Newmark Parameter]WARNING! Simulation Duration is not the integer multiple of Delta t for A, auto-corrected to "<< input_params.simulation_duration  << std::endl;
  }


  double ResultTimeA;
  if (field_parser.search(2, "ResultTimeA", "ResultTime")) {
    ResultTimeA = field_parser.next(ResultTimeA);
    std::cout << "- ResultTimeA: "<< ResultTimeA  << std::endl;
    input_params.result_times_A = std::max((int)(ResultTimeA/input_params.newmark_A.deltat),1);
    if(is_multiple(ResultTimeA,input_params.newmark_A.deltat)==false){
      std::cout << "[Newmark Parameter]WARNING! Result time for A is not the integer multiple of deltat A, auto-corrected to"<< input_params.result_times_A * input_params.newmark_A.deltat << std::endl;
    }
  }else {
    input_params.result_times_A = 1;
    std::cout << "- ResultTimeA: "<< input_params.newmark_A.deltat << "[No entry, default value is taken]" << std::endl;
    std::cout << "[Newmark Parameter]WARNING! Result time isn't set for system A,  default value is taken as:"<< input_params.newmark_A.deltat << std::endl;
  }

  double ResultTimeB;
  if (field_parser.search(2, "ResultTimeB","ResultTime")) {
    ResultTimeB = field_parser.next(ResultTimeB);
    std::cout << "- ResultTimeB: "<< ResultTimeB  << std::endl;
    input_params.result_times_B = std::max((int)(ResultTimeB/input_params.newmark_B.deltat),1);
    if(is_multiple(ResultTimeB,input_params.newmark_B.deltat)==false){
      std::cout << "[Newmark Parameter]WARNING! Result time for B is not the integer multiple of deltat B, auto-corrected to"<< input_params.result_times_B * input_params.newmark_B.deltat << std::endl;
    }
    
  } else {
    input_params.result_times_B = 1;
    std::cout << "- ResultTimeB: "<< input_params.newmark_B.deltat << "[No entry, default value is taken]" << std::endl;
    std::cout << "[Newmark Parameter]WARNING! Result time isn't set for system B,  default value is taken as:"<< input_params.newmark_B.deltat << std::endl;
  }

  
  // if (field_parser.search(1,"UseRigidBodyModesB"))
  // {
  //   input_params.bUseRigidBodyModes = true;

  //   if (field_parser.search(1, "RBVectorBase")) {
  //     input_params.RB_vectors_base = field_parser.next(
  //         input_params.RB_vectors_base);
  //     std::cout << input_params.RB_vectors_base << std::endl;
  //   } else {
  //     homemade_error_msg("Missing the system B's rigid body mode vectors!");
  //   }
    
  //   if (field_parser.search(1, "NbOfRBVectors")) {
  //     input_params.nb_of_rb_vectors = field_parser.next(
  //         input_params.nb_of_rb_vectors);
  //   } else {
  //     input_params.nb_of_rb_vectors = 6;
  //   }
  // }
  // else
  // {
  //   input_params.bUseRigidBodyModes = false;
  // }

  // // Set CG coupling solver convergence
  // input_params.CG_coupled_conv_abs = 1e-20;
  // input_params.CG_coupled_conv_rel = 1e-5;
  // input_params.CG_coupled_div = 1e5;
  // input_params.CG_coupled_conv_max = 1e4;
  // input_params.CG_coupled_conv_corr =1e-5;

  // if( field_parser.search(1,"CoupledConvAbs") )
  // {
  //   input_params.CG_coupled_conv_abs = field_parser.next(input_params.CG_coupled_conv_abs);
  // }
  // if( field_parser.search(1,"CoupledConvRel") )
  // {
  //   input_params.CG_coupled_conv_rel = field_parser.next(input_params.CG_coupled_conv_rel);
  // }
  // if( field_parser.search(1,"CoupledCorrConvRel") )
  // {
  //   input_params.CG_coupled_conv_corr = field_parser.next(input_params.CG_coupled_conv_corr);
  // }
  // if( field_parser.search(1,"CoupledDiv") )
  // {
  //   input_params.CG_coupled_div = field_parser.next(input_params.CG_coupled_div);
  // }
  // if( field_parser.search(1,"CoupledIterMax") )
  // {
  //   input_params.CG_coupled_conv_max = field_parser.next(input_params.CG_coupled_conv_max);
  // }
  
  // input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_OPERATOR;
  // if ( field_parser.search(1, "CGPreconditionerType") )
  // {
  //   std::string CG_precond_type_string = field_parser.next(CG_precond_type_string);
  //   if(CG_precond_type_string == "NONE")
  //     input_params.CG_precond_type = carl::BaseCGPrecondType::NO_PRECONDITIONER;
  //   else if(CG_precond_type_string == "Coupling_operator")
  //     input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_OPERATOR;
  //   else if(CG_precond_type_string == "Coupling_operator_jacobi")
  //     input_params.CG_precond_type = carl::BaseCGPrecondType::COUPLING_JACOBI;
  // }

  std::cout << " ...... Reading input files: FINISH!  " << std::endl;

  };
}