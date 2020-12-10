/*
 * carl_feti_setup_init_input_parser.h
 *
 *  Created on: Feb 14, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "carl_feti_setup_init_input_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
    feti_setup_init_params& input_params) {

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

  if (field_parser.search(1, "ExtSolverA")) {
    input_params.ext_solver_BIG = field_parser.next(
        input_params.ext_solver_BIG);
    std::cout << input_params.ext_solver_BIG << std::endl;
  } else {
    homemade_error_msg("Missing the external solver A command line!");
  }

  if (field_parser.search(1, "ExtSolverB")) {
    input_params.ext_solver_micro = field_parser.next(
        input_params.ext_solver_micro);
    std::cout << input_params.ext_solver_micro << std::endl;
  } else {
    homemade_error_msg("Missing the external solver B command line!");
  }

  std::string ext_solver_type;
  if (field_parser.search(1, "ExtSolverAType")) {
    ext_solver_type = field_parser.next(
        ext_solver_type);
    std::cout << ext_solver_type << std::endl;
    if(ext_solver_type == "LIBMESH_LINEAR")
    {
      input_params.ext_solver_BIG_type = carl::ExtSolverType::LIBMESH_LINEAR;
    } else if(ext_solver_type == "LIBMESH_DYNAMIC")
      {
        input_params.ext_solver_BIG_type = carl::ExtSolverType::LIBMESH_DYNAMIC;
      } else if(ext_solver_type == "LIBMESH_NONLINGEOM")
      {
        input_params.ext_solver_BIG_type = carl::ExtSolverType::LIBMESH_NONLINGEOM;
    } else { 
      homemade_error_msg("Invalid external solver A type!");
    }
  } else {
    homemade_error_msg("Missing the external solver A type!");
  }

  if (field_parser.search(1, "ExtSolverBType")) {
    ext_solver_type = field_parser.next(
        ext_solver_type);
    std::cout << ext_solver_type << std::endl;
    if(ext_solver_type == "LIBMESH_LINEAR")
    {
      input_params.ext_solver_BIG_type = carl::ExtSolverType::LIBMESH_LINEAR;
    } else if(ext_solver_type == "LIBMESH_DYNAMIC")
      {
        input_params.ext_solver_BIG_type = carl::ExtSolverType::LIBMESH_DYNAMIC;
      } else if(ext_solver_type == "LIBMESH_NONLINGEOM")
      {
        input_params.ext_solver_BIG_type = carl::ExtSolverType::LIBMESH_NONLINGEOM;
    } else { 
      homemade_error_msg("Invalid external solver B type!");
    }
  } else {
    homemade_error_msg("Missing the external solver B type!");
  }

  if (field_parser.search(1, "ExtSolverAInput")) {
    input_params.ext_solver_BIG_input = field_parser.next(
        input_params.ext_solver_BIG_input);
    std::cout << input_params.ext_solver_BIG_input << std::endl;
  } else {
    homemade_error_msg("Missing the input file for the external solver A!");
  }

  if (field_parser.search(1, "ExtSolverBInput")) {
    input_params.ext_solver_micro_input = field_parser.next(
        input_params.ext_solver_micro_input);
    std::cout << input_params.ext_solver_micro_input << std::endl;
  } else {
    homemade_error_msg("Missing the input file for the external solver B!");
  }

  if (field_parser.search(1, "ScratchFolderPath")) {
    input_params.scratch_folder_path = field_parser.next(
        input_params.scratch_folder_path);
    std::cout << input_params.scratch_folder_path << std::endl;
  } else {
    homemade_error_msg("Missing the external scratch folder path!");
  }

  if (field_parser.search(1, "CouplingMatricesFolder")) {
    input_params.coupling_folder_path = field_parser.next(
        input_params.coupling_folder_path);
    std::cout << input_params.coupling_folder_path << std::endl;
  } else {
    homemade_error_msg("Missing the coupling matrices path!");
  }

  if (field_parser.search(1,"UseRigidBodyModesB"))
  {
    input_params.bUseRigidBodyModes = true;
    if (field_parser.search(1, "ExtForceSystemB")) {
      input_params.force_micro_path = field_parser.next(
          input_params.force_micro_path);
      std::cout << input_params.force_micro_path << std::endl;
    } else {
      homemade_error_msg("Missing the system B's external force file path!");
    }

    if (field_parser.search(1, "RBVectorBase")) {
      input_params.RB_vectors_base = field_parser.next(
          input_params.RB_vectors_base);
      std::cout << input_params.RB_vectors_base << std::endl;
    } else {
      homemade_error_msg("Missing the system B's rigid body mode vectors!");
    }
    
    if (field_parser.search(1, "NbOfRBVectors")) {
      input_params.nb_of_rb_vectors = field_parser.next(
          input_params.nb_of_rb_vectors);
    } else {
      input_params.nb_of_rb_vectors = 6;
    }
  }
  else
  {
    input_params.bUseRigidBodyModes = false;
  }

  // Set CG coupling solver convergence
  input_params.CG_coupled_conv_abs = 1e-20;
  input_params.CG_coupled_conv_rel = 1e-5;
  input_params.CG_coupled_div = 1e5;
  input_params.CG_coupled_conv_max = 1e4;
  input_params.CG_coupled_conv_corr =1e-5;

  if( field_parser.search(1,"CoupledConvAbs") )
  {
    input_params.CG_coupled_conv_abs = field_parser.next(input_params.CG_coupled_conv_abs);
  }
  if( field_parser.search(1,"CoupledConvRel") )
  {
    input_params.CG_coupled_conv_rel = field_parser.next(input_params.CG_coupled_conv_rel);
  }
  if( field_parser.search(1,"CoupledCorrConvRel") )
  {
    input_params.CG_coupled_conv_corr = field_parser.next(input_params.CG_coupled_conv_corr);
  }
  if( field_parser.search(1,"CoupledDiv") )
  {
    input_params.CG_coupled_div = field_parser.next(input_params.CG_coupled_div);
  }
  if( field_parser.search(1,"CoupledIterMax") )
  {
    input_params.CG_coupled_conv_max = field_parser.next(input_params.CG_coupled_conv_max);
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

  if (field_parser.search(1,"OutputFolder")) {
    input_params.output_folder = field_parser.next(
        input_params.output_folder);
  } else {
    homemade_error_msg("Missing the output filename base!");
  }
};

};
