/*
 * libmesh_assemble_system_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "libmesh_assemble_system_input_parser.h"

void get_input_params(GetPot& field_parser,libmesh_assemble_input_params& input_params) 
{
  // Set mesh files
  if (field_parser.search(1, "Mesh")) {
    input_params.mesh_file = field_parser.next(input_params.mesh_file);
  } else {
    homemade_error_msg("[CArl Parameters]Missing the system mesh file!");
  }

  // Set constant parameters
  if ( field_parser.search(1, "PhysicalParameters") )
  {
    input_params.physical_params_file = field_parser.next(input_params.physical_params_file);
  }
  else
  {
    homemade_error_msg("[CArl Parameters]Missing the physical parameters file!");
  }

  // Set weight function
  std::string sys_type;
  if ( field_parser.search(1, "SystemType") )
  {
    sys_type = field_parser.next(sys_type);
    if(sys_type == "Macro" || sys_type == "MACRO" || sys_type == "macro")
      input_params.system_type = WeightFunctionSystemType::MACRO;
    else if(sys_type == "Micro" || sys_type == "MICRO" || sys_type == "micro")
      input_params.system_type = WeightFunctionSystemType::MICRO;
    else if(sys_type == "NoWeight" || sys_type == "NOWEIGHT" || sys_type == "noweight")
    {
      input_params.system_type = WeightFunctionSystemType::NO_WEIGHT;
      std::cout << " >> [CArl Parameters]Warning: Will not use the weight parameters!" << std::endl;
    }
    else
      homemade_error_msg("[CArl Parameters]Invalid system type (must be either Macro, Micro or NoWeight)!");
  } else {
    homemade_error_msg("[CArl Parameters]Missing the system type (must be either Macro, Micro or NoWeight)!");
  }

  if ( field_parser.search(1, "MeshWeight") )
  {
    input_params.mesh_weight_file = field_parser.next(input_params.mesh_weight_file);
  } else {
    homemade_error_msg("[CArl Parameters]Missing the weight mesh file!");
  }

  if( field_parser.search(1, "WeightIndexes") )
  {
    input_params.weight_domain_idx_file = field_parser.next(input_params.weight_domain_idx_file);
  } else {
    homemade_error_msg("[CArl Parameters]Missing the weight value file!");
  }

  // Output
  if (field_parser.search(2, "--output", "OutputBase"))
  {
    input_params.output_base = field_parser.next(
      input_params.output_base);
  } else {
    input_params.output_base = "test_system";
  }

  if (field_parser.search(1, "ExportRBVectors")) {
    input_params.bCalculateRBVectors = true;
  } else {
    input_params.bCalculateRBVectors = false;
  }
//}

  if (field_parser.search(1,"Dynamic")) {
    input_params.dynamic_analysis = true;
  } else {
    input_params.dynamic_analysis = false;
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
    input_params.deltatB = 0.25;
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
//  if (field_parser.search(1, "transient")) {
//    input_params.transient = true;
//  } else {
//    input_params.transient = false;
//  }
//
//  if (field_parser.search(1, "n_timesteps")) {
//    input_params.n_timesteps = field_parser.next(
//        input_params.n_timesteps);
//  } else {
//    input_params.n_timesteps = 1;
//  }

//  if (field_parser.search(1, "write_interval")) {
//    input_params.write_interval = field_parser.next(
//        input_params.write_interval);
//  } else {
//    input_params.write_interval = 1;
//  }

  // if (field_parser.search(1, "max_nonlinear_iterations")) {
  //   input_params.max_nonlinear_iterations  = field_parser.next(
  //       input_params.max_nonlinear_iterations );
  // } else {
  //   input_params.max_nonlinear_iterations = 15;
  // }

  // if (field_parser.search(1, "max_linear_iterations")) {
  //   input_params.max_linear_iterations  = field_parser.next(
  //       input_params.max_linear_iterations );
  // }else {
  //   input_params.max_linear_iterations = 50000;
  // }

  // if (field_parser.search(1, "initial_linear_tolerance")) {
  //   input_params.initial_linear_tolerance  = field_parser.next(
  //       input_params.initial_linear_tolerance );
  // }else {
  //   input_params.initial_linear_tolerance = 1.e-3;
  // }

  // if (field_parser.search(1, "absolute_residual_tolerance")) {
  //   input_params.absolute_residual_tolerance  = field_parser.next(
  //       input_params.absolute_residual_tolerance );
  // }else {
  //   input_params.absolute_residual_tolerance = 0.0;
  // }
  // if (field_parser.search(1, "solver_quite")) {
  //   input_params.solver_quiet = true;
  // } else {
  //   input_params.solver_quiet = false;
  // }

  // if (field_parser.search(1, "relative_step_tolerance")) {
  //   input_params.relative_step_tolerance = field_parser.next(
  //       input_params.relative_step_tolerance);
  // } else {
  //   input_params.relative_step_tolerance = 1.e-3;
  // }

  // if (field_parser.search(1, "relative_residual_tolerance")) {
  //   input_params.relative_residual_tolerance = field_parser.next(
  //       input_params.relative_residual_tolerance);
  // } else {
  //   input_params.relative_residual_tolerance = 0.0;
  // }
}

// void get_material_params(GetPot& field_parser, 
//  material_params& input_params)
// {
//  if (field_parser.search(1, "rho")) {
//    input_params.rho = field_parser.next(
//        input_params.rho);
//  } else {
//    input_params.rho = 1000;
//    homemade_error_msg("[WARNING][Material Parameters]Missing the paramater rho");
//  }

//  if (field_parser.search(1, "poisson_ratio")) {
//    input_params.poisson_ratio = field_parser.next(
//        input_params.poisson_ratio);
//  } else {
//    input_params.poisson_ratio = 0.2;
//    homemade_error_msg("[WARNING][Material Parameters]Missing the paramater poisson_ratio");
//  }

//  if (field_parser.search(1, "young_modulus")) {
//    input_params.young_modulus = field_parser.next(
//        input_params.young_modulus);
//  } else {
//    input_params.young_modulus = 3e6;
//    homemade_error_msg("[WARNING][Material Parameters]Missing the paramater young_modulus");
//  }
// }


// void get_file_params(GetPot command_line,
//          libmesh_assemble_input_params& carl_params,  
//          material_params& input_material_params)
// {
//  // File parser
//   GetPot field_parser;
//   // If there is an input file, parse it to get the parameters. Else, parse the command line
//   std::string input_filename;
//   if (command_line.search(1, "-i")) 
//   {
//     input_filename = command_line.next(input_filename);
//     field_parser.parse_input_file(input_filename, "#", "\n", "\t\n");
//     get_input_params(field_parser,carl_params);
//   } 
//   else 
//   {
//     field_parser = command_line;
//   }

//  //If there is a file in which dynamic parameter are passed
//  if (command_line.search(1, "d"))
//  {
//    input_filename = command_line.next(input_filename);
//    field_parser.parse_input_file(input_filename, "#", "\n", "\t\n");
//    get_dynamic_params(field_parser,input_dynamic_params);
//  } 
//    else 
//  {
//      printf("No dynamic file passed");
//  }
//
//  //If there is a file in which material parameter are passed
//  if (command_line.search(2, "-m"))
//  {
//      input_filename = command_line.next(input_filename);
//      field_parser.parse_input_file(input_filename, "#", "\n", "\t\n");
//      get_material_params(field_parser, input_material_params);
//  } 
//    else 
//  {
//      printf("No material file passed");
//  }
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */
