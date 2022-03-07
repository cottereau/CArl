/*
 * libmesh_build_interpolation_matrix_input_parser.cpp
 *
 *  Created on: Sep 26, 2021
 *      Author: Chensheng Luo
 */

#include "libmesh_build_interpolation_matrix_input_parser.h"

void get_input_params(GetPot& field_parser,libmesh_build_interpolation_params& input_params) 
{
  //Set mass tilde matrix for A
  if (field_parser.search(1, "mTildeA")) {
    input_params.path_tilde_matrix_A = field_parser.next(input_params.path_tilde_matrix_A);
  } else {
    homemade_error_msg("[CArl Parameters]Missing the mass tilde file!");
  }

  //Set mass tilde matrix for B
  if (field_parser.search(1, "mTildeB")) {
    input_params.path_tilde_matrix_B = field_parser.next(input_params.path_tilde_matrix_B);
  } else {
    homemade_error_msg("[CArl Parameters]Missing the mass tilde file!");
  }

    //Set coupling matrix for A
  if (field_parser.search(1, "couplingMatrixA")) {
    input_params.path_macro_coupling_matrix = field_parser.next(input_params.path_macro_coupling_matrix);
  } else {
    homemade_error_msg("[CArl Parameters]Missing the coupling matrix for A file!");
  }

    //Set coupling matrix for B
  if (field_parser.search(1, "couplingMatrixB")) {
    input_params.path_micro_coupling_matrix = field_parser.next(input_params.path_micro_coupling_matrix);
  } else {
    homemade_error_msg("[CArl Parameters]Missing the coupling matrix for B file!");
  }

  //Set output folder
  if (field_parser.search(1, "OutputFolder")) {
    input_params.output_base = field_parser.next(input_params.output_base);
  } else {
    homemade_error_msg("[CArl Parameters]Missing output base!");
  }

  if (field_parser.search(1, "AlphaA")) {
    input_params.alpha_A = field_parser.next(input_params.alpha_A);
  } else {
    input_params.alpha_A = 0.25;
  }

  if (field_parser.search(1, "AlphaB")) {
    input_params.alpha_B = field_parser.next(input_params.alpha_B);
  } else {
    input_params.alpha_B = 0.25;
  }

  if (field_parser.search(1, "InverseNameA")) {
    input_params.path_inv_matrix_A = field_parser.next(input_params.path_inv_matrix_B);
  } else {
    input_params.path_inv_matrix_A = 0.25;
  }

  if (field_parser.search(1, "InverseNameB")) {
    input_params.path_inv_matrix_B = field_parser.next(input_params.path_inv_matrix_B);
  } else {
    input_params.path_inv_matrix_B = 0.25;
  }

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