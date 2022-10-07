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



}

/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */