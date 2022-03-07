/*
 * libmesh_build_interpolation_matrix_input_parser.h
 *
 *  Created on: Sep 26, 2021
 *      Author: Chensheng Luo
 */

#ifndef LIBMESH_BUILD_INTERPOLATION_MATRIX_INPUT_PARSER_H_
#define LIBMESH_BUILD_INTERPOLATION_MATRIX_INPUT_PARSER_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"


struct libmesh_build_interpolation_params{
  std::string path_tilde_matrix_A;          ///< Path to mass tilde matrix for A
  std::string path_macro_coupling_matrix;   ///< coupling matrix for A
  std::string path_inv_matrix_A;            ///< invert mass matrix for A

  std::string path_tilde_matrix_B;          ///< Path to mass tilde matrix for B
  std::string path_micro_coupling_matrix;   ///< coupling matrix for B
  std::string path_inv_matrix_B;            ///< invert mass matrix for B
  
  double alpha_A; //step time of calculation
  double alpha_B; 
  std::string output_base;
};


void get_input_params(GetPot& field_parser, 
	libmesh_build_interpolation_params& input_params);



#endif /* LIBMESH_BUILD_INTERPOLATION_MATRIX_INPUT_PARSER_H_ */
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */