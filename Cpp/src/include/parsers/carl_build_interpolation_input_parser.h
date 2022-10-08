/*
 * carl_build_interpolation_input_parser.h
 *
 *  Created on: Sep 26, 2021
 *      Author: Chensheng Luo
 */

#ifndef CARL_BUILD_INTERPOLATION_INPUT_PARSER_H_
#define CARL_BUILD_INTERPOLATION_INPUT_PARSER_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"
#include "newmark_param_parser.h"

/// Structure containing the parameters for the interpolation matrix calculation.
struct carl_build_interpolation_params{
  std::string path_tilde_matrix_A;          ///< Path to mass tilde matrix for A
  std::string path_macro_coupling_matrix;   ///< Path to coupling matrix for A
  std::string path_inv_matrix_A;            ///< Path to store the invert mass matrix for A

  std::string path_tilde_matrix_B;          ///< Path to mass tilde matrix for B
  std::string path_micro_coupling_matrix;   ///< Path to coupling matrix for B
  std::string path_inv_matrix_B;            ///< Path to store the invert mass matrix for B

  std::string output_base;                  ///< Output base of results
  
  carl::NewmarkParams newmark_A;            ///< Newmark parameter for A
  carl::NewmarkParams newmark_B;            ///< Newmark parameter for B
  
  double invert_threshold;                  ///< Acceptable threshold for \f$ || M^k (M^k)^{-1}-I||_{\inf} \f$, default to be inf(check manually)

};

/** \brief **DYN-DI**Structure containing the parameters for the interpolation matrix calculation
 *  
 
Parameters:

  + **Matrix input/output**:
    - `mTildeA` : Sum of Mass matrix and Newmark matrix for part A, (i.e. \f$ \tilde{M}^A \f$)obtained by external solver
    - `mTildeB` : Sum of Mass matrix and Newmark matrix for part B, (i.e. \f$ \tilde{M}^B \f$)obtained by external solver
    - `couplingMatrixA` : Coupling matrix for part A, obtained in coupling step
    - `couplingMatrixB` : Coupling matrix for part B, obtained in coupling step
    - `OutputFolder` : Folder where the result is stored

  + **Newmark parameters**:
    - `NewmarkParameters` *or* (`NewmarkParametersA` and `NewmarkParametersB`) : Path to the Newmark parameter of two model, see get_newmark_params(GetPot& field_parser,NewmarkParams& newmark) for its redaction format.

  + (*OPTIONAL WITH DEFAULT*)**Matrix ouput**:
    - `InverseNameA` : Path in which the inverse matrix of Mass for part A will be stored. *Default*:`inv_dense_A`
    - `InverseNameB` : Path in which the inverse matrix of Mass for part B will be stored. *Default*:`inv_dense_B`

  + (*OPTIONAL*)Invert examination threshold
    - `InvertThreshold` : Acceptable threshold for \f$ || M^k (M^k)^{-1}-I||_{\inf} \f$. If no value entered, please check manually result_test_inversion.txt file in output folder.

 * 
 */

void get_input_params(GetPot& field_parser, 
	carl_build_interpolation_params& input_params);



#endif /* CARL_BUILD_INTERPOLATION_INPUT_PARSER_H_*/
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */