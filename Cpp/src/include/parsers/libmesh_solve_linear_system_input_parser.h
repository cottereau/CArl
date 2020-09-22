/*
 * libmesh_solve_linear_system_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef LIBMESH_SOLVE_LINEAR_SYSTEM_INPUT_PARSER_H_
#define LIBMESH_SOLVE_LINEAR_SYSTEM_INPUT_PARSER_H_

#include "carl_headers.h"
#include "ext_solver_libmesh_enums.h"

namespace carl
{

struct libmesh_solve_linear_system_input_params {
  std::string sys_matrix_file;  ///< Path to the system matrix file.
  std::string sys_rhs_vec_file; ///< Path to the system RHS vector file.
  std::string output_base;  ///< Output filename base.

  double sys_eps;   ///< Relative convergence parameter.
  int sys_iter_div; ///< Maximum number of iterations.

  std::string path_to_rb_vectors; ///< Path to a folder containing the rigid body mode vectors.
  int nb_of_rb_vectors; ///< Number of RB mode vectors.
  bool bUseRBVectors;   ///< Use rigid body mode vectors? *Default*: `false`.
  double deltatA;       ///< time_step A
  double betaA;         ///< beta Newmark coefficient for A   
  double gammaA;        ///< gamma Newmark coefficient for A   
  double deltatB;       ///< time_step A
  double betaB;         ///< beta Newmark coefficient for A   
  double gammaB;        ///< gamma Newmark coefficient for A   
//  bool transient;
  unsigned int n_timesteps; // number of time step
//  unsigned int write_interval;
//  bool solver_quiet;
//  double relative_step_tolerance;
//  double relative_residual_tolerance;
//  unsigned int max_nonlinear_iterations;
//  unsigned int max_linear_iterations;
//  double initial_linear_tolerance;
//  double absolute_residual_tolerance;
//
};

/** \brief Parser function for the coupled solver test programs.
 *  
 *  Required parameters:
 *    - `SysMatrix` : path to the system matrix file.
 *    - `SysRHSVector` : path to the system RHS vector file.
 *    - `OutputBase` : output filename base.
 *
 *  Optional parameters:
 *    - `SysEps` : relative convergence parameter. *Default*: 1e-5.
 *    - `SysIterDiv` : maximum number of iterations. *Default*: 1e3.
 *    - `RBVectorBase` : filename base to the rigid body mode vectors. The program expects that these vectors will be named as `[RBVectorBase]_rb_vector_XYZ_n_[NbOfRBVectors].petscvec`, where `XYZ` is an integer index going from `0` to `NbOfRBVectors - 1`. Setting this parameter sets `bUseRBVectors` to `true` - else, it is set to `false`.
 *    - `NbOfRBVectors` : number of RB mode vectors. *Default*: 6.
 */
void get_input_params(GetPot& field_parser,
    libmesh_solve_linear_system_input_params& input_params);

/// Function used to generate a solver input file from "input_params"
void print_input_params(const std::string& output_filename,
    libmesh_solve_linear_system_input_params& input_params);
}


#endif /* LIBMESH_SOLVE_LINEAR_SYSTEM_INPUT_PARSER_H_ */
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */
