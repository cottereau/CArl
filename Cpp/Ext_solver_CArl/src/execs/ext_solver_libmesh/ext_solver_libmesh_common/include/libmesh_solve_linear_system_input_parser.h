/*
 * libmesh_solve_linear_system_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef LIBMESH_SOLVE_LINEAR_SYSTEM_INPUT_PARSER_H_
#define LIBMESH_SOLVE_LINEAR_SYSTEM_INPUT_PARSER_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"

struct libmesh_solve_linear_system_input_params {
	std::string sys_matrix_file;	///< Path to the system matrix file.
	std::string sys_rhs_vec_file;	///< Path to the system RHS vector file.
	std::string output_base; 	///< Output filename base.

	double sys_eps;		///< Relative convergence parameter.
	int sys_iter_div;	///< Maximum number of iterations.

	std::string path_to_rb_vectors;	///< Path to a folder containing the rigid body mode vectors.
	int nb_of_rb_vectors;	///< Number of RB mode vectors.
	bool bUseRBVectors;		///< Use rigid body mode vectors? *Default*: `false`.
};

/**	\brief Parser function for the coupled solver test programs.
 *	
 *	Required parameters:
 *	  - `SysMatrix` : path to the system matrix file.
 *    - `SysRHSVector` : path to the system RHS vector file.
 *    - `OutputBase` : output filename base.
 *
 *  Optional parameters:
 *	  - `SysEps` : relative convergence parameter. *Default*: 1e-5.
 *	  - `SysIterDiv` : maximum number of iterations. *Default*: 1e3.
 *    - `RBVectorBase` : filename base to the rigid body mode vectors. The program expects that these vectors will be named as `[RBVectorBase]_rb_vector_XYZ_n_[NbOfRBVectors].petscvec`, where `XYZ` is an integer index going from `0` to `NbOfRBVectors - 1`. Setting this parameter sets `bUseRBVectors` to `true` - else, it is set to `false`.
 *    - `NbOfRBVectors` : number of RB mode vectors. *Default*: 6.
 */
void get_input_params(GetPot& field_parser,
		libmesh_solve_linear_system_input_params& input_params) {

	if (field_parser.search(1, "SysMatrix")) {
		input_params.sys_matrix_file = field_parser.next(
				input_params.sys_matrix_file);
	} else {
		homemade_error_msg("Missing the system matrix file!");
	}

	if (field_parser.search(1, "SysRHSVector")) {
		input_params.sys_rhs_vec_file = field_parser.next(
				input_params.sys_rhs_vec_file);
	} else {
		homemade_error_msg("Missing the system RHS vector file!");
	}

	if (field_parser.search(1, "OutputBase")) {
		input_params.output_base = field_parser.next(
				input_params.output_base);
	} else {
		homemade_error_msg("Missing the output filename base!");
	}

	if (field_parser.search(1, "SysEps")) {
		input_params.sys_eps = field_parser.next(
				input_params.sys_eps);
	} else {
		input_params.sys_eps = 1e-5;
	}

	if (field_parser.search(1, "SysIterDiv")) {
		input_params.sys_iter_div = field_parser.next(
				input_params.sys_iter_div);
	} else {
		input_params.sys_iter_div = 1000;
	}

	if (field_parser.search(1, "RBVectorBase")) {
		input_params.path_to_rb_vectors = field_parser.next(
				input_params.path_to_rb_vectors);
		input_params.bUseRBVectors = true;

		if (field_parser.search(1, "NbOfRBVectors")) {
			input_params.nb_of_rb_vectors = field_parser.next(
					input_params.nb_of_rb_vectors);
		} else {
			input_params.nb_of_rb_vectors = 6;
		}

	} else {
		input_params.bUseRBVectors = false;
	}
};
#endif /* LIBMESH_SOLVE_LINEAR_SYSTEM_INPUT_PARSER_H_ */
