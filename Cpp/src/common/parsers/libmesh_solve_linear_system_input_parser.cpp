/*
 * libmesh_solve_linear_system_input_parser.h
 *
 *  Created on: Apr 16, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#include "libmesh_solve_linear_system_input_parser.h"

namespace carl
{
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

    if (field_parser.search(1,"Newmark")) {
		input_params.alpha1 = field_parser.next(
				input_params.alpha1);
		input_params.beta1 = field_parser.next(
				input_params.beta1);
		input_params.gamma1 = field_parser.next(
				input_params.gamma1);
        //
		input_params.alpha2 = field_parser.next(
				input_params.alpha2);
		input_params.beta2 = field_parser.next(
				input_params.beta2);
		input_params.gamma2 = field_parser.next(
				input_params.gamma2);

    } else {
        input_params.alpha1 = 1.00;  
        input_params.alpha2 = 1.00;  
        input_params.beta1  = 0.25;  
        input_params.beta2  = 0.25;
        input_params.gamma1 = 0.50; 
        input_params.gamma2 = 0.50;  
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

void print_input_params(const std::string& output_filename,
		libmesh_solve_linear_system_input_params& input_params) {

	std::ofstream output_file(output_filename);

	output_file << "SysMatrix " << input_params.sys_matrix_file << std::endl;
	output_file << "SysRHSVector " << input_params.sys_rhs_vec_file << std::endl;
	output_file << "OutputBase " << input_params.output_base << std::endl;

	output_file << "SysEps " << input_params.sys_eps << std::endl;
	output_file << "SysIterDiv " << input_params.sys_iter_div << std::endl;

	if(input_params.bUseRBVectors)
	{
		output_file << "RBVectorBase " << input_params.path_to_rb_vectors << std::endl;
		output_file << "NbOfRBVectors " << input_params.nb_of_rb_vectors << std::endl;
	}

	output_file.close();
};

};
