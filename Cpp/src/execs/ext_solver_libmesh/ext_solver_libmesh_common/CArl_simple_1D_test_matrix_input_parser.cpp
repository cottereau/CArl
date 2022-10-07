/*
 * \file CArl_simple_1D_test_matrix_input_parser.cpp
 *
 *  Created on: Aug. 13, 2022
 *      Author: Chensheng Luo
 *
 * \brief  **TEST** This is the parser used for simple 1D code performance test case
 * 
 */

#include "CArl_simple_1D_test_matrix_input_parser.h"

namespace carl
{
void get_input_params(GetPot& field_parser,
		simple_1D_test_matrix_input_params& input_params) {

    if (field_parser.search(1, "N")) {
    	input_params.N = field_parser.next(input_params.N);
  	} else {
    	input_params.N = 70;
  	}
  	std::cout << input_params.N << std::endl;

  	if (field_parser.search(1, "e")) {
    	input_params.e = field_parser.next(input_params.e);
  	} else {
    	input_params.e = 0.5;
  	}
  	std::cout << input_params.e << std::endl;


  	if (field_parser.search(1, "L")) {
    	input_params.L = field_parser.next(input_params.L);
  	} else {
    	input_params.L = 3.5;
  	}
  	std::cout << input_params.L << std::endl;

  	if (field_parser.search(1, "rho")) {
    	input_params.rho = field_parser.next(input_params.rho);
  	} else {
    	input_params.rho = 1000;
  	}
  	std::cout << input_params.rho << std::endl;

  	if (field_parser.search(1, "Young")) {
    	input_params.E = field_parser.next(input_params.E);
  	} else {
    	input_params.E = 200000;
  	}
  	std::cout << input_params.E << std::endl;

  	if (field_parser.search(1, "kappa")) {
    	input_params.kappa = field_parser.next(input_params.kappa);
  	} else {
    	input_params.kappa = input_params.E;
  	}
  	std::cout << input_params.kappa << std::endl;

  	if (field_parser.search(1, "es")) {
    	input_params.es = field_parser.next(input_params.es);
  	} else {
    	input_params.es = input_params.e;
  	}
  	std::cout << input_params.es << std::endl;

  	if (field_parser.search(1, "CM")) {
    	input_params.CM = field_parser.next(input_params.CM);
  	} else {
    	input_params.CM = 0;
  	}
  	std::cout << input_params.CM << std::endl;

  	if (field_parser.search(1, "CK")) {
    	input_params.CK = field_parser.next(input_params.CK);
  	} else {
    	input_params.CK = 0;
  	}
  	std::cout << input_params.CK << std::endl;


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
    std::cout << "Newmark Parameter (A): "<< input_params.newmark_A.alpha << "," << input_params.newmark_A.beta << "," << input_params.newmark_A.gamma << "," << input_params.newmark_A.deltat << std::endl;
    std::cout << "Newmark Parameter (B): "<< input_params.newmark_B.alpha << "," << input_params.newmark_B.beta << "," << input_params.newmark_B.gamma << "," << input_params.newmark_B.deltat << std::endl;

    if (field_parser.search(1, "SystemBaseA")) {
    	input_params.output_baseA = field_parser.next(input_params.output_baseA);
  	} else {
    	homemade_error_msg("[CArl Parameters]Missing A system matrix output base!");
  	}
  	std::cout << input_params.output_baseA << std::endl;

  	if (field_parser.search(1, "SystemBaseB")) {
    	input_params.output_baseB = field_parser.next(input_params.output_baseB);
  	} else {
    	homemade_error_msg("[CArl Parameters]Missing B system matrix output base!");
  	}
  	std::cout << input_params.output_baseB << std::endl;


  	if (field_parser.search(1, "CouplingBase")) {
    	input_params.coupling_base = field_parser.next(input_params.coupling_base);
  	} else {
    	homemade_error_msg("[CArl Parameters]Missing coupling matrix output base!");
  	}
  	std::cout << input_params.coupling_base << std::endl;


};

};
