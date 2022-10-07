/*
 * \file CArl_simple_1D_test_matrix_input_parser.h
 *
 *  Created on: Aug. 13, 2022
 *      Author: Chensheng Luo
 *
 * \brief  **TEST** This is the parser used for simple 1D code performance test case
 * 
 */

#ifndef CARL_TEST_SIMPLE_ONED_MATRIX_INPUT_PARSER_H_
#define CARL_TEST_SIMPLE_ONED_MATRIX_INPUT_PARSER_H_

#include "carl_headers.h"
#include "ext_solver_libmesh_enums.h"
#include "newmark_param_parser.h"

namespace carl
{

    struct simple_1D_test_matrix_input_params{
        int N;
        double e;
        double L;

        double rho;
        double E;
        double kappa;
        double es;

        double CM;
        double CK;

        carl::NewmarkParams newmark_A;
        carl::NewmarkParams newmark_B;

        std::string output_baseA;
        std::string output_baseB;

        std::string coupling_base;
    };

    void get_input_params(GetPot& field_parser,
        simple_1D_test_matrix_input_params& input_params);
}

#endif /* CARL_TEST_SIMPLE_ONED_MATRIX_INPUT_PARSER_H_ */