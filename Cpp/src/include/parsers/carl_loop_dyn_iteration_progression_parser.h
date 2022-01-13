/*
 * carl_loop_dyn_iteration_progression_parser.h 
 *
 *  Created on: Dec 9, 2021
 *      Author: Chensheng Luo
 */

#ifndef CARL_LOOP_DYN_PROGRE_PARSER_H_
#define CARL_LOOP_DYN_PROGRE_PARSER_H_

#include "carl_headers.h"

namespace carl
{
/// Structure containing the parameters for begining the dynamic loop
struct feti_loop_dyn_iteration_progression_params {
    // --- Parameters used directly for knowing the progression
    int inner_loop_progression;
    int outer_loop_progression;
};

void get_input_params(GetPot& field_parser,
        feti_loop_dyn_iteration_progression_params& input_params);


void print_input_params(const std::string& output_filename,
    feti_loop_dyn_iteration_progression_params& input_params);

}; 
#endif /* CARL_LOOP_DYN_PROGRE_PARSER_H_ */