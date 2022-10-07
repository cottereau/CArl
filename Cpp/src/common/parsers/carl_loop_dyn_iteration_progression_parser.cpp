/*
 * carl_loop_dyn_iteration_progression_parser.cpp
 *
 *  Created on: Dec 9, 2021
 *      Author: Chensheng Luo
 */

#include "carl_loop_dyn_iteration_progression_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
    feti_loop_dyn_iteration_progression_params& input_params) {

  if (field_parser.search(1, "InnerProgression")) {
    input_params.inner_loop_progression = field_parser.next(input_params.inner_loop_progression);
  } else {
    homemade_error_msg("[ITERATION COUNTER]Unable to know inner loop progression!");
  }

  if (field_parser.search(1, "OuterProgression")) {
    input_params.outer_loop_progression = field_parser.next(input_params.outer_loop_progression);
  } else {
    homemade_error_msg("[ITERATION COUNTER]Unable to know outer loop progression!");
  }
}

void print_input_params(const std::string& output_filename,
    feti_loop_dyn_iteration_progression_params& input_params) {

  std::ofstream output_file(output_filename);

  output_file << "InnerProgression " << input_params.inner_loop_progression << std::endl;
  output_file << "OuterProgression " << input_params.outer_loop_progression << std::endl;

  output_file.close();
}
}