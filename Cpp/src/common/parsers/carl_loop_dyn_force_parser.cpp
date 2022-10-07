/*
 * carl_loop_dyn_force_parser.cpp
 *
 *  Created on: June 24, 2022
 *      Author: Chensheng Luo
 */

#include "carl_loop_dyn_force_parser.h"

namespace carl
{

void get_input_params(GetPot& field_parser,
    int force_prepare_mode,
    dyn_force_params& input_params) {

  // Common files
  if (field_parser.search(1, "ModalVectorA")) {
    input_params.modal_A = field_parser.next(
        input_params.modal_A);
    std::cout << input_params.modal_A << std::endl;
  } else {
    homemade_error_msg("Missing modal force A command line!");
  }

  if (field_parser.search(1, "ModalVectorB")) {
    input_params.modal_B = field_parser.next(
        input_params.modal_B);
    std::cout << input_params.modal_B << std::endl;
  } else {
    homemade_error_msg("Missing modal force B command line!");
  }

  if(force_prepare_mode == ForcePrepareMethod::MODAL_SINUS){
          if (field_parser.search(1, "AmplitudeA")) {
            input_params.amplitude_A = field_parser.next(
                input_params.amplitude_A);
            std::cout << input_params.amplitude_A << std::endl;
          } else if(field_parser.search(1, "Amplitude")) {
            input_params.amplitude_A = field_parser.next(
                input_params.amplitude_A);
            std::cout << input_params.amplitude_A << std::endl;
          } else{
            input_params.amplitude_A = 1;
          }

          if (field_parser.search(1, "AmplitudeB")) {
            input_params.amplitude_B = field_parser.next(
                input_params.amplitude_B);
            std::cout << input_params.amplitude_B << std::endl;
          } else if(field_parser.search(1, "Amplitude")) {
            input_params.amplitude_B = field_parser.next(
                input_params.amplitude_B);
            std::cout << input_params.amplitude_B << std::endl;
          } else{
            input_params.amplitude_B = 1;
          }


          if (field_parser.search(1, "FrequencyA")) {
            input_params.frequency_A = field_parser.next(
                input_params.frequency_A);
            std::cout << input_params.frequency_A << std::endl;
          } else if(field_parser.search(1, "Frequency")) {
            input_params.frequency_A = field_parser.next(
                input_params.frequency_A);
            std::cout << input_params.amplitude_A << std::endl;
          }else{
            input_params.frequency_A = 1;
          }

          if (field_parser.search(1, "FrequencyB")) {
            input_params.frequency_B = field_parser.next(
                input_params.frequency_B);
            std::cout << input_params.frequency_B << std::endl;
          } else if(field_parser.search(1, "Frequency")) {
            input_params.frequency_B = field_parser.next(
                input_params.frequency_B);
            std::cout << input_params.frequency_B << std::endl;
          } else{
            input_params.frequency_B = 1;
          }

          if (field_parser.search(1, "InitialPhaseA")) {
            input_params.initialPhase_A = field_parser.next(
                input_params.initialPhase_A);
            std::cout << input_params.initialPhase_A << std::endl;
          } else if(field_parser.search(1, "IniitalPhase")) {
            input_params.initialPhase_A = field_parser.next(
                input_params.initialPhase_A);
            std::cout << input_params.amplitude_A << std::endl;
          }else{
            input_params.initialPhase_A = 0;
          }

          if (field_parser.search(1, "IniitalPhaseB")) {
            input_params.initialPhase_B = field_parser.next(
                input_params.initialPhase_B);
            std::cout << input_params.initialPhase_B << std::endl;
          } else if(field_parser.search(1, "IniitalPhase")) {
            input_params.initialPhase_B = field_parser.next(
                input_params.initialPhase_B);
            std::cout << input_params.initialPhase_B << std::endl;
          } else{
            input_params.initialPhase_B = 0;
          }
  }else if (force_prepare_mode == ForcePrepareMethod::MODAL_CONSTANT){

          if (field_parser.search(1, "AmplitudeA")) {
            input_params.amplitude_A = field_parser.next(
                input_params.amplitude_A);
            std::cout << input_params.amplitude_A << std::endl;
          } else if(field_parser.search(1, "Amplitude")) {
            input_params.amplitude_A = field_parser.next(
                input_params.amplitude_A);
            std::cout << input_params.amplitude_A << std::endl;
          } else{
            input_params.amplitude_A = 1;
          }

          if (field_parser.search(1, "AmplitudeB")) {
            input_params.amplitude_B = field_parser.next(
                input_params.amplitude_B);
            std::cout << input_params.amplitude_B << std::endl;
          } else if(field_parser.search(1, "Amplitude")) {
            input_params.amplitude_B = field_parser.next(
                input_params.amplitude_B);
            std::cout << input_params.amplitude_B << std::endl;
          } else{
            input_params.amplitude_B = 1;
          }

  }else if (force_prepare_mode == ForcePrepareMethod::MODAL_LINEAR){
          if (field_parser.search(1, "SlopeA")) {
            input_params.slope_A = field_parser.next(
                input_params.slope_A);
            std::cout << input_params.slope_A << std::endl;
          } else if(field_parser.search(1, "Slope")) {
            input_params.slope_A = field_parser.next(
                input_params.slope_A);
            std::cout << input_params.slope_A << std::endl;
          } else{
            input_params.slope_A = 1;
          }

          if (field_parser.search(1, "SlopeB")) {
            input_params.slope_B = field_parser.next(
                input_params.slope_B);
            std::cout << input_params.slope_B << std::endl;
          } else if(field_parser.search(1, "Slope")) {
            input_params.slope_B = field_parser.next(
                input_params.slope_B);
            std::cout << input_params.slope_B << std::endl;
          } else{
            input_params.slope_B = 1;
          }
    }else if (force_prepare_mode == ForcePrepareMethod::MODAL_PRODUCT){
        homemade_error_msg("Not Implemented");
    }

    std::cout << "Force parser finish" << std::endl;

};
}