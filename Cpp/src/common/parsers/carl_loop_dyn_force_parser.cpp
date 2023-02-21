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

  std::cout << " -- Reading force prapare file Begin ......" << std::endl;
  // Common files
  if (field_parser.search(1, "ModalVectorA")) {
    input_params.modal_A = field_parser.next(
        input_params.modal_A);
    std::cout << " -- Model Force A: " << input_params.modal_A << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing modal force A command line!");
  }

  if (field_parser.search(1, "ModalVectorB")) {
    input_params.modal_B = field_parser.next(
        input_params.modal_B);
    std::cout << " -- Model Force B: " << input_params.modal_B << std::endl;
  } else {
    homemade_error_msg("[CArl Parameters]ERROR! Missing modal force B command line!");
  }

  if(force_prepare_mode == ForcePrepareMethod::MODAL_SINUS){
          std::cout << " -- Method: Sinusoidal" << std::endl;
          std::cout << " #### F(t) = Amplitude * sin(2*pi*Frequency*t + InitialPhase) * Modal #####" << std::endl;
          if (field_parser.search(2, "AmplitudeA","Amplitude")) {
            input_params.amplitude_A = field_parser.next(
                input_params.amplitude_A);
            std::cout << " -- Amplitude A: " << input_params.amplitude_A << std::endl;
          }else{
            input_params.amplitude_A = 1;
            std::cout << " -- Amplitude A: " << input_params.amplitude_A << "[No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "AmplitudeB","Amplitude")) {
            input_params.amplitude_B = field_parser.next(
                input_params.amplitude_B);
            std::cout << " -- Amplitude B: "<< input_params.amplitude_B << std::endl;
          }else{
            input_params.amplitude_B = 1;
            std::cout << " -- Amplitude B: "<< input_params.amplitude_B << "[No entry, default value is taken]" << std::endl;
          }


          if (field_parser.search(2, "FrequencyA", "Frequency")) {
            input_params.frequency_A = field_parser.next(
                input_params.frequency_A);
            std::cout << " -- Frequency A: "<<input_params.frequency_A << std::endl;
          } else{
            input_params.frequency_A = 1;
            std::cout << " -- Frequency A: "<<input_params.frequency_A << "[No entry, default value is taken]" << std::endl;
          }

          if (field_parser.search(2, "FrequencyB", "Frequency")) {
            input_params.frequency_B = field_parser.next(
                input_params.frequency_B);
            std::cout << " -- Frequency B: " << input_params.frequency_B << std::endl;
          } else{
            input_params.frequency_B = 1;
            std::cout << " -- Frequency B: " << input_params.frequency_B << "[No entry, default value is taken]" << std::endl;
          }

          if (field_parser.search(2, "InitialPhaseA", "IniitalPhase")) {
            input_params.initialPhase_A = field_parser.next(
                input_params.initialPhase_A);
            std::cout << " -- InitialPhase A: " << input_params.initialPhase_A << std::endl;
          } else{
            input_params.initialPhase_A = 0;
            std::cout << " -- InitialPhase A: " << input_params.initialPhase_A << "[No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "IniitalPhaseB","IniitalPhase")) {
            input_params.initialPhase_B = field_parser.next(
                input_params.initialPhase_B);
            std::cout << " -- InitialPhase B: "<< input_params.initialPhase_B << std::endl;
          } else{
            input_params.initialPhase_B = 0;
            std::cout << " -- InitialPhase B: " << input_params.initialPhase_B  << "[No entry, default value is taken]"<< std::endl;
          }
  }else if (force_prepare_mode == ForcePrepareMethod::MODAL_CONSTANT){

          std::cout << " -- Method: Constant " << std::endl;
          std::cout << " #### F(t) = Constant * Modal #####" << std::endl;

          if (field_parser.search(2, "AmplitudeA","Amplitude")) {
            input_params.amplitude_A = field_parser.next(
                input_params.amplitude_A);
            std::cout << " -- Amplitude A: " << input_params.amplitude_A << std::endl;
          }else{
            input_params.amplitude_A = 1;
            std::cout << " -- Amplitude A: " << input_params.amplitude_A << "[No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "AmplitudeB","Amplitude")) {
            input_params.amplitude_B = field_parser.next(
                input_params.amplitude_B);
            std::cout << " -- Amplitude B: "<< input_params.amplitude_B << std::endl;
          }else{
            input_params.amplitude_B = 1;
            std::cout << " -- Amplitude B: "<< input_params.amplitude_B << "[No entry, default value is taken]" << std::endl;
          }

  }else if (force_prepare_mode == ForcePrepareMethod::MODAL_LINEAR){
          
          std::cout << " -- Method: Linear " << std::endl;
          std::cout << " #### F(t) = min(Slope * t + Offset, Saturation) * Modal #####" << std::endl;

          if (field_parser.search(2, "SlopeA", "Slope")) {
            input_params.slope_A = field_parser.next(
                input_params.slope_A);
            std::cout << " -- Slope A: "<< input_params.slope_A << std::endl;
          } else{
            input_params.slope_A = 1;
            std::cout << " -- Slope A: "<< input_params.slope_A << "[No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "SlopeB", "Slope")) {
            input_params.slope_B = field_parser.next(
                input_params.slope_B);
            std::cout << " -- Slope B: "<< input_params.slope_B << std::endl;
          } else{
            input_params.slope_B = 1;
            std::cout << " -- Slope B: "<< input_params.slope_B << "[No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "SaturationA", "Saturation")) {
            input_params.saturation_A = field_parser.next(
                input_params.saturation_A);
            std::cout << " -- Saturation A: "<< input_params.saturation_A << std::endl;
          } else{
            input_params.saturation_A = INFINITY;
            std::cout << " -- Saturation A: INFINITY [No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "SaturationB", "Saturation")) {
            input_params.saturation_B = field_parser.next(
                input_params.saturation_B);
            std::cout << " -- Saturation B: "<< input_params.saturation_B << std::endl;
          } else{
            input_params.saturation_B = INFINITY;
            std::cout << " -- Saturation B: INFINITY [No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "OffsetA", "Offset")) {
            input_params.offset_A = field_parser.next(
                input_params.offset_A);
            std::cout << " -- Offset A: "<< input_params.offset_A << std::endl;
          } else{
            input_params.offset_A = 0;
            std::cout << " -- Offset A: "<< input_params.offset_A << "[No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "OffsetB", "Offset")) {
            input_params.offset_B = field_parser.next(
                input_params.offset_B);
            std::cout << " -- Offset B: "<< input_params.offset_B << std::endl;
          } else{
            input_params.offset_B = 0;
            std::cout << " -- Offset B: "<< input_params.offset_B << "[No entry, default value is taken]"<< std::endl;
          }


    }else if (force_prepare_mode == ForcePrepareMethod::MODAL_PRODUCT){
          std::cout << " -- Method: Product " << std::endl;
          std::cout << " #### F(t) = TimeSeries(t) * Modal #####" << std::endl;

          if (field_parser.search(2, "TimeSeriesA", "TimeSeries")) {
            input_params.time_series_A = field_parser.next(
                input_params.time_series_A);
            std::cout << " -- TimeSeries A: "<< input_params.time_series_A << std::endl;
          } else{
           homemade_error_msg("[CArl Parameters]ERROR! Missing time series file for A (needed for the PRODUCT force preparation method)!");
          }

          if (field_parser.search(2, "TimeSeriesB", "TimeSeries")) {
            input_params.time_series_B = field_parser.next(
                input_params.time_series_B);
            std::cout << " -- TimeSeries B: "<< input_params.time_series_B << std::endl;
          } else{
           homemade_error_msg("[CArl Parameters]ERROR! Missing time series file for B (needed for the PRODUCT force preparation method)!");
          }

          if (field_parser.search(2, "InterpolationA", "Interpolation")) {
            std::string interpolation_method;
            interpolation_method = field_parser.next(interpolation_method);
            std::transform(interpolation_method.begin(),interpolation_method.end(),interpolation_method.begin(),[](unsigned char c){ return std::tolower(c); });
            if(interpolation_method == "linear"){
              input_params.interpolation_method_A = carl::ForceInterpolationMethod::LINEAR;
              std::cout << "- InterpolationMethodA: Linear"<< std::endl;
            }else if(interpolation_method == "nearest"){
              input_params.interpolation_method_A = carl::ForceInterpolationMethod::NEAREST;
              std::cout << "- InterpolationMethodA: Nearest"<< std::endl;
            }else{
             homemade_error_msg("[CArl Parameters]ERROR! Invalid interpolation method for A!");
            }
          } else{
            input_params.interpolation_method_A = carl::ForceInterpolationMethod::OFF;
            std::cout << " -- InterpolationMethodA: Off[No entry, default value is taken]"<< std::endl;
          }

          if (field_parser.search(2, "InterpolationB", "Interpolation")) {
            std::string interpolation_method;
            interpolation_method = field_parser.next(interpolation_method);
            std::transform(interpolation_method.begin(),interpolation_method.end(),interpolation_method.begin(),[](unsigned char c){ return std::tolower(c); });
            if(interpolation_method == "linear"){
              input_params.interpolation_method_B = carl::ForceInterpolationMethod::LINEAR;
              std::cout << "- InterpolationMethodA: Linear"<< std::endl;
            }else if(interpolation_method == "nearest"){
              input_params.interpolation_method_B = carl::ForceInterpolationMethod::NEAREST;
              std::cout << "- InterpolationMethodA: Nearest"<< std::endl;
            }else{
              homemade_error_msg("[CArl Parameters]ERROR! Invalid interpolation method for A!");
            }
          } else{
            input_params.interpolation_method_B = carl::ForceInterpolationMethod::OFF;
            std::cout << " -- InterpolationMethodB: Off[No entry, default value is taken]"<< std::endl;
          }
    }

    std::cout << " ...... Reading force files: FINISH!  " << std::endl;

};
}