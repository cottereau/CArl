/*
 * carl_loop_dyn_force_parser.h
 *
 *  Created on: June 24, 2022
 *      Author: Chensheng Luo
 */

#ifndef CARL_LOOP_DYN_FORCE_PARSER_H_
#define CARL_LOOP_DYN_FORCE_PARSER_H_

#include "carl_headers.h"

namespace carl
{

    struct dyn_force_params{

        // Common parameters
        std::string modal_A;
        std::string modal_B;

        // Optional parameters
        double amplitude_A;
        double amplitude_B;
        double frequency_A;
        double frequency_B;
        double initialPhase_A;
        double initialPhase_B;
        double slope_A;
        double slope_B;

        double offset_A;
        double offset_B;
        double saturation_A;
        double saturation_B;

        std::string time_series_A;
        std::string time_series_B;
        int interpolation_method_A;
        int interpolation_method_B;
    };

/** \brief **DYN** Parser function for force preparation input.
 *  
Required parameters:

 + **Modal Force**: 
     - `ModalVectorA` : Modal Force of object A
     - `ModalVectorB` : Modal Force of object B

Following paramters are based on cases, according to `ForcePrepareMethod` in dynmaic solver input.   
If no object is specified, it will be applied to all parameters!

+ `ForcePrepareMethod` = `ModalSinus`
    - `Amplitude` OR (`AmplitudeA` and `AmplitudeB`), *Default*: 1
    - `Frequency` OR (`FrequencyA` and `FrequencyB`), *Default*: 1
    - `IniitalPhase` OR (`IniitalPhaseA` and `IniitalPhaseB`), *Default*: 0

    It will create: \f$ \text{modal} \times \text{Amplitude}\sin(2\pi\times\text{Frequency}\times t+\text{IniitalPhase}) \f$

+ `ForcePrepareMethod` = `ModalConstant`
    - `Amplitude` OR (`AmplitudeA` and `AmplitudeB`), *Default*: 1

    It will create: \f$ \text{modal} \times \text{Amplitude} \f$

+ `ForcePrepareMethod` = `ModalLinear`
    - `Slope` OR (`SlopeA` and `SlopeB`), *Default*: 1
    - `Saturation` OR (`SaturationA` and `SaturationB`), *Default*: INFINITY
    - `Offset` OR (`OffsetA` and `OffsetB`), *Default*: 0
    
    It will create: \f$ \text{modal} \times min(\text{Slope} \times t + Offset, Saturation) \f$

+ `ForcePrepareMethod` = `ModalProduct`
    - `TimeSeries` OR (`TimeSeriesA` and `TimeSeriesB`), **ATTENTION** In the time series file, time and acceleration but be given, and it should be in chronological order
    - `Interpolation` OR (`InterpolationA` and `InterpolationB`), *Default*: Off, *chosen from* `Linear`, `nearest`

*/
void get_input_params(GetPot& field_parser,int force_prepare_mode,dyn_force_params& input_params);
};
#endif /* CARL_LOOP_DYN_FORCE_PARSER_H_ */