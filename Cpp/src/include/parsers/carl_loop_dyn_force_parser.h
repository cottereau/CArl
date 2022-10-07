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


};

/** \brief **DYN-DI/DYN-CG** Parser function for force preparation input.
 *  
Required parameters:

 + **Modal Force**: 
     - `ModalVectorA` : Modal Force of object A
     - `ModalVectorB` : Modal Force of object B

Following paramters are based on cases, according to `ForcePrepareMethod` in dynmaic solver input.   
If no object is specified, it will be applied to all parameters!

+ `ForcePrepareMethod` = `ModalSinus`
    - `Amplitude` OR (`AmplitudeA` and `AmplitudeB`)
    - `Frequency` OR (`FrequencyA` and `FrequencyB`)
    - `IniitalPhase` OR (`IniitalPhaseA` and `IniitalPhaseB`)  

    It will create: \f$ \text{modal} \times \text{Amplitude}\sin(2\pi\times\text{Frequency}\times t+\text{IniitalPhase}) \f$

+ `ForcePrepareMethod` = `ModalConstant`
    - `Amplitude` OR (`AmplitudeA` and `AmplitudeB`)  

    It will create: \f$ \text{modal} \times \text{Amplitude} \f$

+ `ForcePrepareMethod` = `ModalLinear`
    - `Slope` OR (`SlopeA` and `SlopeB`)  
    
    It will create: \f$ \text{modal} \times \text{Slope} \times t \f$

+ `ForcePrepareMethod` = `ModalProduct`
    - NOT IMPLEMENTED

*/
void get_input_params(GetPot& field_parser,int force_prepare_mode,dyn_force_params& input_params);
};
#endif /* CARL_LOOP_DYN_FORCE_PARSER_H_ */