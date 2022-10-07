/*
 * newmark_param_parser.h
 *
 *  Created on: August 3,2022
 *      Author: Chensheng Luo
 */

#ifndef NEWMARK_PARAM_PARSER_H_
#define NEWMARK_PARAM_PARSER_H_

#include "carl_headers.h"

namespace carl
{

struct NewmarkParams {
    double alpha;
    double beta;
    double gamma;
    double deltat;

};

/** \brief **DYN-DI/DYN-CG** Parser function for newmark parameter.
 *  
Required parameters:
	- `deltat`, `Deltat` or `DELTAT`, default to be 0.001
	- `alpha`, `Alpha` or `ALPHA`, default to be 0
	- `beta`, `Beta` or `BETA`, default to be 0.25
	- `gamma`, `Gamma` or `GAMMA`, default to be 0.5

*/
void get_newmark_params(GetPot& field_parser,NewmarkParams& newmark);
};
#endif /* NEWMARK_PARAM_PARSER_H_ */