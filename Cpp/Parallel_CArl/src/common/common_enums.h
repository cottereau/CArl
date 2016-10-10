/*
 * common_enums.h
 *
 *  Created on: Oct 6, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_ENUMS_H_
#define COMMON_ENUMS_H_

namespace carl
{
// Coupled system solver type
enum CoupledSolverType {
	LATIN_MODIFIED_STIFFNESS = 0, // Use LATIN's modified stiffness
	LATIN_ORIGINAL_STIFFNESS = 1 // Keep the stiffness unchanged
};
}





#endif /* COMMON_ENUMS_H_ */
