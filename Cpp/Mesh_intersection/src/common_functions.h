/*
 * common_functions.h
 *
 *  Created on: Jan 26, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_COMMON_FUNCTIONS_H_
#define COMMON_COMMON_FUNCTIONS_H_

#include "common_header_libmesh.h"

libMesh::Real kronecker_delta(unsigned int i,
				   unsigned int j);

void clear_line();
#endif /* COMMON_COMMON_FUNCTIONS_H_ */
