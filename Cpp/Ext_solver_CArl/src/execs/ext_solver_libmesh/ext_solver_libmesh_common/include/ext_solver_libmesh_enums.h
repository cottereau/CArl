/*
 * ext_solver_libmesh_enums.h
 *
 *  Created on: Apr 17, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef EXT_SOLVER_LIBMESH_ENUMS_H_
#define EXT_SOLVER_LIBMESH_ENUMS_H_

#include "common_header_ext_solver_libmesh.h"

/// Enumerate used to define which weight function must be used to assemble the system
enum   WeightFunctionSystemType {
	NO_WEIGHT = 0,
	MACRO = 1,
	MICRO = 2
};

/// Small enumerate defining a cube's boundary ID's
enum boundary_id_cube
{
		MIN_Z = 1,
		MIN_Y = 2,
		MAX_X = 3,
		MAX_Y = 4,
		MIN_X = 5,
		MAX_Z = 6
};

#endif /* EXT_SOLVER_LIBMESH_ENUMS_H_ */