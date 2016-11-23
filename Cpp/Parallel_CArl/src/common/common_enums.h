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
	LATIN_MODIFIED_STIFFNESS = 0, // Original LATIN solver
	LATIN_ORIGINAL_STIFFNESS = 1, // Original LATIN solver, with unmodified K's
	CG    = 2  // Conjugate gradient solver
};

enum IntersectionMeshingMethod {
	LIBMESH_TETGEN = 0, // libMesh Tetgen algorithm, problematic with
						// Intel compilers
	CGAL = 1 			// Intersection meshing algorithm using
						// CGAL's Triangulation_3
};
}





#endif /* COMMON_ENUMS_H_ */
