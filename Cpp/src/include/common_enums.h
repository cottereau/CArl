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

enum ClusterSchedulerType {
	LOCAL = 0, // No scheduler present, will use carl::exec_command
	PBS = 1, // PBS / Torque
	SLURM = 2 // SLURM, not implemented right now
};


// External solver types
enum ExtSolverType {
	LIBMESH_LINEAR = 0,
	DUMMY = 1,
	LIBMESH_DYNAMIC = 2,
	LIBMESH_NONLINGEOM = 3
};

// Mediator space type
enum MediatorType {
	USE_MACRO = 0,
	USE_MICRO = 1,
	USE_EXTERNAL = 2
};

// CG preconditioner type
enum BaseCGPrecondType {
	NO_PRECONDITIONER = 0, // Identity matrix
	COUPLING_OPERATOR = 1, // C_RR
	COUPLING_JACOBI = 2 // diagonal(C_RR)
};

enum IterationStatus {
	ITERATING = 0,
	DIVERGED = 1,
	CONVERGED = 2
};

enum IntersectionMeshingMethod {
	LIBMESH_TETGEN = 0, // libMesh Tetgen algorithm, problematic with
						// Intel compilers
	CGAL = 1 			// Intersection meshing algorithm using
						// CGAL's Triangulation_3
};

enum SearchMethod
{
	BRUTE = 0,
	FRONT = 1,
	BOTH = 2
};

enum RBModesSystem {
	MACRO = 0,
	MICRO = 1
};
}





#endif /* COMMON_ENUMS_H_ */
