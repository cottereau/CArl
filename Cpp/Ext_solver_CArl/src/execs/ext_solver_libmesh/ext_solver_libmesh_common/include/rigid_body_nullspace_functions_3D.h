/*
 * rigid_body_nullspace_functions_3D.h
 *
 *  Created on: Avr 19, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef RIGID_BODY_NULLSPACE_FUNCTIONS_3D_H_
#define RIGID_BODY_NULLSPACE_FUNCTIONS_3D_H_

#include "common_header_ext_solver_libmesh.h"
#include "PETSC_matrix_operations.h"

/// Build the rigid body modes associated to a given system
void build_rigid_body_vectors(libMesh::ImplicitSystem&  input_system, MatNullSpace& nullsp_sys);

/// Export the rigid body mode vectors to a folder
void write_rigid_body_vectors(MatNullSpace& nullsp_sys, const std::string output_base);

#endif /* RIGID_BODY_NULLSPACE_FUNCTIONS_3D_H_ */
