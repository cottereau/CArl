/*
 * PETSC_matrix_operations.h
 *
 *  Created on: Jan 29, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_PETSC_MATRIX_OPERATIONS_H_
#define COMMON_LIBMESH_CODE_PETSC_MATRIX_OPERATIONS_H_

#include "common_header.h"
#include "common_header_libmesh.h"
#include "common_functions.h"

namespace carl
{

void lump_matrix(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
						libMesh::PetscMatrix<libMesh::Number>& matrixOutput);

void lump_matrix_and_invert(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
									libMesh::PetscMatrix<libMesh::Number>& matrixOutput);

void lump_matrix_and_invert(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
									libMesh::PetscVector<libMesh::Number>& vecOutput);

void print_matrix(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix);

void print_matrix_dim(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix);

void print_matrix_col_line_sum(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix, const std::string name_base);

void print_matrix_matlab(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix, const std::string name_base);

void solve_linear_PETSC(	libMesh::PetscMatrix<libMesh::Number>& A,
							libMesh::PetscVector<libMesh::Number>& b,
							libMesh::PetscVector<libMesh::Number>& x,
							KSP& ksp, PC& pc);

}

#endif /* COMMON_LIBMESH_CODE_PETSC_MATRIX_OPERATIONS_H_ */
