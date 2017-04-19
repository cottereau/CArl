/*
 * PETSC_matrix_operations.h
 *
 *  Created on: Jan 29, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_PETSC_MATRIX_OPERATIONS_H_
#define COMMON_LIBMESH_CODE_PETSC_MATRIX_OPERATIONS_H_

#include "carl_headers.h"

namespace carl
{

void lump_matrix(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
						libMesh::PetscMatrix<libMesh::Number>& matrixOutput);

void lump_matrix_and_invert(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
									libMesh::PetscMatrix<libMesh::Number>& matrixOutput);

void lump_matrix_and_invert(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
									libMesh::PetscVector<libMesh::Number>& vecOutput);

void print_matrix(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix);

void print_matrix_dim(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix, bool bDetailed = false);

void print_matrix_info(libMesh::PetscMatrix<libMesh::Number>& InputMatrix, std::ostream & os=libMesh::out);

void print_matrix_col_line_sum(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix, const std::string name_base);

void print_matrix_matlab(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix, const std::string name_base);

void solve_linear_PETSC(	libMesh::PetscMatrix<libMesh::Number>& A,
							libMesh::PetscVector<libMesh::Number>& b,
							libMesh::PetscVector<libMesh::Number>& x,
							KSP& ksp, PC& pc);

void check_coupling_matrix( libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix ,
		libMesh::Mesh& IntersectionMesh,
		libMesh::Real CouplingScale,
		const std::string matrixtype,
		int n_var = 3);

void write_PETSC_vector(	libMesh::PetscVector<libMesh::Number>& input_vec,
										const std::string& filename);

void read_PETSC_vector(	libMesh::PetscVector<libMesh::Number>& input_vec,
										const std::string& filename);

void print_PETSC_vector(	libMesh::PetscVector<libMesh::Number>& input_vec,
		const std::string& filename);

void write_PETSC_vector( Vec input_vec, const std::string& filename, MPI_Comm comm = PETSC_COMM_WORLD);

void read_PETSC_vector(	Vec input_vec, const std::string& filename, MPI_Comm comm = PETSC_COMM_WORLD);

void write_PETSC_matrix( Mat input_mat, const std::string& filename, MPI_Comm comm = PETSC_COMM_WORLD);

void read_PETSC_matrix(	Mat input_mat, const std::string& filename, MPI_Comm comm = PETSC_COMM_WORLD);

void create_PETSC_dense_matrix_from_vectors(const Vec *vecs_in, int nb_vecs, Mat& matrix_out);

void PETSC_invert_dense_matrix(Mat& matrix_in, Mat& matrix_out);
}



#endif /* COMMON_LIBMESH_CODE_PETSC_MATRIX_OPERATIONS_H_ */
