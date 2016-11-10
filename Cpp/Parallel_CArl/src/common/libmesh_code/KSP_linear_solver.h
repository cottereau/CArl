/*
 * KSP_linear_solver.h
 *
 *  Created on: Nov 10, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_KSP_LINEAR_SOLVER_H_
#define COMMON_LIBMESH_CODE_KSP_LINEAR_SOLVER_H_

#include "carl_headers.h"

#include "PETSC_matrix_operations.h"

#include "generic_solver_interface.h"

namespace carl
{

class KSP_linear_solver : public generic_solver_interface
{

protected:
	std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> > m_KSP_solver;

	double 			m_KSP_eps;
	unsigned int 	m_KSP_iter_max;

	libMesh::PetscMatrix<libMesh::Number> * m_Matrix;


public:

	// Constructor
	KSP_linear_solver(	const libMesh::Parallel::Communicator& comm) :
		m_comm { &comm },

		m_KSP_eps { 1E-8 },
		m_KSP_iter_max { 10000 },

		m_Matrix { NULL }
	{
		m_KSP_solver = std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> >
			(new libMesh::PetscLinearSolver<libMesh::Number>(*m_comm));
	};

	// Specific methods
	void set_solver(libMesh::PetscMatrix<libMesh::Number>& matrix, const std::string name = "");

	void set_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix);

	// Implementations of virtual methods
	void solve(libMesh::PetscVector<libMesh::Number>& sol, libMesh::PetscVector<libMesh::Number>& rhs);

	void print_type();
};

}




#endif /* COMMON_LIBMESH_CODE_KSP_LINEAR_SOLVER_H_ */
