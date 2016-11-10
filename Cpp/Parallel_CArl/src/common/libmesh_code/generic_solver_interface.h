/*
 * generic_solver_interface.h
 *
 *  Created on: Nov 10, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_GENERIC_SOLVER_INTERFACE_H_
#define COMMON_LIBMESH_CODE_GENERIC_SOLVER_INTERFACE_H_

#include "carl_headers.h"

#include "PETSC_matrix_operations.h"

namespace carl
{

class generic_solver_interface
{

protected:

	// Communicator
	const libMesh::Parallel::Communicator * m_comm;

public:

	// Constructor
	generic_solver_interface(	const libMesh::Parallel::Communicator& comm) :
		m_comm { &comm }
	{

	};

	virtual void solve(libMesh::PetscVector<libMesh::Number>& sol, libMesh::PetscVector<libMesh::Number>& rhs);

	virtual void print_type();
};

}



#endif /* COMMON_LIBMESH_CODE_GENERIC_SOLVER_INTERFACE_H_ */
