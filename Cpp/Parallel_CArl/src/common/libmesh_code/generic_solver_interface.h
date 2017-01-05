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
	generic_solver_interface();

public:

	// Constructor
	generic_solver_interface(	const libMesh::Parallel::Communicator& comm) :
		m_comm { &comm }
	{

	};

	virtual KSPConvergedReason get_converged_reason() = 0;

	virtual void set_solver(libMesh::PetscMatrix<libMesh::Number>& matrix, const std::string name = "") = 0;

	virtual void set_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix) = 0;

	virtual void set_coupling_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix) = 0;

	virtual void set_rhs(libMesh::PetscVector<libMesh::Number>& vector) = 0;

	virtual libMesh::PetscMatrix<libMesh::Number>& get_matrix() = 0;

	virtual libMesh::PetscMatrix<libMesh::Number>& get_coupling_matrix() = 0;

	virtual libMesh::PetscVector<libMesh::Number>& get_rhs() = 0;

	// Get the system's dimensions
	virtual void get_system_dimensions(unsigned int& M_out, unsigned int& M_local_out) = 0;

	virtual void get_coupling_dimensions(unsigned int& M_out, unsigned int& N_out, unsigned int& M_local_out, unsigned int& N_local_out) = 0;

	// Calculate v_out = [ C * ( M^-1 ) *C^t ] * v_in
	virtual void apply_ZMiZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out) = 0;

	// Calculate v_out = [ ( M^-1 ) *C^t ] * v_in
	virtual void apply_MiZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out) = 0;

	// Calculate v_out = [ C * ( M^-1 ) ] * v_in
	void apply_ZMi(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// Calculate v_out = ( M^-1 ) * v_in
	virtual void solve(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out) = 0;

	// Calculate v_out = ( M^-1 ) * rhs
	virtual void solve(libMesh::PetscVector<libMesh::Number>& v_out) = 0;

	virtual void print_type() = 0;

	virtual void calculate_pseudo_inverse() = 0;
};

}



#endif /* COMMON_LIBMESH_CODE_GENERIC_SOLVER_INTERFACE_H_ */
