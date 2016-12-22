/*
 * KSP_linear_solver.h
 *
 *  Created on: Nov 10, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_KSP_LINEAR_SOLVER_H_
#define COMMON_LIBMESH_CODE_KSP_LINEAR_SOLVER_H_

#include "carl_headers.h"
#include "generic_solver_interface.h"

#include "PETSC_matrix_operations.h"

namespace carl
{

class KSP_linear_solver : public generic_solver_interface
{

protected:
	std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> > m_KSP_solver;

	double 			m_KSP_eps;
	unsigned int 	m_KSP_iter_max;

	libMesh::PetscMatrix<libMesh::Number> * m_Matrix;
	libMesh::PetscMatrix<libMesh::Number> * m_Coupling;
	libMesh::PetscVector<libMesh::Number> * m_rhs;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_aux_vec_1;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_aux_vec_2;

	bool m_bRhsSet;
	bool m_bMatrixSet;
	bool m_bCouplingSet;

	KSPConvergedReason m_conv_reason;

public:

	// Constructor
	KSP_linear_solver(	const libMesh::Parallel::Communicator& comm) :
		generic_solver_interface(comm),
		m_KSP_eps { 1E-8 },
		m_KSP_iter_max { 10000 },

		m_Matrix { NULL },
		m_Coupling { NULL },
		m_rhs { NULL },
		m_bRhsSet { false },
		m_bMatrixSet { false },
		m_bCouplingSet { false }
	{
		m_KSP_solver = std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> >
			(new libMesh::PetscLinearSolver<libMesh::Number>(*m_comm));
		m_aux_vec_1 = std::unique_ptr<libMesh::PetscVector<libMesh::Number> >
			(new libMesh::PetscVector<libMesh::Number>(*m_comm));
		m_aux_vec_2 = std::unique_ptr<libMesh::PetscVector<libMesh::Number> >
			(new libMesh::PetscVector<libMesh::Number>(*m_comm));
	};

	// Pointers
	libMesh::PetscLinearSolver<libMesh::Number> * solver_ptr();

	libMesh::PetscMatrix<libMesh::Number> * matrix_ptr();

	// Specific methods
	void set_solver(libMesh::PetscMatrix<libMesh::Number>& matrix, const std::string name = "");

	void set_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix);

	void set_coupling_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix);

	void set_rhs(libMesh::PetscVector<libMesh::Number>& vector);

	KSPConvergedReason get_converged_reason();

	// Implementations of virtual methods

	// Get the system's dimensions
	void get_system_dimensions(unsigned int& M_out, unsigned int& M_local_out);

	void get_coupling_dimensions(unsigned int& M_out, unsigned int& N_out, unsigned int& M_local_out, unsigned int& N_local_out);

	// Calculate v_out = [ C * ( M^-1 ) *C^t ] * v_in
	void apply_ZMiZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// Calculate v_out = [ ( M^-1 ) *C^t ] * v_in
	void apply_MiZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// Calculate v_out = ( M^-1 ) * v_in
	void solve(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// Calculate v_out = ( M^-1 ) * rhs
	void solve(libMesh::PetscVector<libMesh::Number>& v_out);

	void print_type();


};

}




#endif /* COMMON_LIBMESH_CODE_KSP_LINEAR_SOLVER_H_ */
