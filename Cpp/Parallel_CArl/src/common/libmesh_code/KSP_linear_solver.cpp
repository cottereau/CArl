/*
 * KSP_linear_solver.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "KSP_linear_solver.h"
namespace carl
{

libMesh::PetscLinearSolver<libMesh::Number> * KSP_linear_solver::solver_ptr()
{
	return m_KSP_solver.get();
};

libMesh::PetscMatrix<libMesh::Number> * KSP_linear_solver::matrix_ptr()
{
	return m_Matrix;
}
void KSP_linear_solver::set_solver(libMesh::PetscMatrix<libMesh::Number>& matrix, const std::string name)
{
	m_KSP_solver->reuse_preconditioner(true);
	this->set_matrix(matrix);
	m_KSP_solver->init(m_Matrix, name.c_str());
}

void KSP_linear_solver::set_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix)
{
	m_Matrix = &matrix;
	libMesh::PetscVector<libMesh::Number> * dummy_vec_1 = m_aux_vec_1.get();
	libMesh::PetscVector<libMesh::Number> * dummy_vec_2 = m_aux_vec_2.get();

	PetscInt N, local_N;
	MatGetSize(matrix.mat(),&N,NULL);
	MatGetLocalSize(matrix.mat(),&local_N,NULL);

	dummy_vec_1->init((unsigned int)N,(unsigned int)local_N);
	dummy_vec_2->init((unsigned int)N,(unsigned int)local_N);

	m_bMatrixSet = true;
};

void KSP_linear_solver::set_coupling_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix)
{
	m_Coupling = &matrix;
	m_bCouplingSet = true;
};

void KSP_linear_solver::set_rhs(libMesh::PetscVector<libMesh::Number>& vector)
{
	m_rhs = &vector;
	m_bRhsSet = true;
};

KSPConvergedReason KSP_linear_solver::get_converged_reason()
{
	return m_conv_reason;
}

void KSP_linear_solver::get_system_dimensions(unsigned int& M_out, unsigned int& M_local_out)
{
	homemade_assert_msg(m_bMatrixSet, "Must be called after setting the system matrix!\n");
	M_out = m_Matrix->m();
	int silly_local;
	MatGetLocalSize(m_Matrix->mat(),&silly_local,NULL);
	M_local_out = silly_local;
}

void KSP_linear_solver::get_coupling_dimensions(unsigned int& M_out, unsigned int& N_out, unsigned int& M_local_out, unsigned int& N_local_out)
{
	homemade_assert_msg(m_bMatrixSet, "Must be called after setting the coupling matrix!\n");
	M_out = m_Coupling->m();
	N_out = m_Coupling->n();
	int silly_local_N;
	int silly_local_M;
	MatGetLocalSize(m_Coupling->mat(),&silly_local_M,&silly_local_N);
	M_local_out = silly_local_M;
	N_local_out = silly_local_N;
}

void KSP_linear_solver::apply_ZMiZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	homemade_assert_msg(m_bCouplingSet, "Coupling matrix not set yet!\n");

	libMesh::PetscVector<libMesh::Number> * aux_vec_rhs_ptr = m_aux_vec_1.get();
	libMesh::PetscVector<libMesh::Number> * aux_vec_sol_ptr = m_aux_vec_2.get();

	MatMultTranspose(m_Coupling->mat(),v_in.vec(),aux_vec_rhs_ptr->vec());
	this->solve(*aux_vec_sol_ptr,*aux_vec_rhs_ptr);
	m_Coupling->vector_mult(v_out,*aux_vec_sol_ptr);
}

void KSP_linear_solver::apply_MiZt(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	homemade_assert_msg(m_bCouplingSet, "Coupling matrix not set yet!\n");

	libMesh::PetscVector<libMesh::Number> * aux_vec_rhs_ptr = m_aux_vec_1.get();

	MatMultTranspose(m_Coupling->mat(),v_in.vec(),aux_vec_rhs_ptr->vec());
	this->solve(v_out,*aux_vec_rhs_ptr);
}

// Implementations of virtual methods
void KSP_linear_solver::solve(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	m_KSP_solver->solve(*m_Matrix,v_out,v_in,m_KSP_eps,m_KSP_iter_max);
	KSPGetConvergedReason(m_KSP_solver->ksp(),&m_conv_reason);
}

void KSP_linear_solver::solve(libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_assert_msg(m_bRhsSet, "Solver RHS not set yet!\n");
	homemade_assert_msg(m_bMatrixSet, "Solver matrix not set yet!\n");
	m_KSP_solver->solve(*m_Matrix,v_out,*m_rhs,m_KSP_eps,m_KSP_iter_max);
	KSPGetConvergedReason(m_KSP_solver->ksp(),&m_conv_reason);
}

void KSP_linear_solver::print_type()
{
	KSPType solver_type_string;
	KSPGetType(m_KSP_solver->ksp(),&solver_type_string);
	std::cout 	<< "|        Solver type : " << solver_type_string << std::endl;
	std::cout 	<< "|        PC     type : " << m_KSP_solver->preconditioner_type() << std::endl;
	std::cout   << "|" << std::endl;
}

}


