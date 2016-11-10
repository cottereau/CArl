/*
 * KSP_linear_solver.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "KSP_linear_solver.h"
namespace carl
{

void KSP_linear_solver::set_solver(libMesh::PetscMatrix<libMesh::Number>& matrix, const std::string name)
{
	m_KSP_solver->reuse_preconditioner(true);
	m_Matrix = &matrix;
	m_KSP_solver->init(m_Matrix, name.c_str());
}

void KSP_linear_solver::set_matrix(libMesh::PetscMatrix<libMesh::Number>& matrix)
{
	m_Matrix = &matrix;
};

// Implementations of virtual methods
void KSP_linear_solver::solve(libMesh::PetscVector<libMesh::Number>& sol, libMesh::PetscVector<libMesh::Number>& rhs)
{
	m_KSP_solver->solve(*m_Matrix,sol,rhs,m_KSP_eps,m_KSP_iter_max);
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


