/*
 * coupled_solver.cpp
 *
 *  Created on: Nov 3, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "coupled_solver.h"

void carl::coupled_solver::set_sys_names(const std::string& name_A, const std::string& name_B)
{
	m_ksp_name_A = name_A + "_";
	m_ksp_name_B = name_B + "_";
};

void carl::coupled_solver::set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
					libMesh::PetscMatrix<libMesh::Number>& M_B,
					libMesh::PetscMatrix<libMesh::Number>& C_RA,
					libMesh::PetscMatrix<libMesh::Number>& C_RB,
					libMesh::PetscMatrix<libMesh::Number>& C_RR)
{
	m_C_RA = &C_RA;
	m_C_RB = &C_RB;
	m_C_RR = &C_RR;

	m_M_A = &M_A;
	m_M_B = &M_B;

	std::cout << "| K_A " << std::endl;
	print_matrix_dim(*m_M_A);

	std::cout << "| K_B " << std::endl;
	print_matrix_dim(*m_M_B);

	m_bMatricesSetUp = true;
};

void carl::coupled_solver::set_forces(	libMesh::PetscVector<libMesh::Number>& F_A,
					libMesh::PetscVector<libMesh::Number>& F_B)
{
	m_F_A = &F_A;
	m_F_B = &F_B;

	m_bForcesSetUp = true;
};

void carl::coupled_solver::set_restart( 	bool bUseRestart,
												bool bPrintRestart,
												const std::string& restart_base_filename)
{
	m_bUseRestart = bUseRestart;
	m_bPrintRestart = bPrintRestart;
	if(m_bUseRestart || m_bPrintRestart)
	{
		m_sol_A_filename 	= restart_base_filename + "_sol_A.dat";
		m_sol_B_filename 	= restart_base_filename + "_sol_B.dat";
	}
}

void carl::coupled_solver::set_info(	bool bSavePartitionInfo,
											const std::string& info_base_filename)
{
	m_bSavePartitionInfo = bSavePartitionInfo;
	if(m_bSavePartitionInfo)
	{
		m_info_matrix_M_A_filename 	= info_base_filename + "_matrix_M_A.dat";
		m_info_matrix_M_B_filename  = info_base_filename + "_matrix_M_B.dat";
		m_matrix_M_A_filename 	= info_base_filename + "_matrix_M_A.m";
		m_matrix_M_B_filename  = info_base_filename + "_matrix_M_B.m";
		m_matrix_C_A_filename 	= info_base_filename + "_matrix_C_A.m";
		m_matrix_C_B_filename  = info_base_filename + "_matrix_C_B.m";
	}
}

libMesh::NumericVector<libMesh::Number>& carl::coupled_solver::get_solution_micro()
{
	return *m_sol_B.get();
}

libMesh::NumericVector<libMesh::Number>& carl::coupled_solver::get_solution_BIG()
{
	return *m_sol_A.get();
}

void carl::coupled_solver::check_dimensions()
{
	int M_A_mmm, M_A_nnn, M_A_local_mmm, M_A_local_nnn,
		M_B_mmm, M_B_nnn, M_B_local_mmm, M_B_local_nnn,

		C_A_mmm, C_A_nnn, C_A_local_mmm, C_A_local_nnn,
		C_B_mmm, C_B_nnn, C_B_local_mmm, C_B_local_nnn,

		F_A_size, F_A_local_size,
		F_B_size, F_B_local_size;


	M_A_mmm = m_M_A->m(); M_A_nnn = m_M_A->n();
	MatGetLocalSize(m_M_A->mat(),&M_A_local_mmm,&M_A_local_nnn);

	M_B_mmm = m_M_B->m(); M_B_nnn = m_M_B->n();
	MatGetLocalSize(m_M_B->mat(),&M_B_local_mmm,&M_B_local_nnn);

	C_A_mmm = m_C_RA->m(); C_A_nnn = m_C_RA->n();
	MatGetLocalSize(m_C_RA->mat(),&C_A_local_mmm,&C_A_local_nnn);

	C_B_mmm = m_C_RB->m(); C_B_nnn = m_C_RB->n();
	MatGetLocalSize(m_C_RB->mat(),&C_B_local_mmm,&C_B_local_nnn);

	F_A_size = m_F_A->size(); F_A_local_size = m_F_A->local_size();
	F_B_size = m_F_B->size(); F_B_local_size = m_F_B->local_size();

	// Test if the matrices are squared
	homemade_assert_msg( M_A_mmm == M_A_nnn , "   check_dimensions : M_A.m() != M_A.n() !");
	homemade_assert_msg( M_B_mmm == M_B_nnn , "   check_dimensions : M_B.m() != M_B.n() !");

	// Test the couplings
	homemade_assert_msg( M_A_nnn == C_A_nnn , "   check_dimensions : M_A.n() != C_A.n() !");
	homemade_assert_msg( M_B_nnn == C_B_nnn , "   check_dimensions : M_B.n() != C_B.n() !");
	homemade_assert_msg( C_A_mmm == C_B_mmm , "   check_dimensions : C_A.m() != C_B.m() !");

	// Test the forces
	homemade_assert_msg( M_A_nnn == F_A_size , "   check_dimensions : M_A.n() != F_A.size() !");
	homemade_assert_msg( M_B_nnn == F_B_size , "   check_dimensions : M_B.n() != F_B.size() !");


	// -> Now local !

	// Test if the matrices are squared
	homemade_assert_msg( M_A_local_mmm == M_A_local_nnn , "   check_dimensions : M_A.m() != M_A.n() (local) !");
	homemade_assert_msg( M_B_local_mmm == M_B_local_nnn , "   check_dimensions : M_B.m() != M_B.n() (local) !");

	// Test the couplings
	homemade_assert_msg( M_A_local_nnn == C_A_local_nnn , "   check_dimensions : M_A.n() != C_A.n() (local) !");
	homemade_assert_msg( M_B_local_nnn == C_B_local_nnn , "   check_dimensions : M_B.n() != C_B.n() (local) !");
	homemade_assert_msg( C_A_local_mmm == C_B_local_mmm , "   check_dimensions : C_A.m() != C_B.m() (local) !");

	// Test the forces
	homemade_assert_msg( M_A_local_nnn == F_A_local_size , "   check_dimensions : M_A.n() != F_A.size() (local) !");
	homemade_assert_msg( M_B_local_nnn == F_B_local_size , "   check_dimensions : M_B.n() != F_B.size() (local) !");

	m_bCheckDimensions = true;
};
