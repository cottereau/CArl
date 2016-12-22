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

void carl::coupled_solver::build_null_space_projection_matrices(libMesh::PetscMatrix<libMesh::Number>& M_sys, libMesh::PetscMatrix<libMesh::Number>& C_sys)
{
	homemade_assert_msg(m_bMatricesSetUp,"Must have the matrices ready!");

	// Get the rigid body modes vector and build the corresponding matrix
	MatNullSpace nullsp_sys;

	MatGetNullSpace(M_sys.mat(),&nullsp_sys);

	if(nullsp_sys)
	{
		PetscBool 	null_has_cte;
		PetscInt  	null_nb_vecs;
		const Vec*	null_vecs;
		PetscInt	C_sys_M, C_sys_N, C_sys_M_local, C_sys_N_local;
		PetscInt	R_mat_M, R_mat_N, R_mat_M_local, R_mat_N_local;
		Mat			RITRI_mat, inv_RITRI_mat;

		// Get the input matrix's dimensions
		MatGetLocalSize(C_sys.mat(),&C_sys_M_local,&C_sys_N_local);
		MatGetSize(C_sys.mat(),&C_sys_M,&C_sys_N);

		// -> Create the matrices!
		// R_mat      : n_sys   x nb_vecs ( 3 for 2D, 6 for 3D )
		MatNullSpaceGetVecs(nullsp_sys,&null_has_cte,&null_nb_vecs,&null_vecs);
		create_PETSC_dense_matrix_from_vectors(null_vecs,null_nb_vecs,m_null_R);
		MatTranspose(m_null_R,MAT_INITIAL_MATRIX,&m_null_RT);

		MatGetLocalSize(m_null_R,&R_mat_M_local,&R_mat_N_local);
		MatGetSize(m_null_R,&R_mat_M,&R_mat_N);

		//             M       x N
		// C_sys     : n_coupl x n_sys
		// RI_mat    : n_coupl x nb_vecs ( 3 for 2D, 6 for 3D )
		// RITRI_mat : nb_vecs x nb_vecs ( 3 for 2D, 6 for 3D )
		MatCreateDense(PETSC_COMM_WORLD,C_sys_M_local,R_mat_N_local,C_sys_M,R_mat_N,NULL,&RI_mat);
		MatCreateDense(PETSC_COMM_WORLD,R_mat_N_local,R_mat_N_local,R_mat_N,R_mat_N,NULL,&RITRI_mat);

		// RI = C_sys * R_mat;
		MatMatMult(C_sys.mat(),m_null_R,MAT_REUSE_MATRIX,PETSC_DECIDE,&RI_mat);
		MatTranspose(RI_mat,MAT_INITIAL_MATRIX,&RI_T_mat);

		// Cannot use MatTransposeMatMult with dense matrices ...
		MatMatMult(RI_T_mat,RI_mat,MAT_REUSE_MATRIX,PETSC_DECIDE,&RITRI_mat);

		// Invert (fortunately, only a 6x6 matrix ...)
		PETSC_invert_dense_matrix(RITRI_mat,inv_RITRI_mat);

		// Calculate the projectors!

		// PI_mat    : n_coupl x n_coupl
		// F_mat     : n_coupl x n_sys		(same as C_sys)
		MatCreateDense(PETSC_COMM_WORLD,C_sys_M_local,C_sys_M_local,C_sys_M,C_sys_M,NULL,&m_null_PI);
		MatCreateDense(PETSC_COMM_WORLD,C_sys_M_local,C_sys_N_local,C_sys_M,C_sys_N,NULL,&m_null_F);

		// PI_mat = Id - RI_mat * ( inv_RITRI_mat ) * RI_T_mat
		// ... but MatMatMatMult is not supported for dense matrices ...

		// aux_matrix = RI_mat * inv_RITRI_mat
		// aux_matrix : n_coupl x nb_vecs (same as RI_mat)
		Mat aux_matrix;
		MatDuplicate(RI_mat,MAT_DO_NOT_COPY_VALUES,&aux_matrix);
		MatMatMult(RI_mat,inv_RITRI_mat,MAT_REUSE_MATRIX,PETSC_DECIDE,&aux_matrix);

		// PI_mat = Id - aux_matrix * RI_T_mat
		MatMatMult(aux_matrix,RI_T_mat,MAT_REUSE_MATRIX,PETSC_DECIDE,&m_null_PI);
		MatScale(m_null_PI,-1);
		MatShift(m_null_PI,1);

		// F_mat = aux_matrix * R_T_mat
		MatMatMult(aux_matrix,m_null_RT,MAT_REUSE_MATRIX,PETSC_DECIDE,&m_null_F);

		// Set up flag
		m_bCreatedRigidBodyProjectors = true;

		// Cleanup
		MatDestroy(&aux_matrix);
//		MatDestroy(&RI_mat);
//		MatDestroy(&RI_T_mat);
		MatDestroy(&RITRI_mat);
	}
};
