#include "LATIN_solver.h"

void carl::PETSC_LATIN_solver::use_exact_inverse()
{
	m_bUseLumping = false;
};

void carl::PETSC_LATIN_solver::use_lumped_inverse()
{
	m_bUseLumping = true;
};

void carl::PETSC_LATIN_solver::set_params(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB)
{
	m_k_dA = i_k_dA;
	m_k_dB = i_k_dB;
	m_k_cA = i_k_cA;
	m_k_cB = i_k_cB;
	m_bParamsSetUp = true;
};

void carl::PETSC_LATIN_solver::set_sys_names(const std::string& name_A, const std::string& name_B)
{
	m_ksp_name_A = name_A + "_";
	m_ksp_name_B = name_B + "_";
};

void carl::PETSC_LATIN_solver::set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
					libMesh::PetscMatrix<libMesh::Number>& M_B,
					libMesh::PetscMatrix<libMesh::Number>& C_RA,
					libMesh::PetscMatrix<libMesh::Number>& C_RB,
					libMesh::PetscMatrix<libMesh::Number>& C_RR,
					double product_prealloc_P_A,
					double product_prealloc_P_B,
					double product_prealloc_H_A,
					double product_prealloc_H_B)
{
	libMesh::PerfLog perf_log("Matrix setup",MASTER_bPerfLog_LATIN_solver_matrix_assemble);

	// -> Will need k_d/cI
	homemade_assert_msg( m_bParamsSetUp , "LATIN parameters not set up!");

	m_M_A = &M_A;
	m_M_B = &M_B;

	m_C_RA = &C_RA;
	m_C_RB = &C_RB;
	m_C_RR = &C_RR;

	// -> Invert the C_RR matrix
	if(m_bUseLumping)
	{
		perf_log.push("Lumping");
		lump_matrix_and_invert(* m_C_RR,* m_invC_RR_vec);
		perf_log.pop("Lumping");
	}
	else
	{
		// TODO : implement inverse matrix
		libmesh_error_msg("   set_matrices : Exact inverse matrix not implemented yet!!!");
	}

	/*
	 * 	Since the definition of the PETSc matrix of a libMesh::PetscMatrix
	 *	object is done at the declaration, we must use a dynamically
	 *	allocated libMesh::PetscMatrix, and remember to de-allocate it
	 *	during the destruction.
	 */

	m_bDeallocateMatrices = true;

	// -> Calculate P_I = invC_RR * C_I
	perf_log.push("P_I");
	MatConvert(m_C_RA->mat(),MATSAME,MAT_INITIAL_MATRIX,&m_PETSC_P_A);
	MatConvert(m_C_RB->mat(),MATSAME,MAT_INITIAL_MATRIX,&m_PETSC_P_B);

	MatDiagonalScale(m_PETSC_P_A,m_invC_RR_vec->vec(),PETSC_NULL);
	MatDiagonalScale(m_PETSC_P_B,m_invC_RR_vec->vec(),PETSC_NULL);
//		MatMatMult(m_invC_RR.mat(), m_C_RA->mat(), MAT_INITIAL_MATRIX, product_prealloc_P_A, &m_PETSC_P_A);
//		MatMatMult(m_invC_RR.mat(), m_C_RB->mat(), MAT_INITIAL_MATRIX, product_prealloc_P_B, &m_PETSC_P_B);

	m_P_A = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_P_A, *m_comm);
	m_P_B = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_P_B, *m_comm);
	perf_log.pop("P_I");

	std::cout << "| P_A " << std::endl;
	print_matrix_dim(*m_P_A);

	std::cout << "| P_B " << std::endl;
	print_matrix_dim(*m_P_B);

	// -> Calculate H_I = k_dI * C_I^t * P_I + M_I

	// C_I^t * P_I
	perf_log.push("Product","H_I");
	MatTransposeMatMult(m_C_RA->mat(), m_PETSC_P_A, MAT_INITIAL_MATRIX, product_prealloc_H_A, &m_PETSC_H_A);
	MatTransposeMatMult(m_C_RB->mat(), m_PETSC_P_B, MAT_INITIAL_MATRIX, product_prealloc_H_B, &m_PETSC_H_B);
	perf_log.pop("Product","H_I");

	// k_dI * ( C_I^t * P_I ) + M_I
	perf_log.push("Sum","H_I");
	MatAYPX(m_PETSC_H_A, m_k_dA, m_M_A->mat(), DIFFERENT_NONZERO_PATTERN);
	MatAYPX(m_PETSC_H_B, m_k_dB, m_M_B->mat(), DIFFERENT_NONZERO_PATTERN);
	perf_log.pop("Sum","H_I");

	perf_log.push("Build","H_I");
	m_H_A = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_H_A, *m_comm);
	m_H_B = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_H_B, *m_comm);
	perf_log.pop("Build","H_I");

	std::cout << "| H_A " << std::endl;
	print_matrix_dim(*m_H_A);

	std::cout << "| H_B " << std::endl;
	print_matrix_dim(*m_H_B);

	 std::cout 	<< " L1 : MA " << m_M_A->l1_norm() << " | "
			   	<< " HA " << m_H_A->l1_norm() << " | "
	 			<< " HB " << m_H_B->l1_norm() << " | "
	 			<< " CR " << m_C_RR->l1_norm() << " | "
	 			<< " invCR " << m_invC_RR_vec->l1_norm() << " | "
	 			<< " PA " << m_P_A->l1_norm() << " | "
	 			<< " PB " << m_P_B->l1_norm() << std::endl << std::endl;

	 std::cout 	<< " linfty : MA " << m_M_A->linfty_norm() << " | "
			 	<< " HA " << m_H_A->linfty_norm() << " | "
	 			<< " HB " << m_H_B->linfty_norm() << " | "
	 			<< " CR " << m_C_RR->linfty_norm() << " | "
	 			<< " invCR " << m_invC_RR_vec->linfty_norm() << " | "
	 			<< " PA " << m_P_A->linfty_norm() << " | "
	 			<< " PB " << m_P_B->linfty_norm() << std::endl << std::endl;

	m_bMatricesSetUp = true;
};

void carl::PETSC_LATIN_solver::set_matrices_nonlinear(	libMesh::PetscMatrix<libMesh::Number>& M_A,
					libMesh::PetscMatrix<libMesh::Number>& C_RA,
					libMesh::PetscMatrix<libMesh::Number>& C_RB,
					libMesh::PetscMatrix<libMesh::Number>& C_RR,
					double product_prealloc_P_A,
					double product_prealloc_P_B,
					double product_prealloc_H_A,
					double product_prealloc_H_B)
{
	libMesh::PerfLog perf_log("Matrix setup",MASTER_bPerfLog_LATIN_solver_matrix_assemble);

	// -> Will need k_d/cI
	homemade_assert_msg( m_bParamsSetUp , "LATIN parameters not set up!");

	m_M_A = &M_A;

	m_C_RA = &C_RA;
	m_C_RB = &C_RB;
	m_C_RR = &C_RR;

	// -> Invert the C_RR matrix
	if(m_bUseLumping)
	{
		perf_log.push("Lumping");
		lump_matrix_and_invert(* m_C_RR,* m_invC_RR_vec);
		perf_log.pop("Lumping");
	}
	else
	{
		// TODO : implement inverse matrix
		libmesh_error_msg("   set_matrices : Exact inverse matrix not implemented yet!!!");
	}

	/*
	 * 	Since the definition of the PETSc matrix of a libMesh::PetscMatrix
	 *	object is done at the declaration, we must use a dynamically
	 *	allocated libMesh::PetscMatrix, and remember to de-allocate it
	 *	during the destruction.
	 */

	m_bDeallocateMatrices = true;

	// -> Calculate P_I = invC_RR * C_I
	perf_log.push("P_I");
	MatConvert(m_C_RA->mat(),MATSAME,MAT_INITIAL_MATRIX,&m_PETSC_P_A);
	MatConvert(m_C_RB->mat(),MATSAME,MAT_INITIAL_MATRIX,&m_PETSC_P_B);

	MatDiagonalScale(m_PETSC_P_A,m_invC_RR_vec->vec(),PETSC_NULL);
	MatDiagonalScale(m_PETSC_P_B,m_invC_RR_vec->vec(),PETSC_NULL);
//		MatMatMult(m_invC_RR.mat(), m_C_RA->mat(), MAT_INITIAL_MATRIX, product_prealloc_P_A, &m_PETSC_P_A);
//		MatMatMult(m_invC_RR.mat(), m_C_RB->mat(), MAT_INITIAL_MATRIX, product_prealloc_P_B, &m_PETSC_P_B);

	m_P_A = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_P_A, *m_comm);
	m_P_B = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_P_B, *m_comm);
	perf_log.pop("P_I");

	std::cout << "| P_A " << std::endl;
	print_matrix_dim(*m_P_A);

	std::cout << "| P_B " << std::endl;
	print_matrix_dim(*m_P_B);

	// -> Calculate H_I = k_dI * C_I^t * P_I + M_I

	// C_I^t * P_I
	perf_log.push("Product","H_I");
	MatTransposeMatMult(m_C_RA->mat(), m_PETSC_P_A, MAT_INITIAL_MATRIX, product_prealloc_H_A, &m_PETSC_H_A);
	MatTransposeMatMult(m_C_RB->mat(), m_PETSC_P_B, MAT_INITIAL_MATRIX, product_prealloc_H_B, &m_PETSC_Extra_M_B);
	perf_log.pop("Product","H_I");

	// k_dI * ( C_I^t * P_I ) + M_I
	perf_log.push("Sum","H_I");
	MatAYPX(m_PETSC_H_A, m_k_dA, m_M_A->mat(), DIFFERENT_NONZERO_PATTERN);
	// Extra = k_dI * ( C_I^t * P_I )
	MatScale(m_PETSC_Extra_M_B, m_k_dB);
	perf_log.pop("Sum","H_I");

	perf_log.push("Build","H_I");
	m_H_A = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_H_A, *m_comm);
	m_Extra_M_B = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_Extra_M_B, *m_comm);
	perf_log.pop("Build","H_I");

	std::cout << "| H_A " << std::endl;
	print_matrix_dim(*m_H_A);

	std::cout << "| E_B " << std::endl;
	print_matrix_dim(*m_Extra_M_B);

	 std::cout 	<< " L1 : MA " << m_M_A->l1_norm() << " | "
			   	<< " HA " << m_H_A->l1_norm() << " | "
	 			<< " EB " << m_Extra_M_B->l1_norm() << " | "
	 			<< " CR " << m_C_RR->l1_norm() << " | "
	 			<< " invCR " << m_invC_RR_vec->l1_norm() << " | "
	 			<< " PA " << m_P_A->l1_norm() << " | "
	 			<< " PB " << m_P_B->l1_norm() << std::endl << std::endl;

	 std::cout 	<< " linfty : MA " << m_M_A->linfty_norm() << " | "
			 	<< " HA " << m_H_A->linfty_norm() << " | "
	 			<< " EB " << m_Extra_M_B->linfty_norm() << " | "
	 			<< " CR " << m_C_RR->linfty_norm() << " | "
	 			<< " invCR " << m_invC_RR_vec->linfty_norm() << " | "
	 			<< " PA " << m_P_A->linfty_norm() << " | "
	 			<< " PB " << m_P_B->linfty_norm() << std::endl << std::endl;

	m_bMatricesSetUp = true;
};

void carl::PETSC_LATIN_solver::set_forces(	libMesh::PetscVector<libMesh::Number>& F_A,
					libMesh::PetscVector<libMesh::Number>& F_B)
{
	m_F_A = &F_A;
	m_F_B = &F_B;

	m_bForcesSetUp = true;
};

void carl::PETSC_LATIN_solver::set_forces_nonlinear(	libMesh::PetscVector<libMesh::Number>& F_A)
{
	m_F_A = &F_A;

	m_bForcesSetUp = true;
};

void carl::PETSC_LATIN_solver::set_convergence_limits(double eps, int convIter)
{
	m_LATIN_conv_eps = eps;
	m_LATIN_conv_max_n = convIter;
};

void carl::PETSC_LATIN_solver::set_relaxation(double relax)
{
	m_LATIN_relax = relax;
}

void carl::PETSC_LATIN_solver::solve()
{
	libMesh::PerfLog perf_log("Solve",MASTER_bPerfLog_LATIN_solver_solve);

	std::cout << "| LATIN solver: " << std::endl;
	std::cout << "|     Initialization ..." << std::endl; std::cout.flush();

	// -> Test if the parameters are set up
	homemade_assert_msg( m_bParamsSetUp , "   solve : LATIN parameters not set up!");
	homemade_assert_msg( m_bMatricesSetUp , "   solve : Matrices not set up!");
	homemade_assert_msg( m_bForcesSetUp , "   solve : Forces not set up!");

	// -> Matrix dimensions

	// Check them beforehand
	this->check_dimensions();

	unsigned int 	dim_A, dim_A_local,
					dim_B, dim_B_local,
					dim_R, dim_R_local;

	int silly_local = 0;

	dim_A = m_H_A->m();
	dim_B = m_H_B->m();
	dim_R = m_P_A->m();

	MatGetLocalSize(m_H_A->mat(),&silly_local,NULL);
	dim_A_local = silly_local;
	MatGetLocalSize(m_H_B->mat(),&silly_local,NULL);
	dim_B_local = silly_local;
	MatGetLocalSize(m_P_A->mat(),&silly_local,NULL);
	dim_R_local = silly_local;

	// -> Declarations
	double iter_eps = 1 + m_LATIN_conv_eps;
	int    iter_nb = 0;

	double norm_diff = 0;
	double norm_sum = 0;

	perf_log.push("Vector declarations","Initialization");
	m_sol_A->init(dim_A,dim_A_local);
	m_sol_B->init(dim_B,dim_B_local);

	// Temporary solutions
	libMesh::PetscVector<libMesh::Number> w_A(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> w_B(*m_comm, dim_R, dim_R_local);
	w_A.zero();
	w_B.zero();

	libMesh::PetscVector<libMesh::Number> u_A(*m_comm, dim_A, dim_A_local);
	libMesh::PetscVector<libMesh::Number> u_B(*m_comm, dim_B, dim_B_local);
	u_A.zero();
	u_B.zero();

	// Effective forces
	libMesh::PetscVector<libMesh::Number> f_eff_A(*m_comm, dim_A, dim_A_local);
	libMesh::PetscVector<libMesh::Number> f_eff_B(*m_comm, dim_B, dim_B_local);
	f_eff_A.zero();
	f_eff_B.zero();

	// Auxiliary vectors
	libMesh::PetscVector<libMesh::Number> aux_A(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> aux_B(*m_comm, dim_R, dim_R_local);
	aux_A.zero();
	aux_B.zero();

	libMesh::PetscVector<libMesh::Number> phi_cA(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> phi_cB(*m_comm, dim_R, dim_R_local);
	phi_cA.zero();
	phi_cB.zero();

	libMesh::PetscVector<libMesh::Number> phi_dA(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> phi_dB(*m_comm, dim_R, dim_R_local);
	phi_dA.zero();
	phi_dB.zero();

	libMesh::PetscVector<libMesh::Number> phi_diff_A(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> phi_diff_B(*m_comm, dim_R, dim_R_local);
	phi_diff_A.zero();
	phi_diff_B.zero();

	libMesh::PetscVector<libMesh::Number> phi_sum_A(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> phi_sum_B(*m_comm, dim_R, dim_R_local);
	phi_sum_A.zero();
	phi_sum_B.zero();

	perf_log.pop("Vector declarations","Initialization");

	// Solvers
	perf_log.push("KSP solvers setup","Initialization");
	libMesh::PetscLinearSolver<libMesh::Number> KSP_Solver_H_A(*m_comm);
	libMesh::PetscLinearSolver<libMesh::Number> KSP_Solver_H_B(*m_comm);
	KSP_Solver_H_A.reuse_preconditioner(true);
	KSP_Solver_H_B.reuse_preconditioner(true);

	KSP_Solver_H_A.init(m_H_A, m_ksp_name_A.c_str());
	KSP_Solver_H_B.init(m_H_B, m_ksp_name_B.c_str());
	perf_log.pop("KSP solvers setup","Initialization");

	// -> Initialize the vectors
	// u_0,I = H_I^-1 * F_I (KSP SOLVER!)

	perf_log.push("KSP solver - A","Initialization");
	KSP_Solver_H_A.solve ( *m_H_A, *m_sol_A, *m_F_A, m_KSP_A_eps, m_KSP_A_iter_max);
	perf_log.pop("KSP solver - A","Initialization");
	perf_log.push("KSP solver - B","Initialization");
	KSP_Solver_H_B.solve ( *m_H_B, *m_sol_B, *m_F_B, m_KSP_B_eps, m_KSP_B_iter_max);
	perf_log.pop("KSP solver - B","Initialization");

	// w_0,I = P_I * u_0,I
	m_P_A->vector_mult(w_A,*m_sol_A);
	m_P_B->vector_mult(w_B,*m_sol_B);

	// phi_d0,I = - k_dI * w_0,I
	phi_dA.add(-m_k_dA,w_A);
	phi_dB.add(-m_k_dB,w_B);

	while (iter_eps > m_LATIN_conv_eps && iter_nb < m_LATIN_conv_max_n)
	{
//		clear_line();
		std::cout << "|     Iter no. " << iter_nb << " " << iter_eps << std::endl; std::cout.flush();
		// -> Coupled step
		perf_log.push("Coupled iterations");

		// aux_i,I = k_cI * w_i-1,I - phi_i-1,I
		aux_A = phi_dA;
		aux_A.scale(-1);
		aux_A.add(m_k_cA,w_A);

		aux_B = phi_dB;
		aux_B.scale(-1);
		aux_B.add(m_k_cB,w_B);

		// w_i,I = ( aux_i,A + aux_i,B ) / ( k_cA + k_cB )
		w_A = aux_A;
		w_A.add(aux_B);
		w_A.scale( 1./(m_k_cA + m_k_cB) );
		w_B = w_A;

		// phi_ci,I = k_cI * w_i,I - aux_i,I
		phi_cA = aux_A;
		phi_cA.scale(-1);
		phi_cA.add(m_k_cA,w_A);

		phi_cB = aux_B;
		phi_cB.scale(-1);
		phi_cB.add(m_k_cB,w_B);

		perf_log.pop("Coupled - iterations");

		// -> Decoupled step

		// aux_i,I = - k_dI * w_i,I - phi_i,I
		aux_A = phi_cA;
		aux_A.scale(-1);
		aux_A.add( - m_k_dA,w_A);

		aux_B = phi_cB;
		aux_B.scale(-1);
		aux_B.add( - m_k_dB,w_B);

		perf_log.push("Effective force","Decoupled - iterations");
		// f_eff_i,I = F_A - C_I^t aux_i,I
		MatMultTranspose(m_C_RA->mat(), aux_A.vec(), f_eff_A.vec());
		f_eff_A.scale(-1);
		f_eff_A.add(*m_F_A);
		MatMultTranspose(m_C_RB->mat(), aux_B.vec(), f_eff_B.vec());
		f_eff_B.scale(-1);
		f_eff_B.add(*m_F_B);
		perf_log.pop("Effective force","Decoupled - iterations");

		// u_i,I = H_I^-1 * f_eff_i,I (KSP SOLVER!)
		perf_log.push("KSP solver - A","Decoupled - iterations");
		KSP_Solver_H_A.solve ( *m_H_A, u_A, f_eff_A, m_KSP_A_eps, m_KSP_A_iter_max);
		perf_log.pop("KSP solver - A","Decoupled - iterations");
		perf_log.push("KSP solver - B","Decoupled - iterations");
		KSP_Solver_H_B.solve ( *m_H_B, u_B, f_eff_B, m_KSP_B_eps, m_KSP_B_iter_max);
		perf_log.pop("KSP solver - A","Decoupled - iterations");

		// w_i,I = P_I * u_i,I
		perf_log.push("Projection","Decoupled - iterations");
		m_P_A->vector_mult(w_A,u_A);
		m_P_B->vector_mult(w_B,u_B);
		perf_log.pop("Projection","Decoupled - iterations");

		// phi_i,I = - k_dI * w_i,I - aux_i,I
		phi_dA = aux_A;
		phi_dA.scale(-1);
		phi_dA.add( - m_k_dA,w_A);

		phi_dB = aux_B;
		phi_dB.scale(-1);
		phi_dB.add( - m_k_dB,w_B);

		// Relaxation : sol_i,I = LATINRelax*u_i,I + (1 - LATINRelax)*sol_i-1,I;
		m_sol_A->scale(1 - m_LATIN_relax);
		m_sol_A->add(m_LATIN_relax,u_A);

		m_sol_B->scale(1 - m_LATIN_relax);
		m_sol_B->add(m_LATIN_relax,u_B);

		// -> Test convergence over the phi's
		phi_diff_A = phi_dA;
		phi_diff_A.add(-1,phi_cA);

		phi_diff_B = phi_dB;
		phi_diff_B.add(-1,phi_cB);

		phi_sum_A = phi_dA;
		phi_sum_A.add(phi_cA);

		phi_sum_B = phi_dB;
		phi_sum_B.add(phi_cB);

		norm_diff = phi_diff_A.l2_norm() + phi_diff_B.l2_norm();
		norm_sum = phi_sum_A.l2_norm() + phi_sum_B.l2_norm();

		if(norm_sum < 1E-18)
		{
			break;
		}

		iter_eps = 2*norm_diff / norm_sum;
		m_LATIN_Index[iter_nb] = iter_eps;
		++iter_nb;
	}

	m_LATIN_conv_n = iter_nb;
//	clear_line();
	std::cout << "|     nb. of iterations : " << iter_nb;
	if(iter_nb == m_LATIN_conv_max_n)
	{
		std::cout << " (MAX!)";
	}
	std::cout << std::endl;
	std::cout << "|     eps               : " << iter_eps << std::endl << std::endl;

	m_bSolved = true;

}

void carl::PETSC_LATIN_solver::solve_nonlinear(libMesh::EquationSystems& EqSys_micro, const std::string type_name_micro)
{
	libMesh::PerfLog perf_log("Solve",MASTER_bPerfLog_LATIN_solver_solve);

	std::cout << "| LATIN solver (nonlinear): " << std::endl;
	std::cout << "|     Initialization ..." << std::endl; std::cout.flush();

	// -> Test if the parameters are set up
	homemade_assert_msg( m_bParamsSetUp , "   solve : LATIN parameters not set up!");
	homemade_assert_msg( m_bMatricesSetUp , "   solve : Matrices not set up!");
	homemade_assert_msg( m_bForcesSetUp , "   solve : Forces not set up!");

	// -> Matrix dimensions

	// Check them beforehand
//	this->check_dimensions();

	unsigned int 	dim_A, dim_A_local,
					dim_B, dim_B_local,
					dim_R, dim_R_local;

	int silly_local = 0;

	dim_A = m_H_A->m();
	std::cout << dim_A << std::endl;
	dim_B = m_Extra_M_B->m();

	std::cout << dim_B << std::endl;
	dim_R = m_P_A->m();

	MatGetLocalSize(m_H_A->mat(),&silly_local,NULL);
	dim_A_local = silly_local;
	MatGetLocalSize(m_Extra_M_B->mat(),&silly_local,NULL);
	dim_B_local = silly_local;
	MatGetLocalSize(m_P_A->mat(),&silly_local,NULL);
	dim_R_local = silly_local;

	// -> Declarations
	double iter_eps = 1 + m_LATIN_conv_eps;
	int    iter_nb = 0;

	double norm_diff = 0;
	double norm_sum = 0;

	perf_log.push("Vector declarations","Initialization");
	m_sol_A->init(dim_A,dim_A_local);
	m_sol_B->init(dim_B,dim_B_local);
	libMesh::PetscVector<libMesh::Number> Extra_F_B(*m_comm,dim_B,dim_B_local);
	Extra_F_B.zero();

	// Temporary solutions
	libMesh::PetscVector<libMesh::Number> w_A(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> w_B(*m_comm, dim_R, dim_R_local);
	w_A.zero();
	w_B.zero();

	libMesh::PetscVector<libMesh::Number> u_A(*m_comm, dim_A, dim_A_local);
	libMesh::PetscVector<libMesh::Number> u_B(*m_comm, dim_B, dim_B_local);
	u_A.zero();
	u_B.zero();

	// Effective forces
	libMesh::PetscVector<libMesh::Number> f_eff_A(*m_comm, dim_A, dim_A_local);
	libMesh::PetscVector<libMesh::Number> f_eff_B(*m_comm, dim_B, dim_B_local);
	f_eff_A.zero();
	f_eff_B.zero();

	// Auxiliary vectors
	libMesh::PetscVector<libMesh::Number> aux_A(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> aux_B(*m_comm, dim_R, dim_R_local);
	aux_A.zero();
	aux_B.zero();

	libMesh::PetscVector<libMesh::Number> phi_cA(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> phi_cB(*m_comm, dim_R, dim_R_local);
	phi_cA.zero();
	phi_cB.zero();

	libMesh::PetscVector<libMesh::Number> phi_dA(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> phi_dB(*m_comm, dim_R, dim_R_local);
	phi_dA.zero();
	phi_dB.zero();

	libMesh::PetscVector<libMesh::Number> phi_diff_A(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> phi_diff_B(*m_comm, dim_R, dim_R_local);
	phi_diff_A.zero();
	phi_diff_B.zero();

	libMesh::PetscVector<libMesh::Number> phi_sum_A(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> phi_sum_B(*m_comm, dim_R, dim_R_local);
	phi_sum_A.zero();
	phi_sum_B.zero();

	perf_log.pop("Vector declarations","Initialization");

	// Solvers
	perf_log.push("KSP solvers setup","Initialization");
	libMesh::PetscLinearSolver<libMesh::Number> KSP_Solver_H_A(*m_comm);
	KSP_Solver_H_A.reuse_preconditioner(true);

	KSP_Solver_H_A.init(m_H_A, m_ksp_name_A.c_str());
	perf_log.pop("KSP solvers setup","Initialization");

	// -> Initialize the vectors
	// u_0,I = H_I^-1 * F_I (KSP SOLVER!)

	perf_log.push("KSP solver","Initialization");
	KSP_Solver_H_A.solve ( *m_H_A, *m_sol_A, *m_F_A, m_KSP_A_eps, m_KSP_A_iter_max);
	perf_log.pop("KSP solver","Initialization");

	// Set up non-linear system
	libMesh::NonlinearImplicitSystem& Sys_micro = libMesh::cast_ref<libMesh::NonlinearImplicitSystem&>(EqSys_micro.get_system(type_name_micro));

	/*
	 * 		What do I have to do here:
	 * 		- Extract the force multiplier and number of solves
	 * 		- For loop:
	 * 			- Assemble the residual and jacobian
	 * 			- Add the coupling term to the residual and the jacobian
	 * 			- Solve
	 */

	int n_solves = EqSys_micro.parameters.get<int>("nb_nonlinear_steps");
	double coupling_residual_multiplier = 1./(n_solves);

	std::cout << EqSys_micro.parameters.get<unsigned int> ("nonlinear solver maximum iterations") << " "
			  << EqSys_micro.parameters.get<libMesh::Real>         ("nonlinear solver absolute residual tolerance") << " "
			  << EqSys_micro.parameters.get<libMesh::Real>         ("nonlinear solver relative residual tolerance") << " " << std::endl;

	LargeDeformationElasticity * micro_assemble_obj = libMesh::cast_ptr<LargeDeformationElasticity *>
			(Sys_micro.nonlinear_solver->residual_object);

	// Associate the extra terms to the residual and jacobian assemble object
	// J_extra = C_extra
	// r_extra : depends on the previous solution
	micro_assemble_obj->set_coupling_term_matrix(m_Extra_M_B);
	micro_assemble_obj->set_coupling_term_vector(&Extra_F_B);

	libMesh::PetscNonlinearSolver<libMesh::Real> * dummySolver =
			libMesh::cast_ptr<libMesh::PetscNonlinearSolver<libMesh::Real> * >(Sys_micro.nonlinear_solver.get());

	KSP nonlinear_ksp;

	SNESGetKSP(dummySolver->snes(),&nonlinear_ksp);

//	for(int nnn = 1; nnn <= n_solves; ++nnn)
	{
		// - Add the residual term
		// r = r_0 - C_extra * u_(n-1)
		// r_extra = - C_extra * u_(n-1)
		m_Extra_M_B->vector_mult(Extra_F_B, *Sys_micro.solution.get());
		Extra_F_B.scale(-1);

		// - Solve
//		Sys_micro.assembly(true,true);
//		std::cout << Sys_micro.matrix->l1_norm() << " " << Sys_micro.rhs->l1_norm() << std::endl;
		Sys_micro.solve();

		Sys_micro.nonlinear_solver->print_converged_reason();

		SNESConvergedReason whyyyyyy;
		SNESGetConvergedReason(dummySolver->snes(),&whyyyyyy);

	    std::cout << "System solved at nonlinear iteration "
	      << Sys_micro.n_nonlinear_iterations()
	      << " , final nonlinear residual norm: "
	      << Sys_micro.final_nonlinear_residual()
	      << std::endl << std::endl;
	}

	m_sol_B = Sys_micro.solution.get();

	// w_0,I = P_I * u_0,I
	m_P_A->vector_mult(w_A,*m_sol_A);
	m_P_B->vector_mult(w_B,*m_sol_B);

	// phi_d0,I = - k_dI * w_0,I
	phi_dA.add(-m_k_dA,w_A);
	phi_dB.add(-m_k_dB,w_B);

	std::cout << m_LATIN_conv_eps << " " << m_LATIN_conv_max_n << std::endl;
	while (iter_eps > m_LATIN_conv_eps && iter_nb < m_LATIN_conv_max_n)
	{
//		clear_line();
		std::cout << "|     Iter no. " << iter_nb << " " << iter_eps << std::endl; std::cout.flush();
		// -> Coupled step
		perf_log.push("Coupled iterations");

		// aux_i,I = k_cI * w_i-1,I - phi_i-1,I
		aux_A = phi_dA;
		aux_A.scale(-1);
		aux_A.add(m_k_cA,w_A);

		aux_B = phi_dB;
		aux_B.scale(-1);
		aux_B.add(m_k_cB,w_B);

		// w_i,I = ( aux_i,A + aux_i,B ) / ( k_cA + k_cB )
		w_A = aux_A;
		w_A.add(aux_B);
		w_A.scale( 1./(m_k_cA + m_k_cB) );
		w_B = w_A;

		// phi_ci,I = k_cI * w_i,I - aux_i,I
		phi_cA = aux_A;
		phi_cA.scale(-1);
		phi_cA.add(m_k_cA,w_A);

		phi_cB = aux_B;
		phi_cB.scale(-1);
		phi_cB.add(m_k_cB,w_B);

		perf_log.pop("Coupled - iterations");

		// -> Decoupled step

		// aux_i,I = - k_dI * w_i,I - phi_i,I
		aux_A = phi_cA;
		aux_A.scale(-1);
		aux_A.add( - m_k_dA,w_A);

		aux_B = phi_cB;
		aux_B.scale(-1);
		aux_B.add( - m_k_dB,w_B);

		perf_log.push("Effective force","Decoupled - iterations");
		// f_eff_A,I = F_A - C_I^t aux_B,I
		MatMultTranspose(m_C_RA->mat(), aux_A.vec(), f_eff_A.vec());
		f_eff_A.scale(-1);
		f_eff_A.add(*m_F_A);
		// f_eff_B,I = C_I^t aux_B,I
		MatMultTranspose(m_C_RB->mat(), aux_B.vec(), f_eff_B.vec());
		perf_log.pop("Effective force","Decoupled - iterations");

		perf_log.push("KSP solver","Decoupled - iterations");
		// u_i,I = H_I^-1 * f_eff_i,I (KSP SOLVER!)
		KSP_Solver_H_A.solve ( *m_H_A, u_A, f_eff_A, m_KSP_A_eps, m_KSP_A_iter_max);
		perf_log.pop("KSP solver","Decoupled - iterations");

		// Non-linear solver
		for(int nnn = 1; nnn <= n_solves; ++nnn)
		{
			// - Add the residual term

			// r = r_0 - C_extra * u_(n-1) - step * f_eff_B
			// r_extra = - C_extra * u_(n-1) - step * f_eff_B
			m_Extra_M_B->vector_mult(Extra_F_B, *Sys_micro.solution.get());
			Extra_F_B.add(nnn*coupling_residual_multiplier,f_eff_B);
			Extra_F_B.scale(-1.);

//			std::cout << Extra_F_B.l1_norm() << " " << Extra_F_B.l2_norm() << std::endl;
//			Sys_micro.assembly(true,true);
//			std::cout << Sys_micro.matrix->l1_norm() << " " << std::endl;
//			std::cout << Sys_micro.rhs->l1_norm() << " " << Sys_micro.rhs->l2_norm()  << std::endl;
			// - Solve
			Sys_micro.nonlinear_solver->init();
			Sys_micro.solve();

			m_comm->barrier();

//			SNESConvergedReason whyyyyyy;
//			SNESGetConvergedReason(dummySolver->snes(),&whyyyyyy);
//			KSPConvergedReason whaaaa;
//			KSPGetConvergedReason(nonlinear_ksp,&whaaaa);
//			int iterations_why = 0;
//			KSPGetIterationNumber(nonlinear_ksp,&iterations_why);
//			std::cout << whyyyyyy << " " << whaaaa << " "<< iterations_why << std::endl;

			Sys_micro.nonlinear_solver->print_converged_reason();

			std::cout << Sys_micro.nonlinear_solver->get_total_linear_iterations() << std::endl;
		    std::cout << "System solved at nonlinear iteration "
		      << Sys_micro.n_nonlinear_iterations()
		      << " , final nonlinear residual norm: "
		      << Sys_micro.final_nonlinear_residual()
		      << std::endl << std::endl;
		}

		// w_i,I = P_I * u_i,I
		perf_log.push("Projection","Decoupled - iterations");
		m_P_A->vector_mult(w_A,u_A);
		m_P_B->vector_mult(w_B,*Sys_micro.solution.get());
		perf_log.pop("Projection","Decoupled - iterations");

		// phi_i,I = - k_dI * w_i,I - aux_i,I
		phi_dA = aux_A;
		phi_dA.scale(-1);
		phi_dA.add( - m_k_dA,w_A);

		phi_dB = aux_B;
		phi_dB.scale(-1);
		phi_dB.add( - m_k_dB,w_B);

		// Relaxation : sol_i,I = LATINRelax*u_i,I + (1 - LATINRelax)*sol_i-1,I;
		m_sol_A->scale(1 - m_LATIN_relax);
		m_sol_A->add(m_LATIN_relax,u_A);

		m_sol_B->scale(1 - m_LATIN_relax);
		m_sol_B->add(m_LATIN_relax,u_B);

		// -> Test convergence over the phi's
		phi_diff_A = phi_dA;
		phi_diff_A.add(-1,phi_cA);

		phi_diff_B = phi_dB;
		phi_diff_B.add(-1,phi_cB);

		phi_sum_A = phi_dA;
		phi_sum_A.add(phi_cA);

		phi_sum_B = phi_dB;
		phi_sum_B.add(phi_cB);

		norm_diff = phi_diff_A.l2_norm() + phi_diff_B.l2_norm();
		norm_sum = phi_sum_A.l2_norm() + phi_sum_B.l2_norm();

		std::cout << norm_sum << std::endl; std::cout.flush();

		if(norm_sum < 1E-18)
		{
			break;
		}

		iter_eps = 2*norm_diff / norm_sum;
		m_LATIN_Index[iter_nb] = iter_eps;
		++iter_nb;
	}

	m_LATIN_conv_n = iter_nb;
//	clear_line();
	std::cout << "|     nb. of iterations : " << iter_nb;
	if(iter_nb == m_LATIN_conv_max_n)
	{
		std::cout << " (MAX!)";
	}
	std::cout << std::endl;
	std::cout << "|     eps               : " << iter_eps << std::endl << std::endl;

	m_bSolved = true;

}

void carl::PETSC_LATIN_solver::check_dimensions()
{
	int H_A_mmm, H_A_nnn, H_A_local_mmm, H_A_local_nnn,
		H_B_mmm, H_B_nnn, H_B_local_mmm, H_B_local_nnn,

		C_A_mmm, C_A_nnn, C_A_local_mmm, C_A_local_nnn,
		C_B_mmm, C_B_nnn, C_B_local_mmm, C_B_local_nnn,

		P_A_mmm, P_A_nnn, P_A_local_mmm, P_A_local_nnn,
		P_B_mmm, P_B_nnn, P_B_local_mmm, P_B_local_nnn,

		F_A_size, F_A_local_size,
		F_B_size, F_B_local_size;


	H_A_mmm = m_H_A->m(); H_A_nnn = m_H_A->n();
	MatGetLocalSize(m_H_A->mat(),&H_A_local_mmm,&H_A_local_nnn);

	H_B_mmm = m_H_B->m(); H_B_nnn = m_H_B->n();
	MatGetLocalSize(m_H_B->mat(),&H_B_local_mmm,&H_B_local_nnn);

	C_A_mmm = m_C_RA->m(); C_A_nnn = m_C_RA->n();
	MatGetLocalSize(m_C_RA->mat(),&C_A_local_mmm,&C_A_local_nnn);

	C_B_mmm = m_C_RB->m(); C_B_nnn = m_C_RB->n();
	MatGetLocalSize(m_C_RB->mat(),&C_B_local_mmm,&C_B_local_nnn);

	P_A_mmm = m_P_A->m(); P_A_nnn = m_P_A->n();
	MatGetLocalSize(m_P_A->mat(),&P_A_local_mmm,&P_A_local_nnn);

	P_B_mmm = m_P_B->m(); P_B_nnn = m_P_B->n();
	MatGetLocalSize(m_P_B->mat(),&P_B_local_mmm,&P_B_local_nnn);

	F_A_size = m_F_A->size(); F_A_local_size = m_F_A->local_size();
	F_B_size = m_F_B->size(); F_B_local_size = m_F_B->local_size();

	// Test if the matrices are squared
	homemade_assert_msg( H_A_mmm == H_A_nnn , "   check_dimensions : H_A.m() != H_A.n() !");
	homemade_assert_msg( H_B_mmm == H_B_nnn , "   check_dimensions : H_B.m() != H_B.n() !");

	// Test the projections
	homemade_assert_msg( H_A_nnn == P_A_nnn , "   check_dimensions : H_A.n() != P_A.n() !");
	homemade_assert_msg( H_B_nnn == P_B_nnn , "   check_dimensions : H_B.n() != P_B.n() !");
	homemade_assert_msg( P_A_mmm == P_B_mmm , "   check_dimensions : P_A.m() != P_B.m() !");

	// Test the couplings
	homemade_assert_msg( H_A_nnn == C_A_nnn , "   check_dimensions : H_A.n() != C_A.n() !");
	homemade_assert_msg( H_B_nnn == C_B_nnn , "   check_dimensions : H_B.n() != C_B.n() !");
	homemade_assert_msg( C_A_mmm == C_B_mmm , "   check_dimensions : C_A.m() != C_B.m() !");

	// Test the forces
	homemade_assert_msg( H_A_nnn == F_A_size , "   check_dimensions : H_A.n() != F_A.size() !");
	homemade_assert_msg( H_B_nnn == F_B_size , "   check_dimensions : H_B.n() != F_B.size() !");


	// -> Now local !

	// Test if the matrices are squared
	homemade_assert_msg( H_A_local_mmm == H_A_local_nnn , "   check_dimensions : H_A.m() != H_A.n() (local) !");
	homemade_assert_msg( H_B_local_mmm == H_B_local_nnn , "   check_dimensions : H_B.m() != H_B.n() (local) !");

	// Test the projections
	homemade_assert_msg( H_A_local_nnn == P_A_local_nnn , "   check_dimensions : H_A.n() != P_A.n() (local) !");
	homemade_assert_msg( H_B_local_nnn == P_B_local_nnn , "   check_dimensions : H_B.n() != P_B.n() (local) !");
	homemade_assert_msg( P_A_local_mmm == P_B_local_mmm , "   check_dimensions : P_A.m() != P_B.m() (local) !");

	// Test the couplings
	homemade_assert_msg( H_A_local_nnn == C_A_local_nnn , "   check_dimensions : H_A.n() != C_A.n() (local) !");
	homemade_assert_msg( H_B_local_nnn == C_B_local_nnn , "   check_dimensions : H_B.n() != C_B.n() (local) !");
	homemade_assert_msg( C_A_local_mmm == C_B_local_mmm , "   check_dimensions : C_A.m() != C_B.m() (local) !");

	// Test the forces
	homemade_assert_msg( H_A_local_nnn == F_A_local_size , "   check_dimensions : H_A.n() != F_A.size() (local) !");
	homemade_assert_msg( H_B_local_nnn == F_B_local_size , "   check_dimensions : H_B.n() != F_B.size() (local) !");

	m_bCheckDimensions = true;
}

libMesh::NumericVector<libMesh::Number>& carl::PETSC_LATIN_solver::get_solution_micro()
{
	return *m_sol_B;
}

libMesh::NumericVector<libMesh::Number>& carl::PETSC_LATIN_solver::get_solution_BIG()
{
	return *m_sol_A;
}


void carl::PETSC_LATIN_solver::print_convergence(std::ostream& convergenceOut)
{
	if(m_bSolved)
	{
		for(int iii = 0; iii < m_LATIN_conv_n; ++iii)
		{
			convergenceOut << iii << " " << m_LATIN_Index[iii] << std::endl;
		}
	}
}
