#include <petscksp.h>
#include <petsc/private/kspimpl.h>
#include "CG_coupled_solver.h"

void carl::PETSC_CG_coupled_solver::set_preconditioner_type(BaseCGPrecondType type_input)
{
	m_precond_type = type_input;
};

void carl::PETSC_CG_coupled_solver::set_info(	bool bSavePartitionInfo,
											const std::string& info_base_filename)
{
	// Call the coupled_solver function
	coupled_solver::set_info(bSavePartitionInfo,info_base_filename);

	if(m_bSavePartitionInfo)
	{
		m_info_matrix_PC_filename 	= info_base_filename + "_matrix_PC.dat";
		m_matrix_PC_filename 	= info_base_filename + "_matrix_PC.m";
	}
}

void carl::PETSC_CG_coupled_solver::set_solvers(generic_solver_interface * solver_A, generic_solver_interface * solver_B)
{
	// First, check if the matrices were set
	homemade_assert_msg(m_bMatricesSetUp, "Must set up the matrices beforehand!\n");

	// Set up the system solvers
	m_sys_A_solver = solver_A;
	m_sys_B_solver = solver_B;

	m_sys_A_solver->set_solver(*m_M_A,m_ksp_name_A.c_str());
	m_sys_A_solver->set_coupling_matrix(*m_C_RA);

	m_sys_B_solver->set_solver(*m_M_B,m_ksp_name_B.c_str());
	m_sys_B_solver->set_coupling_matrix(*m_C_RB);

	m_coupling_solver.set_preconditioner_type(m_precond_type);

	// Set up the preconditioner - if needed
	if(m_precond_type == BaseCGPrecondType::COUPLING_OPERATOR)
	{
		m_coupling_solver.set_inverse_precond(*m_C_RR);
	}

	// Set up the coupled equation solver
	m_coupling_solver.set_solver_CG(*m_sys_A_solver,*m_sys_B_solver);

	m_coupling_solver.set_convergence_limits(m_CG_conv_eps_abs,m_CG_conv_eps_rel,m_CG_conv_max_n,m_CG_div_tol);
};

void carl::PETSC_CG_coupled_solver::set_restart( 	bool bUseRestart,
												bool bPrintRestart,
												const std::string& restart_base_filename)
{
	// Call the coupled_solver function
	coupled_solver::set_restart(bUseRestart,bPrintRestart,restart_base_filename);

	if(m_bUseRestart || m_bPrintRestart)
	{
		m_conv_filename 	= restart_base_filename + "_conv.dat";
		m_u0_A_filename     = restart_base_filename + "_u0_A.dat";
		m_u0_B_filename     = restart_base_filename + "_u0_B.dat";
		m_p_i_filename    	= restart_base_filename + "_p_i.dat";
		m_r_i_filename 	    = restart_base_filename + "_r_i.dat";
		m_lambda_i_filename = restart_base_filename + "_lambda_i.dat";
		m_rho_filename 	    = restart_base_filename + "_rho.dat";
	}
}

void carl::PETSC_CG_coupled_solver::set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
					libMesh::PetscMatrix<libMesh::Number>& M_B,
					libMesh::PetscMatrix<libMesh::Number>& C_RA,
					libMesh::PetscMatrix<libMesh::Number>& C_RB,
					libMesh::PetscMatrix<libMesh::Number>& C_RR)
{
	libMesh::PerfLog perf_log("Matrix setup",MASTER_bPerfLog_CG_solver_matrix_assemble);

	coupled_solver::set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

	m_bMatricesSetUp = true;

	std::cout << "| -> Using CG " << std::endl;

//	// Calculate the preconditioner - if needed
//	if(m_bUsePreconditioner)
//	{
//		this->build_preconditioner();
//	}
//	else
//	{
//		// Use the identity matrix
//		int M = C_RR.m();
//		int N = C_RR.n();
//
//		PetscInt local_M, local_N;
//
//		MatGetLocalSize(C_RR.mat(),&local_M,&local_N);
//
//		m_PC = std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> >(new libMesh::PetscMatrix<libMesh::Number>(*m_comm));
//		m_PC->init(M,N,local_M,local_M,1,0);
//		m_PC->close();
//		MatShift(m_PC->mat(),1.0);
//	}
};

void carl::PETSC_CG_coupled_solver::set_convergence_limits(double eps_abs, double eps_rel, int convIter, double div_tol)
{
	m_CG_conv_eps_abs = eps_abs;
	m_CG_conv_eps_rel = eps_rel;
	m_CG_conv_max_n   = convIter;
	m_CG_div_tol      = div_tol;
};

void carl::PETSC_CG_coupled_solver::solve()
{
	std::cout << "| CG solver: " << std::endl;
	std::cout << "|     Initialization ..." << std::endl; std::cout.flush();
	std::cout << "|        eps abs. = " << m_CG_conv_eps_abs <<
			          ", max. iter. = " << m_CG_conv_eps_rel <<
					  ", max. iter. = " << m_CG_conv_max_n <<
			          ", div. tol. = " << m_CG_div_tol << std::endl;

	// -> Test if the parameters are set up
	homemade_assert_msg( m_bParamsSetUp , "   solve : parameters not set up!");
	homemade_assert_msg( m_bMatricesSetUp , "   solve : Matrices not set up!");
	homemade_assert_msg( m_bForcesSetUp , "   solve : Forces not set up!");

	libMesh::PerfLog perf_log("Solve",MASTER_bPerfLog_CG_solver_solve);

	int rank = m_comm->rank();

	// -> Matrix dimensions
	// Check them beforehand
	this->check_dimensions();

	unsigned int 	dim_A, dim_A_local,
					dim_B, dim_B_local,
					dim_R, dim_R_local;

	int silly_local = 0;

	dim_A = m_M_A->m();
	dim_B = m_M_B->m();
	dim_R = m_C_RA->m();

	MatGetLocalSize(m_M_A->mat(),&silly_local,NULL);
	dim_A_local = silly_local;
	MatGetLocalSize(m_M_B->mat(),&silly_local,NULL);
	dim_B_local = silly_local;
	MatGetLocalSize(m_C_RA->mat(),&silly_local,NULL);
	dim_R_local = silly_local;

	// -> Declarations - convergence parameters
	perf_log.push("Vector declarations","Initialization");


	// -> Declarations - CG vectors

	/*
	 * 		dim_A vectors :		m_sol_A, u0_A, aux_A, rhs_A
	 * 		dim_B vectors :		m_sol_B, u0_B, aux_B, rhs_B
	 */

	// Final solution
	m_sol_A->init(dim_A,dim_A_local);
	m_sol_B->init(dim_B,dim_B_local);

	// Initial solution
	libMesh::PetscVector<libMesh::Number> vec_u0_A(*m_comm,dim_A,dim_A_local);
	libMesh::PetscVector<libMesh::Number> vec_u0_B(*m_comm,dim_B,dim_B_local);

	// Auxiliary solver rhs
	libMesh::PetscVector<libMesh::Number> vec_aux_A(*m_comm,dim_A,dim_A_local);
	libMesh::PetscVector<libMesh::Number> vec_aux_B(*m_comm,dim_B,dim_B_local);


	// Set initial values
	m_sol_A->zero();
	m_sol_B->zero();

	vec_u0_A.zero();
	vec_u0_B.zero();

	vec_aux_A =  vec_u0_A;
	vec_aux_B =  vec_u0_B;

	// Set names
	PetscObjectSetName((PetscObject)m_sol_A->vec(),"sol_A");
	PetscObjectSetName((PetscObject)m_sol_B->vec(),"sol_B");

	/*
	 * 		dim_R vectors :		lambda, 	pi, 	    	ri, 	zi
	 * 							lambda_old, pi_old, qi_old, ri_old
	 * 										pi_temp
	 *
	 */

	// Coupling space vectors
	libMesh::PetscVector<libMesh::Number> vec_lambda(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> vec_lambda_zero(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> vec_aux_lambda(*m_comm,dim_R, dim_R_local);

	// Coupled system rhs
	libMesh::PetscVector<libMesh::Number> vec_coupled_rhs(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> vec_coupled_rhs_aux(*m_comm, dim_R, dim_R_local);

	// Set initial values
	vec_lambda.zero();
	vec_lambda_zero.zero();
	vec_coupled_rhs.zero();
	vec_coupled_rhs_aux.zero();
	vec_aux_lambda.zero();

	// Calculate the coupled rhs

	// u_0,I = M_I^-1 * F_I (KSP SOLVER!)
	perf_log.push("KSP solver - A","Initialization");
	m_sys_A_solver->solve(*m_F_A,vec_u0_A);
	perf_log.pop("KSP solver - A","Initialization");

	perf_log.push("KSP solver - B","Initialization");
	m_sys_B_solver->solve (*m_F_B,vec_u0_B);
	perf_log.pop("KSP solver - B","Initialization");

//	m_M_A->print_matlab("M_A.m");
//	m_C_RA->print_matlab("C_A.m");
//
//	m_M_B->print_matlab("M_B.m");
//	m_C_RB->print_matlab("C_B.m");
//
//	m_F_A->print_matlab("F_A.m");
//	m_F_B->print_matlab("F_B.m");

//	vec_u0_A.print_matlab("sol_A0.m");
//	vec_u0_B.print_matlab("sol_B0.m");

	m_C_RA->vector_mult(vec_coupled_rhs,vec_u0_A);
	m_C_RB->vector_mult(vec_coupled_rhs_aux,vec_u0_B);
	vec_coupled_rhs.add(-1,vec_coupled_rhs_aux);

	m_coupling_solver.set_system_rhs(vec_coupled_rhs);

	perf_log.pop("Vector declarations","Initialization");

	// -> Initialize the vectors
	if(!m_bUseRestart)
	{
		perf_log.push("CG vector setup","Initialization");
		if(m_bPrintRestart)
		{
			write_PETSC_vector(vec_u0_B  ,m_u0_B_filename);
			write_PETSC_vector(vec_u0_A  ,m_u0_A_filename);
		}

		if(m_bUseNullSpaceB)
		{
			// lambda_0 = F_null * F_B
			m_coupling_solver.apply_CG_nullspace_force_projection(*m_F_B,vec_lambda_zero);
//			MatMult(m_null_F,m_F_B->vec(),vec_lambda_zero.vec());
			//
			//		libMesh::PetscMatrix<libMesh::Number> dummy_Mat(m_null_F,*m_comm);
			//		dummy_Mat.print_matlab("F_null_proj.m");
			//		libMesh::PetscMatrix<libMesh::Number> dummy_Mat_bis(m_null_R,*m_comm);
			//		dummy_Mat_bis.print_matlab("R_null.m");
		}

		perf_log.pop("CG vector setup","Initialization");
	}
	else
	{
		perf_log.push("CG vector setup","Initialization");

		read_PETSC_vector(vec_u0_A  ,m_u0_A_filename);
		read_PETSC_vector(vec_u0_B  ,m_u0_B_filename);

		read_PETSC_vector(vec_lambda_zero,m_lambda_i_filename);

		perf_log.pop("CG vector setup","Initialization");
	}


	// Set initial solution
	m_coupling_solver.set_initial_sol(vec_lambda_zero);

	std::cout << "Vectors setup: end" << std::endl;
	// Solve the coupled system
	m_coupling_solver.solve();

	vec_lambda = m_coupling_solver.get_solution();

	// Create corrected solution
	// U_1 = U_0,1 - A_1^-1 * C_1 * lambda
	perf_log.push("KSP solver - A","Solution");
	m_sys_A_solver->apply_MinvZt(vec_lambda,vec_aux_A);
	*m_sol_A.get() = vec_u0_A;
	m_sol_A->add(-1,vec_aux_A);
	perf_log.pop("KSP solver - A","Solution");

	// U_2 = U_0,2 - A_2^-1 * C_2 * lambda
	perf_log.push("KSP solver - B","Solution");
	m_sys_B_solver->apply_MinvZt(vec_lambda,vec_aux_B);
	*m_sol_B.get() = vec_u0_B;
	m_sol_B->add(vec_aux_B);

	if(m_bUseNullSpaceB)
	{
		// Must add correction
		// U_2_corr = U_2 + null_space_corr * coupling_residual
		m_coupling_solver.get_residual_vector(vec_aux_lambda);
		this->add_nullspace_correction(vec_aux_lambda,*m_sol_B);
	}

	perf_log.pop("KSP solver - B","Solution");

	m_bSolved = true;
}

void carl::PETSC_CG_coupled_solver::print_convergence(std::ostream& convergenceOut)
{
	if(m_bSolved)
	{
		for(int iii = 0; iii < m_CG_conv_n; ++iii)
		{
			convergenceOut << iii << " " << m_CG_Index[iii] << std::endl;
		}
	}
}

//void carl::PETSC_CG_coupled_solver::set_null_space_projector()
//{
//	m_coupling_solver.set_solver_CG_projector(m_null_PI);
//}

void carl::PETSC_CG_coupled_solver::set_preconditioner_matrix(libMesh::PetscMatrix<libMesh::Number>& M_precond)
{
	m_coupling_solver.set_precond_matrix(M_precond);
}

void carl::PETSC_CG_coupled_solver::build_null_space_projection_matrices(libMesh::PetscMatrix<libMesh::Number>& M_sys,
		libMesh::PetscMatrix<libMesh::Number>& C_sys)
{
	homemade_assert_msg(m_bMatricesSetUp,"Must have the matrices ready!");

	m_coupling_solver.build_solver_CG_null_space_projection_matrices(M_sys,C_sys);

	m_bCreatedRigidBodyProjectors = true;
	m_bUseNullSpaceB = true;
};

void carl::PETSC_CG_coupled_solver::add_nullspace_correction(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out)
{
	m_coupling_solver.add_CG_nullspace_correction(vec_in,vec_out);
}
