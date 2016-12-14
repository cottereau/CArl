#include "CG_solver.h"
#include <petscksp.h>
#include <petsc/private/kspimpl.h>

void carl::PETSC_CG_solver::use_preconditioner(bool flag)
{
	m_bUsePreconditioner = false;
};

void carl::PETSC_CG_solver::set_info(	bool bSavePartitionInfo,
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

void carl::PETSC_CG_solver::set_restart( 	bool bUseRestart,
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

//void carl::PETSC_CG_solver::set_coordinates(	libMesh::PetscVector<libMesh::Number>& coord_vect_A,
//					libMesh::PetscVector<libMesh::Number>& coord_vect_B)
//{
//	m_coord_vect_A = &coord_vect_A;
//	m_coord_vect_B = &coord_vect_B;
//
//	m_bCoordsSetup = true;
//}

void carl::PETSC_CG_solver::set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
					libMesh::PetscMatrix<libMesh::Number>& M_B,
					libMesh::PetscMatrix<libMesh::Number>& C_RA,
					libMesh::PetscMatrix<libMesh::Number>& C_RB,
					libMesh::PetscMatrix<libMesh::Number>& C_RR)
{
	libMesh::PerfLog perf_log("Matrix setup",MASTER_bPerfLog_CG_solver_matrix_assemble);

	coupled_solver::set_matrices(M_A,M_B,C_RA,C_RB,C_RR);

	std::cout << "| -> Using CG " << std::endl;

	// Calculate the preconditioner - if needed
	if(m_bUsePreconditioner)
	{
		this->build_preconditioner();
	}
	else
	{
		// Use the identity matrix
		int M = C_RR.m();
		int N = C_RR.n();

		PetscInt local_M, local_N;

		MatGetLocalSize(C_RR.mat(),&local_M,&local_N);

		m_PC = std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> >(new libMesh::PetscMatrix<libMesh::Number>(*m_comm));
		m_PC->init(M,N,local_M,local_M,1,0);
		m_PC->close();
		MatShift(m_PC->mat(),1.0);
	}

	m_bMatricesSetUp = true;
};

void carl::PETSC_CG_solver::set_convergence_limits(double eps_abs, double eps_rel, int convIter, double div_tol)
{
	m_CG_conv_eps_abs = eps_abs;
	m_CG_conv_eps_rel = eps_rel;
	m_CG_conv_max_n   = convIter;
	m_CG_div_tol      = div_tol;
};

void carl::PETSC_CG_solver::solve()
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

	double coupled_residual_norm = 0;
	int    iter_nb = 0;

	if(m_bUseRestart)
	{
		if(rank == 0)
		{
			std::ifstream param_data(m_conv_filename);
			param_data >> iter_nb >> coupled_residual_norm;
			param_data.close();
		}

		m_comm->broadcast(iter_nb);
		m_comm->broadcast(coupled_residual_norm);
	}

	// -> Declarations - CG vectors

	m_sol_A->init(dim_A,dim_A_local);
	m_sol_B->init(dim_B,dim_B_local);

	PetscObjectSetName((PetscObject)m_sol_A->vec(),"sol_A");
	PetscObjectSetName((PetscObject)m_sol_B->vec(),"sol_B");
	m_sol_A->zero();
	m_sol_B->zero();

	// Initial solution vectors
	libMesh::PetscVector<libMesh::Number> u0_A(*m_comm, dim_A, dim_A_local);
	libMesh::PetscVector<libMesh::Number> u0_B(*m_comm, dim_B, dim_B_local);
	u0_A.zero();
	u0_B.zero();

	// Auxiliary solver vectors
	libMesh::PetscVector<libMesh::Number> w_A(*m_comm, dim_A, dim_A_local);
	libMesh::PetscVector<libMesh::Number> w_B(*m_comm, dim_B, dim_B_local);
	w_A.zero();
	w_B.zero();

	// Auxiliary solver rhs
	libMesh::PetscVector<libMesh::Number> rhs_A(*m_comm, dim_A, dim_A_local);
	libMesh::PetscVector<libMesh::Number> rhs_B(*m_comm, dim_B, dim_B_local);
	rhs_A.zero();
	rhs_B.zero();

	// Coupling space vectors
	libMesh::PetscVector<libMesh::Number> lambda_vec(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> lambda_vec_old(*m_comm, dim_R, dim_R_local);
	lambda_vec.zero();
	lambda_vec_old.zero();

	// Auxiliary CG vectors
	libMesh::PetscVector<libMesh::Number> r_i(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> r_i_old(*m_comm, dim_R, dim_R_local);
	r_i.zero();
	r_i_old.zero();
	PetscObjectSetName((PetscObject)r_i_old.vec(),"r_i");

	libMesh::PetscVector<libMesh::Number> z_i(*m_comm, dim_R, dim_R_local);
	z_i.zero();
	PetscObjectSetName((PetscObject)z_i.vec(),"z_i");

	libMesh::PetscVector<libMesh::Number> p_i(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> p_i_old(*m_comm, dim_R, dim_R_local);
	p_i.zero();
	p_i_old.zero();
	PetscObjectSetName((PetscObject)p_i_old.vec(),"p_i");

	// Auxiliary solver rhs
	libMesh::PetscVector<libMesh::Number> q_i(*m_comm, dim_R, dim_R_local);
	q_i.zero();

	// Search parameters
	double alpha_i = 0;
	double beta_ip = 0;
	double rho_i = 0;
	double rho_i_old = 0;
	double rho_0 = 0;
	double dummy_aux = 0;

	perf_log.pop("Vector declarations","Initialization");

	// Set up the null spaces


	// Solvers
	perf_log.push("KSP solvers setup","Initialization");
//	KSP_linear_solver KSP_Solver_M_A(*m_comm);
//	KSP_linear_solver KSP_Solver_M_B(*m_comm);

//	KSP_Solver_M_A.set_solver(*m_M_A, m_ksp_name_A.c_str());
//	KSP_Solver_M_B.set_solver(*m_M_B, m_ksp_name_B.c_str());

//	KSP_Solver_M_A.print_type();
//	KSP_Solver_M_B.print_type();

	KSP PETSc_ksp_A, PETSc_ksp_B;

	KSPCreate(PETSC_COMM_WORLD, &PETSc_ksp_A);
	KSPCreate(PETSC_COMM_WORLD, &PETSc_ksp_B);
	KSPSetOperators(PETSc_ksp_A,m_M_A->mat(),m_M_A->mat());
	KSPSetOperators(PETSc_ksp_B,m_M_B->mat(),m_M_B->mat());

	KSPSetFromOptions(PETSc_ksp_A);
	KSPSetFromOptions(PETSc_ksp_B);


	perf_log.pop("KSP solvers setup","Initialization");

	// -> Initialize the vectors
	if(!m_bUseRestart)
	{
		// u_0,I = M_I^-1 * F_I (KSP SOLVER!)
		perf_log.push("KSP solver - A","Initialization");
//		KSP_Solver_M_A.solve (u0_A, *m_F_A);
		KSPSolve(PETSc_ksp_A, m_F_A->vec(), u0_A.vec());
		write_PETSC_vector(u0_A  ,m_u0_A_filename);
		perf_log.pop("KSP solver - A","Initialization");

		perf_log.push("KSP solver - B","Initialization");
//		KSP_Solver_M_B.solve (u0_B, *m_F_B);
		KSPSolve(PETSc_ksp_B, m_F_B->vec(), u0_B.vec());
		write_PETSC_vector(u0_B  ,m_u0_B_filename);
		perf_log.pop("KSP solver - B","Initialization");

		perf_log.push("CG vector setup","Initialization");

		// lambda_0 = 0
		lambda_vec_old.zero();

		// r_0 = C_A * u_0,A + C_B * u_0,B
		m_C_RA->vector_mult(r_i_old,u0_A);
		m_C_RB->vector_mult_add(r_i_old,u0_B);
		std::cout << r_i_old.l2_norm() << std::endl;

		// z_0 = [ M^-1 ] * r_0 = PC * r_0
		m_PC->vector_mult(z_i,r_i_old);
		std::cout << m_PC->linfty_norm() << " " << m_PC->l1_norm() << std::endl;

		// p_0 = z_0
		p_i_old = z_i;

		// rho_0 = r_0 * z_0
		rho_i_old = r_i_old.dot(z_i);
		rho_0 = rho_i_old;
		perf_log.pop("CG vector setup","Initialization");
	}
	else
	{
		perf_log.push("CG vector setup","Initialization");

		read_PETSC_vector(u0_A  ,m_u0_A_filename);
		read_PETSC_vector(u0_B  ,m_u0_B_filename);

		read_PETSC_vector(lambda_vec_old  ,m_lambda_i_filename);

		read_PETSC_vector(r_i_old  ,m_r_i_filename);
		read_PETSC_vector(p_i_old  ,m_p_i_filename);

		if(rank == 0)
		{
			std::ifstream param_data(m_rho_filename);
			param_data >> rho_0 >> rho_i_old;
			param_data.close();
		}

		m_comm->broadcast(rho_0);
		m_comm->broadcast(rho_i_old);
		perf_log.pop("CG vector setup","Initialization");
	}

	libMesh::PerfData timing_data;
	bool bKeepRunning = true;
	bool bConverged = false;

	KSPConvergedReason conv_reason;

	while (bKeepRunning)
	{
		std::cout << "|     Iter no. " << iter_nb << " " << std::endl; std::cout.flush();

		// -> Generate new correction
		// w_i,I = M_I^-1 * C_I^T* p_i
		MatMultTranspose(m_C_RA->mat(),p_i_old.vec(),rhs_A.vec());
		perf_log.push("KSP solver - A","Coupled CG iterations");
//		KSP_Solver_M_A.solve(w_A,rhs_A);
		KSPSolve(PETSc_ksp_A, rhs_A.vec(), w_A.vec());
		perf_log.pop("KSP solver - A","Coupled CG iterations");
		timing_data = perf_log.get_perf_data("KSP solver - A","Coupled CG iterations");
		KSPGetConvergedReason(PETSc_ksp_A,&conv_reason);
		std::cout << "|        Solver A time : " << timing_data.tot_time/(iter_nb + 1) << std::endl;
		std::cout << "|          Converged ? : " << conv_reason << std::endl;

		MatMultTranspose(m_C_RB->mat(),p_i_old.vec(),rhs_B.vec());
		perf_log.push("KSP solver - B","Coupled CG iterations");
//		KSP_Solver_M_B.solve(w_B,rhs_B);
		KSPSolve(PETSc_ksp_B, rhs_B.vec(), w_B.vec());
		perf_log.pop("KSP solver - B","Coupled CG iterations");
		timing_data = perf_log.get_perf_data("KSP solver - B","Coupled CG iterations");
		KSPGetConvergedReason(PETSc_ksp_B,&conv_reason);
		std::cout << "|        Solver B time : " << timing_data.tot_time/(iter_nb + 1) << std::endl;
		std::cout << "|          Converged ? : " << conv_reason << std::endl;

		perf_log.push("New correction","Coupled CG iterations");

		// q_i = C_A * w_i,A + C_B * w_i,B
		m_C_RA->vector_mult(q_i,w_A);
		m_C_RB->vector_mult_add(q_i,w_B);

		// alpha_i = rho_i / ( p_i * q_i )
		dummy_aux = p_i_old.dot(q_i);

		alpha_i = rho_i_old / dummy_aux;

		std::cout << "|        " << rho_i_old << " " <<  dummy_aux << " " << alpha_i << std::endl;

		lambda_vec = lambda_vec_old;
		lambda_vec.add(alpha_i,p_i_old);

		perf_log.pop("New correction","Coupled CG iterations");

		// -> Update auxiliary vectors
		perf_log.push("Update auxiliary vectors","Coupled CG iterations");
		// r_(i+1) = r_i - alpha_i * p_i
		r_i = r_i_old;
		r_i.add(-alpha_i,p_i_old);

		// z_(i+1) = PC * r_(i+1)
		m_PC->vector_mult(z_i,r_i);

		// rho_(i+1) = r_(i+1) * z_(i+1)
		rho_i = r_i.dot(z_i);

		// beta_(i+1) = rho_(i+1) / rho_i
		std::cout << rho_i << " " <<  rho_i_old << std::endl;
		std::cout << r_i.l2_norm() << " " << z_i.l2_norm() << std::endl;
		beta_ip = rho_i / rho_i_old;

		std::cout << "|        " << beta_ip << std::endl;

		// p_(i+1) = z_(i+1) + beta_(i+1) * p_i
		p_i = z_i;
		p_i.add(beta_ip,p_i_old);
		perf_log.pop("Update auxiliary vectors","Coupled CG iterations");

		// -> Check the convergence
		++iter_nb;

		// Absolute convergence
		if(rho_i < m_CG_conv_eps_abs)
		{
			bKeepRunning = false;
			bConverged = true;
		}

		// Relative convergence
		if(rho_i < m_CG_conv_eps_rel * rho_0)
		{
			bKeepRunning = false;
			bConverged = true;
		}

		// Iteration divergence
		if(iter_nb > m_CG_conv_max_n)
		{
			bKeepRunning = false;
			bConverged = false;
		}

		// Residue divergence
		if(rho_i > m_CG_div_tol * rho_0)
		{
			bKeepRunning = false;
			bConverged = false;
		}

		// Update the values and vectors
		rho_i_old = rho_i;
		p_i_old = p_i;
		lambda_vec_old = lambda_vec;
		r_i_old = r_i;

		if(m_bPrintRestart)
		{
			std::cout << "|        Writing to files " << m_lambda_i_filename << ", etc ..." << std::endl;
			write_PETSC_vector(p_i_old  ,m_p_i_filename);
			write_PETSC_vector(r_i_old  ,m_r_i_filename);
			write_PETSC_vector(lambda_vec_old,m_lambda_i_filename);

			if(rank == 0)
			{
				std::ofstream param_data(m_conv_filename);
				param_data << iter_nb << " " << rho_i_old << std::endl;
				param_data.close();

				std::ofstream rho_data(m_rho_filename);
				param_data << rho_0 << " " << rho_i_old << std::endl;
				param_data.close();
			}
		}
		std::cout << "|" << std::endl;
	}
	m_CG_conv_n = iter_nb;

	std::cout << "|     nb. of iterations : " << iter_nb << std::endl;
	std::cout << "|     || r_n ||_PC : " << rho_i_old << std::endl;

	homemade_assert_msg(bConverged," -> CG solver diverged!");

	// Create corrected solution
	// U_I = U_0,I - A_I^-1 * C_I * lambda
	perf_log.push("KSP solver - A","Solution");

	MatMultTranspose(m_C_RA->mat(),lambda_vec.vec(),rhs_A.vec());
//	KSP_Solver_M_A.solve(w_A,rhs_A);
	KSPSolve(PETSc_ksp_A, rhs_A.vec(), w_A.vec());
	*m_sol_A.get() = u0_A;
	m_sol_A->add(-1,w_A);
	perf_log.pop("KSP solver - A","Solution");

	perf_log.push("KSP solver - B","Solution");
	MatMultTranspose(m_C_RB->mat(),lambda_vec.vec(),rhs_B.vec());
	KSPSolve(PETSc_ksp_B, rhs_B.vec(), w_B.vec());
	*m_sol_B.get() = u0_B;
	m_sol_B->add(1,w_B);
	perf_log.pop("KSP solver - B","Solution");

//	MatNullSpaceDestroy(&M_A_nullsp);
//	MatNullSpaceDestroy(&M_B_nullsp);

	m_bSolved = true;
}

void carl::PETSC_CG_solver::print_convergence(std::ostream& convergenceOut)
{
	if(m_bSolved)
	{
		for(int iii = 0; iii < m_CG_conv_n; ++iii)
		{
			convergenceOut << iii << " " << m_CG_Index[iii] << std::endl;
		}
	}
}
