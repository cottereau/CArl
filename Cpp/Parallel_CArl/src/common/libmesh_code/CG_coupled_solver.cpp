#include <petscksp.h>
#include <petsc/private/kspimpl.h>
#include "CG_coupled_solver.h"

void carl::PETSC_CG_coupled_solver::use_preconditioner(bool flag)
{
	m_bUsePreconditioner = flag;
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

//void carl::PETSC_CG_coupled_solver::set_coordinates(	libMesh::PetscVector<libMesh::Number>& coord_vect_A,
//					libMesh::PetscVector<libMesh::Number>& coord_vect_B)
//{
//	m_coord_vect_A = &coord_vect_A;
//	m_coord_vect_B = &coord_vect_B;
//
//	m_bCoordsSetup = true;
//}

void  carl::PETSC_CG_coupled_solver::build_preconditioner()
{
	Mat dummy_mat;
	MatRARt(m_M_A->mat(),m_C_RA->mat(),MAT_INITIAL_MATRIX,1,&dummy_mat);
	MatRARt(m_M_B->mat(),m_C_RB->mat(),MAT_INITIAL_MATRIX,1,&m_PC_PETSc);

	MatAXPY(m_PC_PETSc,1,dummy_mat,DIFFERENT_NONZERO_PATTERN);
	m_PC = std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> >(new libMesh::PetscMatrix<libMesh::Number>(m_PC_PETSc,*m_comm));

	MatDestroy(&dummy_mat);
};

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
};

void carl::PETSC_CG_coupled_solver::set_convergence_limits(double eps_abs, double eps_rel, int convIter, double div_tol)
{
	m_CG_conv_eps_abs = eps_abs;
	m_CG_conv_eps_rel = eps_rel;
	m_CG_conv_max_n   = convIter;
	m_CG_div_tol      = div_tol;
};

//void carl::PETSC_CG_coupled_solver::solve()
//{
//	std::cout << "| CG solver: " << std::endl;
//	std::cout << "|     Initialization ..." << std::endl; std::cout.flush();
//	std::cout << "|        eps abs. = " << m_CG_conv_eps_abs <<
//			          ", max. iter. = " << m_CG_conv_eps_rel <<
//					  ", max. iter. = " << m_CG_conv_max_n <<
//			          ", div. tol. = " << m_CG_div_tol << std::endl;
//
//	// -> Test if the parameters are set up
//	homemade_assert_msg( m_bParamsSetUp , "   solve : parameters not set up!");
//	homemade_assert_msg( m_bMatricesSetUp , "   solve : Matrices not set up!");
//	homemade_assert_msg( m_bForcesSetUp , "   solve : Forces not set up!");
//
//	libMesh::PerfLog perf_log("Solve",MASTER_bPerfLog_CG_solver_solve);
//
//	int rank = m_comm->rank();
//
//	// -> Matrix dimensions
//	// Check them beforehand
//	this->check_dimensions();
//
//	unsigned int 	dim_A, dim_A_local,
//					dim_B, dim_B_local,
//					dim_R, dim_R_local;
//
//	int silly_local = 0;
//
//	dim_A = m_M_A->m();
//	dim_B = m_M_B->m();
//	dim_R = m_C_RA->m();
//
//	MatGetLocalSize(m_M_A->mat(),&silly_local,NULL);
//	dim_A_local = silly_local;
//	MatGetLocalSize(m_M_B->mat(),&silly_local,NULL);
//	dim_B_local = silly_local;
//	MatGetLocalSize(m_C_RA->mat(),&silly_local,NULL);
//	dim_R_local = silly_local;
//
//	// -> Declarations - convergence parameters
//	perf_log.push("Vector declarations","Initialization");
//
//	double coupled_residual_norm = 0;
//	int    iter_nb = 0;
//
//	if(m_bUseRestart)
//	{
//		if(rank == 0)
//		{
//			std::ifstream param_data(m_conv_filename);
//			param_data >> iter_nb >> coupled_residual_norm;
//			param_data.close();
//		}
//
//		m_comm->broadcast(iter_nb);
//		m_comm->broadcast(coupled_residual_norm);
//	}
//
//	// -> Declarations - CG vectors
//
//	m_sol_A->init(dim_A,dim_A_local);
//	m_sol_B->init(dim_B,dim_B_local);
//
//	PetscObjectSetName((PetscObject)m_sol_A->vec(),"sol_A");
//	PetscObjectSetName((PetscObject)m_sol_B->vec(),"sol_B");
//	m_sol_A->zero();
//	m_sol_B->zero();
//
//	// Initial solution vectors
//	libMesh::PetscVector<libMesh::Number> u0_A(*m_comm, dim_A, dim_A_local);
//	libMesh::PetscVector<libMesh::Number> u0_B(*m_comm, dim_B, dim_B_local);
//	u0_A.zero();
//	u0_B.zero();
//
//	// Auxiliary solver vectors
//	libMesh::PetscVector<libMesh::Number> w_A(*m_comm, dim_A, dim_A_local);
//	libMesh::PetscVector<libMesh::Number> w_B(*m_comm, dim_B, dim_B_local);
//	w_A.zero();
//	w_B.zero();
//
//	// Auxiliary solver rhs
//	libMesh::PetscVector<libMesh::Number> rhs_A(*m_comm, dim_A, dim_A_local);
//	libMesh::PetscVector<libMesh::Number> rhs_B(*m_comm, dim_B, dim_B_local);
//	rhs_A.zero();
//	rhs_B.zero();
//
//	// Coupling space vectors
//	libMesh::PetscVector<libMesh::Number> lambda_vec(*m_comm, dim_R, dim_R_local);
//	libMesh::PetscVector<libMesh::Number> lambda_vec_old(*m_comm, dim_R, dim_R_local);
//	lambda_vec.zero();
//	lambda_vec_old.zero();
//
//	// Auxiliary CG vectors
//	libMesh::PetscVector<libMesh::Number> r_i(*m_comm, dim_R, dim_R_local);
//	libMesh::PetscVector<libMesh::Number> r_i_old(*m_comm, dim_R, dim_R_local);
//	r_i.zero();
//	r_i_old.zero();
//	PetscObjectSetName((PetscObject)r_i_old.vec(),"r_i");
//
//	libMesh::PetscVector<libMesh::Number> z_i(*m_comm, dim_R, dim_R_local);
//	z_i.zero();
//	PetscObjectSetName((PetscObject)z_i.vec(),"z_i");
//
//	libMesh::PetscVector<libMesh::Number> p_i_temp(*m_comm, dim_R, dim_R_local);
//	libMesh::PetscVector<libMesh::Number> p_i(*m_comm, dim_R, dim_R_local);
//	libMesh::PetscVector<libMesh::Number> p_i_old(*m_comm, dim_R, dim_R_local);
//	p_i_temp.zero();
//	p_i.zero();
//	p_i_old.zero();
//	PetscObjectSetName((PetscObject)p_i_old.vec(),"p_i");
//
//	// Auxiliary solver rhs
//	libMesh::PetscVector<libMesh::Number> q_i(*m_comm, dim_R, dim_R_local);
//	q_i.zero();
//
//	// Search parameters
//	double alpha_i = 0;
//	double beta_ip = 0;
//	double rho_i = 0;
//	double rho_i_old = 0;
//	double rho_0 = 0;
//	double dummy_aux = 0;
//
//	perf_log.pop("Vector declarations","Initialization");
//
//	// Solvers
//	perf_log.push("KSP solvers setup","Initialization");
//	KSP_linear_solver KSP_Solver_M_A(*m_comm);
//	KSP_linear_solver KSP_Solver_M_B(*m_comm);
//
//	KSP_Solver_M_A.set_solver(*m_M_A, m_ksp_name_A.c_str());
//	KSP_Solver_M_B.set_solver(*m_M_B, m_ksp_name_B.c_str());
//
//	KSP_Solver_M_A.print_type();
//	KSP_Solver_M_B.print_type();
//
//	perf_log.pop("KSP solvers setup","Initialization");
//
//	// -> Initialize the vectors
//	if(!m_bUseRestart)
//	{
//		// u_0,I = M_I^-1 * F_I (KSP SOLVER!)
//		perf_log.push("KSP solver - A","Initialization");
//		KSP_Solver_M_A.solve (u0_A, *m_F_A);
//		write_PETSC_vector(u0_A  ,m_u0_A_filename);
//		perf_log.pop("KSP solver - A","Initialization");
//
//		perf_log.push("KSP solver - B","Initialization");
//		KSP_Solver_M_B.solve (u0_B, *m_F_B);
//		write_PETSC_vector(u0_B  ,m_u0_B_filename);
//		perf_log.pop("KSP solver - B","Initialization");
//
//		perf_log.push("CG vector setup","Initialization");
//
//		// lambda_0 = F_null * F_B
//		MatMult(m_null_F,m_F_B->vec(),lambda_vec_old.vec());
//
//		// r_0 = C_A * u_0,A + C_B * u_0,B
//		m_C_RA->vector_mult(r_i_old,u0_A);
//		m_C_RB->vector_mult_add(r_i_old,u0_B);
//		std::cout << r_i_old.l2_norm() << std::endl;
//
//		// z_0 = [ M^-1 ] * r_0 = PC * r_0
//		m_PC->vector_mult(z_i,r_i_old);
//		std::cout << m_PC->linfty_norm() << " " << m_PC->l1_norm() << std::endl;
//
//		// p_0 = PI_/RB * z_0
//		MatMult(m_null_PI,z_i.vec(),p_i_old.vec());
//
//		// rho_0 = r_0 * z_0
//		rho_i_old = r_i_old.dot(z_i);
//		rho_0 = rho_i_old;
//		perf_log.pop("CG vector setup","Initialization");
//	}
//	else
//	{
//		perf_log.push("CG vector setup","Initialization");
//
//		read_PETSC_vector(u0_A  ,m_u0_A_filename);
//		read_PETSC_vector(u0_B  ,m_u0_B_filename);
//
//		read_PETSC_vector(lambda_vec_old  ,m_lambda_i_filename);
//
//		read_PETSC_vector(r_i_old  ,m_r_i_filename);
//		read_PETSC_vector(p_i_old  ,m_p_i_filename);
//
//		if(rank == 0)
//		{
//			std::ifstream param_data(m_rho_filename);
//			param_data >> rho_0 >> rho_i_old;
//			param_data.close();
//		}
//
//		m_comm->broadcast(rho_0);
//		m_comm->broadcast(rho_i_old);
//		perf_log.pop("CG vector setup","Initialization");
//	}
//
//	libMesh::PerfData timing_data;
//	bool bKeepRunning = true;
//	bool bConverged = false;
//
//	KSPConvergedReason conv_reason;
//
//	while (bKeepRunning)
//	{
//		std::cout << "|     Iter no. " << iter_nb << " " << std::endl; std::cout.flush();
//
//		// -> Generate new correction
//		// w_i,I = M_I^-1 * C_I^T* p_i
//		MatMultTranspose(m_C_RA->mat(),p_i_old.vec(),rhs_A.vec());
//		perf_log.push("KSP solver - A","Coupled CG iterations");
//		KSP_Solver_M_A.solve(w_A,rhs_A);
//		perf_log.pop("KSP solver - A","Coupled CG iterations");
//		timing_data = perf_log.get_perf_data("KSP solver - A","Coupled CG iterations");
//		KSPGetConvergedReason(KSP_Solver_M_A.solver_ptr()->ksp(),&conv_reason);
//		std::cout << "|        Solver A time : " << timing_data.tot_time/(iter_nb + 1) << std::endl;
//		std::cout << "|          Converged ? : " << conv_reason << std::endl;
//
//		MatMultTranspose(m_C_RB->mat(),p_i_old.vec(),rhs_B.vec());
//		perf_log.push("KSP solver - B","Coupled CG iterations");
//		KSP_Solver_M_B.solve(w_B,rhs_B);
//		perf_log.pop("KSP solver - B","Coupled CG iterations");
//		timing_data = perf_log.get_perf_data("KSP solver - B","Coupled CG iterations");
//		KSPGetConvergedReason(KSP_Solver_M_B.solver_ptr()->ksp(),&conv_reason);
//		std::cout << "|        Solver B time : " << timing_data.tot_time/(iter_nb + 1) << std::endl;
//		std::cout << "|          Converged ? : " << conv_reason << std::endl;
//
//		perf_log.push("New correction","Coupled CG iterations");
//
//		// q_i = C_A * w_i,A + C_B * w_i,B
//		m_C_RA->vector_mult(q_i,w_A);
//		m_C_RB->vector_mult_add(q_i,w_B);
//
//		// alpha_i = rho_i / ( p_i * q_i )
//		dummy_aux = p_i_old.dot(q_i);
//
//		alpha_i = rho_i_old / dummy_aux;
//
//		std::cout << "| !!!    " << rho_i_old << " " <<  dummy_aux << " " << alpha_i << std::endl;
//
//		lambda_vec = lambda_vec_old;
//		lambda_vec.add(alpha_i,p_i_old);
//
//		perf_log.pop("New correction","Coupled CG iterations");
//
//		// -> Update auxiliary vectors
//		perf_log.push("Update auxiliary vectors","Coupled CG iterations");
//		// r_(i+1) = r_i - alpha_i * p_i
//		r_i = r_i_old;
//		r_i.add(-alpha_i,p_i_old);
//
//		// z_(i+1) = PC * r_(i+1)
//		m_PC->vector_mult(z_i,r_i);
//
//		// rho_(i+1) = r_(i+1) * z_(i+1)
//		rho_i = r_i.dot(z_i);
//
//		// beta_(i+1) = rho_(i+1) / rho_i
//		std::cout << rho_i << " " <<  rho_i_old << std::endl;
//		std::cout << r_i.l2_norm() << " " << z_i.l2_norm() << std::endl;
//		beta_ip = rho_i / rho_i_old;
//
//		std::cout << "|        " << beta_ip << std::endl;
//
//		// p_temp_(i+1) = z_(i+1) + beta_(i+1) * p_i
//		p_i_temp = z_i;
//		p_i_temp.add(beta_ip,p_i_old);
//
//		// p_(i+1) = PI_/RB * p_temp_(i+1)
//		MatMult(m_null_PI,p_i_temp.vec(),p_i.vec());
//
//		perf_log.pop("Update auxiliary vectors","Coupled CG iterations");
//
//		// -> Check the convergence
//		++iter_nb;
//
//		// Absolute convergence
//		if(rho_i < m_CG_conv_eps_abs)
//		{
//			bKeepRunning = false;
//			bConverged = true;
//		}
//
//		// Relative convergence
//		if(rho_i < m_CG_conv_eps_rel * rho_0)
//		{
//			bKeepRunning = false;
//			bConverged = true;
//		}
//
//		// Iteration divergence
//		if(iter_nb > m_CG_conv_max_n)
//		{
//			bKeepRunning = false;
//			bConverged = false;
//		}
//
//		// Residue divergence
//		if(rho_i > m_CG_div_tol * rho_0)
//		{
//			bKeepRunning = false;
//			bConverged = false;
//		}
//
//		// Update the values and vectors
//		rho_i_old = rho_i;
//		p_i_old = p_i;
//		lambda_vec_old = lambda_vec;
//		r_i_old = r_i;
//
//		if(m_bPrintRestart)
//		{
//			std::cout << "|        Writing to files " << m_lambda_i_filename << ", etc ..." << std::endl;
//			write_PETSC_vector(p_i_old  ,m_p_i_filename);
//			write_PETSC_vector(r_i_old  ,m_r_i_filename);
//			write_PETSC_vector(lambda_vec_old,m_lambda_i_filename);
//
//			if(rank == 0)
//			{
//				std::ofstream param_data(m_conv_filename);
//				param_data << iter_nb << " " << rho_i_old << std::endl;
//				param_data.close();
//
//				std::ofstream rho_data(m_rho_filename);
//				param_data << rho_0 << " " << rho_i_old << std::endl;
//				param_data.close();
//			}
//		}
//		std::cout << "|" << std::endl;
//	}
//	m_CG_conv_n = iter_nb;
//
//	std::cout << "|     nb. of iterations : " << iter_nb << std::endl;
//	std::cout << "|     || r_n ||_PC : " << rho_i_old << std::endl;
//
//	homemade_assert_msg(bConverged," -> CG solver diverged!");
//
//	// Create corrected solution
//	// U_I = U_0,I - A_I^-1 * C_I * lambda
//	perf_log.push("KSP solver - A","Solution");
//
//	MatMultTranspose(m_C_RA->mat(),lambda_vec.vec(),rhs_A.vec());
//	KSP_Solver_M_A.solve(w_A,rhs_A);
//	*m_sol_A.get() = u0_A;
//	m_sol_A->add(-1,w_A);
//	perf_log.pop("KSP solver - A","Solution");
//
//	perf_log.push("KSP solver - B","Solution");
//	MatMultTranspose(m_C_RB->mat(),lambda_vec.vec(),rhs_B.vec());
//	KSP_Solver_M_B.solve(w_B,rhs_B);
//	*m_sol_B.get() = u0_B;
//	m_sol_B->add(1,w_B);
//	perf_log.pop("KSP solver - B","Solution");
//
//	m_bSolved = true;
//}

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
	libMesh::PetscVector<libMesh::Number> vec_rhs_A(*m_comm,dim_A,dim_A_local);
	libMesh::PetscVector<libMesh::Number> vec_rhs_B(*m_comm,dim_B,dim_B_local);

	// Set initial values
	m_sol_A->zero();
	m_sol_B->zero();

	vec_u0_A.zero();
	vec_u0_B.zero();

	vec_rhs_A =  vec_u0_A;
	vec_rhs_B =  vec_u0_B;

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
	libMesh::PetscVector<libMesh::Number> vec_lambda_old(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> vec_lambda_aux(*m_comm, dim_R, dim_R_local);

	// Auxiliary solver vectors
	libMesh::PetscVector<libMesh::Number> vec_aux_A(*m_comm,dim_R,dim_R_local);
	libMesh::PetscVector<libMesh::Number> vec_aux_B(*m_comm,dim_R,dim_R_local);

	// Auxiliary CG vectors
	libMesh::PetscVector<libMesh::Number> vec_pi(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> vec_pi_old(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> vec_pi_temp(*m_comm, dim_R, dim_R_local);

	libMesh::PetscVector<libMesh::Number> vec_qi_old(*m_comm, dim_R, dim_R_local);

	libMesh::PetscVector<libMesh::Number> vec_ri(*m_comm, dim_R, dim_R_local);
	libMesh::PetscVector<libMesh::Number> vec_ri_old(*m_comm, dim_R, dim_R_local);

	libMesh::PetscVector<libMesh::Number> vec_zi(*m_comm, dim_R, dim_R_local);

	// Set initial values
	vec_lambda.zero();
	vec_lambda_old 	= vec_lambda;
	vec_lambda_aux 	= vec_lambda;
	vec_aux_A 		= vec_lambda;
	vec_aux_B 		= vec_lambda;
	vec_pi 			= vec_lambda; vec_pi_old = vec_lambda; vec_pi_temp = vec_lambda;
	vec_qi_old 		= vec_lambda;
	vec_ri 			= vec_lambda; vec_ri_old = vec_lambda;
	vec_zi 			= vec_lambda;

	// Set names
	PetscObjectSetName((PetscObject)vec_ri_old.vec(),"r_i");
	PetscObjectSetName((PetscObject)vec_pi_old.vec(),"p_i");

	// Search parameters
	double alpha_i_old = 0;
	double beta_i = 0;
	double rho_i = 0;
	double rho_i_old = 0;
	double rho_0 = 0;
	double aux_double = 0;

	perf_log.pop("Vector declarations","Initialization");

	// Solvers
	perf_log.push("KSP solvers setup","Initialization");
	KSP_linear_solver KSP_Solver_M_A(*m_comm);
	KSP_linear_solver KSP_Solver_M_B(*m_comm);

	KSP_Solver_M_A.set_solver(*m_M_A, m_ksp_name_A.c_str());
	KSP_Solver_M_B.set_solver(*m_M_B, m_ksp_name_B.c_str());

	KSP_Solver_M_A.set_coupling_matrix(*m_C_RA);
	KSP_Solver_M_B.set_coupling_matrix(*m_C_RB);

	KSP_Solver_M_A.print_type();
	KSP_Solver_M_B.print_type();

	perf_log.pop("KSP solvers setup","Initialization");

	// -> Initialize the vectors
	if(!m_bUseRestart)
	{
		// u_0,I = M_I^-1 * F_I (KSP SOLVER!)
		perf_log.push("KSP solver - A","Initialization");
		KSP_Solver_M_A.solve (vec_u0_A, *m_F_A);
		write_PETSC_vector(vec_u0_A  ,m_u0_A_filename);
		perf_log.pop("KSP solver - A","Initialization");

		perf_log.push("KSP solver - B","Initialization");
		KSP_Solver_M_B.solve (vec_u0_B, *m_F_B);
		write_PETSC_vector(vec_u0_B  ,m_u0_B_filename);
		perf_log.pop("KSP solver - B","Initialization");

		perf_log.push("CG vector setup","Initialization");

		// lambda_0 = F_null * F_B
		MatMult(m_null_F,m_F_B->vec(),vec_lambda_old.vec());

		// r_0 = C_A * u_0,A - C_B * u_0,B - A^+ lambda_0
		m_C_RA->vector_mult(vec_ri_old,vec_u0_A);
		m_C_RB->vector_mult_add(vec_lambda_aux,vec_u0_B);
		vec_ri_old.add(-1,vec_lambda_aux);

		KSP_Solver_M_A.apply_ZMiZt(vec_lambda_old, vec_lambda_aux);
		vec_ri_old.add(-1,vec_lambda_aux);
		KSP_Solver_M_B.apply_ZMiZt(vec_lambda_old, vec_lambda_aux);
		vec_ri_old.add(-1,vec_lambda_aux);

		// z_0 = M_PC * r_0
		m_PC->vector_mult(vec_zi,vec_ri_old);

		// p_0 = PI_null * z_0
		MatMult(m_null_PI,vec_zi.vec(),vec_pi_old.vec());

		// rho_0 = r_0 * z_0
		rho_i_old = vec_ri_old.dot(vec_zi);
		rho_0 = rho_i_old;

		std::cout << "|     MAtrix norms " << std::endl;
		std::cout << "|        K1 (L1, Linfty) :" << m_M_A->l1_norm() << " " << m_M_A->linfty_norm() << std::endl;
		std::cout << "|        K2 (L1, Linfty) :" << m_M_B->l1_norm() << " " << m_M_B->linfty_norm() << std::endl;
		std::cout << "|        C1 (L1, Linfty) :" << m_C_RA->l1_norm() << " " << m_C_RA->linfty_norm() << std::endl;
		std::cout << "|        C2 (L1, Linfty) :" << m_C_RB->l1_norm() << " " << m_C_RB->linfty_norm() << std::endl;
		std::cout << "|" << std::endl;

//		print_matrix_matlab(*m_M_A,"K1.m");
//		print_matrix_matlab(*m_M_B,"K2.m");
//		print_matrix_matlab(*m_C_RA,"C1.m");
//		print_matrix_matlab(*m_C_RB,"C2.m");

		std::cout << "|     Initial values " << std::endl;
		std::cout << "|        Lambda (L2) :" << vec_lambda_old.l2_norm() << std::endl;
		std::cout << "|        r (L2)      :" << vec_ri_old.l2_norm() << std::endl;
		std::cout << "|        z (L2)      :" << vec_zi.l2_norm() << std::endl;
		std::cout << "|        p (L2)      :" << vec_pi_old.l2_norm() << std::endl;
		std::cout << "|        rho         :" << rho_i_old << std::endl;
		std::cout << "|" << std::endl;

		perf_log.pop("CG vector setup","Initialization");
	}
	else
	{
		perf_log.push("CG vector setup","Initialization");

		read_PETSC_vector(vec_u0_A  ,m_u0_A_filename);
		read_PETSC_vector(vec_u0_B  ,m_u0_B_filename);

		read_PETSC_vector(vec_lambda_old  ,m_lambda_i_filename);

		read_PETSC_vector(vec_ri_old  ,m_r_i_filename);
		read_PETSC_vector(vec_pi_old  ,m_p_i_filename);

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

	while (bKeepRunning)
	{
		std::cout << "|     Iter no. " << iter_nb << " " << std::endl; std::cout.flush();

		// -> Generate new correction
		// w_i,I = C_I * M_I^-1 * C_I^T* p_i
		perf_log.push("KSP solver - A","Coupled CG iterations");
		KSP_Solver_M_A.apply_ZMiZt(vec_pi_old,vec_aux_A);
		perf_log.pop("KSP solver - A","Coupled CG iterations");

		timing_data = perf_log.get_perf_data("KSP solver - A","Coupled CG iterations");
		std::cout << "|        Solver A time : " << timing_data.tot_time/(iter_nb + 1) << std::endl;
		std::cout << "|        Converged ?   : " << KSP_Solver_M_A.get_converged_reason() << std::endl;

		perf_log.push("KSP solver - B","Coupled CG iterations");
		KSP_Solver_M_B.apply_ZMiZt(vec_pi_old,vec_aux_B);
		perf_log.pop("KSP solver - B","Coupled CG iterations");

		timing_data = perf_log.get_perf_data("KSP solver - B","Coupled CG iterations");
		std::cout << "|        Solver B time : " << timing_data.tot_time/(iter_nb + 1) << std::endl;
		std::cout << "|          Converged ? : " << KSP_Solver_M_B.get_converged_reason() << std::endl;

		perf_log.push("New correction","Coupled CG iterations");

		// q_i = w_i,A + w_i,B
		vec_qi_old =  vec_aux_A;
		vec_qi_old += vec_aux_B;

		// alpha_i = rho_i / ( p_i * q_i )
		aux_double = vec_pi_old.dot(vec_qi_old);
		alpha_i_old = rho_i_old / aux_double;
//		aux_double = 0;

		// lambda_(i+1) = lambda_i + alpha_i * p_i
		vec_lambda = vec_lambda_old;
		vec_lambda.add(alpha_i_old,vec_pi_old);

		perf_log.pop("New correction","Coupled CG iterations");

		// -> Update auxiliary vectors
		perf_log.push("Update auxiliary vectors","Coupled CG iterations");

		// r_(i+1) = r_i - alpha_i * p_i
		vec_ri = vec_ri_old;
		vec_ri.add(-alpha_i_old,vec_pi_old);;

		// z_(i+1) = PC * r_(i+1)
		m_PC->vector_mult(vec_zi,vec_ri);

		// rho_(i+1) = r_(i+1) * z_(i+1)
		rho_i = vec_ri.dot(vec_zi);

		// beta_(i+1) = rho_(i+1) / rho_i
		beta_i = rho_i / rho_i_old;

		// p_temp_(i+1) = z_(i+1) + beta_(i+1) * p_i
		vec_pi_temp = vec_zi;
		vec_pi_temp.add(beta_i,vec_pi_old);

		// p_(i+1) = PI_null * p_temp_(i+1)
		MatMult(m_null_PI,vec_pi_temp.vec(),vec_pi.vec());

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
		vec_pi_old = vec_pi;
		vec_lambda_old = vec_lambda;
		vec_ri_old = vec_ri;

		std::cout << "|     Iteration no. " << iter_nb << std::endl;
		std::cout << "|        Lambda (L2) :" << vec_lambda_old.l2_norm() << std::endl;
		std::cout << "|        r (L2)      :" << vec_ri_old.l2_norm() << std::endl;
		std::cout << "|        z (L2)      :" << vec_zi.l2_norm() << std::endl;
		std::cout << "|        p (L2)      :" << vec_pi_old.l2_norm() << std::endl;
		std::cout << "|        q (L2)      :" << vec_qi_old.l2_norm() << std::endl;
		std::cout << "|        rho         :" << rho_i_old << std::endl;
		std::cout << "|        alpha       :" << alpha_i_old << std::endl;
		std::cout << "|        beta        :" << beta_i << std::endl;
		std::cout << "|        aux         :" << aux_double << std::endl;
		std::cout << "|" << std::endl;

		if(m_bPrintRestart)
		{
			std::cout << "|        Writing to files " << m_lambda_i_filename << ", etc ..." << std::endl;
			write_PETSC_vector(vec_pi_old  ,m_p_i_filename);
			write_PETSC_vector(vec_ri_old  ,m_r_i_filename);
			write_PETSC_vector(vec_lambda_old,m_lambda_i_filename);

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
	KSP_Solver_M_A.apply_MiZt(vec_rhs_A,vec_aux_A);
	*m_sol_A.get() = vec_u0_A;
	m_sol_A->add(-1,vec_aux_A);
	perf_log.pop("KSP solver - A","Solution");

	perf_log.push("KSP solver - B","Solution");
	KSP_Solver_M_B.apply_MiZt(vec_rhs_B,vec_aux_B);
	*m_sol_B.get() = vec_u0_B;
	m_sol_B->add(1,vec_aux_B);
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
