/*
 * base_CG_solver.cpp
 *
 *  Created on: Dec 21, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "base_CG_solver.h"

namespace carl
{

void base_CG_solver::set_sol_vectors()
{
	m_sol = std::unique_ptr<libMesh::PetscVector<libMesh::Number> >(new libMesh::PetscVector<libMesh::Number>(*m_comm));
	m_sol->init(m_sys_N,m_sys_local_N);
	m_sol->zero();
	m_sol->close();

	m_initial_sol = std::unique_ptr<libMesh::PetscVector<libMesh::Number> >(new libMesh::PetscVector<libMesh::Number>(*m_comm));
	m_initial_sol->init(m_sys_N,m_sys_local_N);
	m_initial_sol->zero();
	m_initial_sol->close();
};

void base_CG_solver::apply_M(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	m_sys_mat->vector_mult(v_out,v_in);
};

void base_CG_solver::apply_LATIN_operator(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_error_msg("Not tested yet!\n");

	// LATIN operator: v_out = ( M + k * ( C^t * C_R^-1 * C ) )* v_in

	// v_out = k * ( C^t * C_R^-1 * C ) * v_in
	m_solver_A->apply_ZMiZt(v_in,v_out);
	v_out.scale(m_search_k);

	// v_out += M * v_in
	m_sys_mat->vector_mult_add(v_out,v_in);
};

void base_CG_solver::apply_CG_operator(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out)
{
	homemade_error_msg("Not tested yet!\n");

	// CG operator: v_out = ( C1 * M1^-1 * C1^t + C2 * M2^-1 * C2^t ) * v_in

	// temp_sol_I = ( CI * MI^-1 * CI^t ) * v_in
	m_solver_A->apply_MiZt(v_in,*m_temp_sol_A);
	m_solver_B->apply_MiZt(v_in,*m_temp_sol_B);

	// v_out = temp_sol_A + temp_sol_B
	v_out = *m_temp_sol_A;
	v_out += *m_temp_sol_B;
};

void base_CG_solver::build_identity_precond()
{
	m_M_PC = new libMesh::PetscMatrix<libMesh::Number>(*m_comm);
	m_M_PC->init(m_sys_N,m_sys_N,m_sys_local_N,m_sys_local_N);
	m_M_PC->zero();
	m_M_PC->close();
	MatShift(m_M_PC->mat(),1);

	m_bBuiltExplicitPreconditioner  = true;
};

void base_CG_solver::set_matrix_precond()
{
	homemade_error_msg("Not implemented yet!\n");

	// We can use PETSc's preconditioners here
};

void base_CG_solver::set_LATIN_precond()
{
	homemade_error_msg("Not implemented yet!\n");

	// Maybe using simply "m_sys_mat"'s preconditioner is enough?
};

void base_CG_solver::set_CG_precond()
{
	homemade_error_msg("Not implemented yet!\n");

	// Nadia's article speaks about using the inverse of the coupling
	// operator K_C as the preconditioner, but it only defines K_1 and K_2.
	// Is it the C_R matrix we've used for the LATIN method?
};

void base_CG_solver::use_preconditioner(bool bUsePreconditioner)
{
	// Preconditioner to be implemented later
	homemade_error_msg("Preconditioners not implemented yet!\n");
	m_bUsePreconditioner = bUsePreconditioner;
};

void base_CG_solver::set_solver_matrix(libMesh::PetscMatrix<libMesh::Number>& sys_mat_in)
{
	// Set the internal matrix
	m_sys_mat = &sys_mat_in;
	apply_system_matrix = &base_CG_solver::apply_M;

	// Check and set the dimensions
	homemade_assert_msg(m_sys_mat->n() == m_sys_mat->m(), "System matrix must be square!\n");
	m_sys_N = m_sys_mat->m();
	int silly_local;
	MatGetLocalSize(m_sys_mat->mat(),&silly_local,NULL);
	m_sys_local_N = silly_local;

	this->set_sol_vectors();

	// Preconditioner to be implemented later
	{
		this->build_identity_precond();
	}

	m_bSystemOperatorSet = true;
};

void base_CG_solver::set_solver_LATIN(generic_solver_interface& solver_correction, libMesh::PetscVector<libMesh::Number>& sys_mat_in, double search_k_in)
{
	// Set the coupling correction solver
	m_solver_A = &solver_correction;
	apply_system_matrix = &base_CG_solver::apply_LATIN_operator;

	// Check and set the system dimensions
	homemade_assert_msg(m_sys_mat->n() == m_sys_mat->m(), "System matrix must be square!\n");
	m_sys_N = m_sys_mat->m();
	int silly_local;
	MatGetLocalSize(m_sys_mat->mat(),&silly_local,NULL);
	m_sys_local_N = silly_local;

	// Set the coupling dimensions
	m_solver_A->get_coupling_dimensions(m_coupling_M,m_coupling_N,m_coupling_local_M,m_coupling_local_N);

	// Set the search constant
	m_search_k = search_k_in;

	this->set_sol_vectors();

	// Still have to think about a good preconditioner for the LATIN case
	{
		this->build_identity_precond();
	}

	m_bSystemOperatorSet = true;
};

void base_CG_solver::set_solver_CG(generic_solver_interface& solver_in_A, generic_solver_interface& solver_in_B)
{
	// Set the internal solvers
	m_solver_A = &solver_in_A;
	m_solver_B = &solver_in_B;
	apply_system_matrix = &base_CG_solver::apply_CG_operator;

	// Set the dimensions - equal to the coupling matrices' nb. of rows
	m_solver_A->get_coupling_dimensions(m_coupling_M,m_coupling_N,m_coupling_local_N,m_coupling_local_N);
	m_sys_N = m_coupling_M;
	m_sys_local_N = m_coupling_local_M;

	m_temp_sol_A = std::unique_ptr<libMesh::PetscVector<libMesh::Number> >(new libMesh::PetscVector<libMesh::Number>(*m_comm));
	m_temp_sol_A ->init(m_sys_N,m_sys_local_N);
	m_temp_sol_A ->zero();
	m_temp_sol_A ->close();

	m_temp_sol_B = std::unique_ptr<libMesh::PetscVector<libMesh::Number> >(new libMesh::PetscVector<libMesh::Number>(*m_comm));
	m_temp_sol_B ->init(m_sys_N,m_sys_local_N);
	m_temp_sol_B ->zero();
	m_temp_sol_B ->close();

	this->set_sol_vectors();

	// Preconditioner to be implemented later
	{
		this->build_identity_precond();
	}

	m_bSystemOperatorSet = true;
};

void base_CG_solver::set_system_rhs(libMesh::PetscVector<libMesh::Number>& rhs_in)
{
	m_rhs = &rhs_in;
	m_brhsSet = true;
};

void base_CG_solver::solve()
{
	homemade_assert_msg(m_bSystemOperatorSet,"System operator must be set before solving!\n");
	homemade_assert_msg(m_brhsSet,"Right-hand side must be set before solving!\n");

	// Create the iteration vectors
	libMesh::PetscVector<libMesh::Number> m_p(*m_comm), m_p_prev(*m_comm);
	libMesh::PetscVector<libMesh::Number> m_q_prev(*m_comm);
	libMesh::PetscVector<libMesh::Number> m_r(*m_comm), m_r_prev(*m_comm);
	libMesh::PetscVector<libMesh::Number> m_z(*m_comm);
	libMesh::PetscVector<libMesh::Number> m_x(*m_comm), m_x_prev(*m_comm);

	libMesh::PetscVector<libMesh::Number> m_aux(*m_comm);

	m_p.init(m_sys_N,m_sys_local_N);
	m_p.zero();
	m_p.close();

	m_p_prev.init(m_p); m_p_prev.close();
	m_q_prev.init(m_p); m_q_prev.close();
	m_r.init(m_p);      m_r.close();
	m_r_prev.init(m_p); m_r_prev.close();
	m_z.init(m_p);      m_z.close();
	m_x.init(m_p);      m_x.close();
	m_x_prev.init(m_p); m_x_prev.close();
	m_aux.init(m_p);    m_aux.close();
	m_x_prev = *m_initial_sol;

	// Set the iteration search parameters
	double m_rho = 0, m_beta = 0;
	double m_rho_prev = 0, m_alpha_prev = 0;
	double m_rho_zero = 0;
	double aux_double = 0;

	std::cout << "|     Finished setup " << std::endl;

	// Initialize the system
	// r(0) = b - A * x(0)
	m_r_prev = *m_rhs;
	(this->*apply_system_matrix)(m_aux,m_x_prev);
	m_r_prev.add(-1,m_aux);

	// z(0) = M_PC * r(0)
	m_M_PC->vector_mult(m_z,m_r_prev);

	// rho(0) = r(0) . z(0)
	m_rho_prev = m_r_prev.dot(m_z);
	m_rho_zero = m_rho_prev;

	// p(0) = M_proj * z(0)?
	// p(0) = z(0) ?
	if(m_bUseNullSpaceProjector)
	{
		std::cout << "|     Using the projector ... " << std::endl;
		m_M_null_proj->vector_mult(m_p_prev,m_z);
	}
	else
	{
		std::cout << "|     NOT using the projector ... " << std::endl;
		m_p_prev = m_z;
	}

	std::cout << "|     Finished preamble: " << std::endl;

	// Iteration parameters
	unsigned int kkk = 0;
	bool bKeepIterating = true;
	bool bConverged = false;
	bool bDiverged = false;

	std::cout << "|" << std::endl;
	std::cout << "|        rho(0)        :" << m_rho_prev << std::endl;
	std::cout << "|" << std::endl;

	while(bKeepIterating)
	{
		// q(k) = A * p(k)
		(this->*apply_system_matrix)(m_p_prev,m_q_prev);

		// aux_double = p(k) . q(k) = p(k) * A * p(k)
		aux_double = m_p_prev.dot(m_q_prev);

		// alpha(k) = r(k) . z(k) / ( p(k) * A * p(k) )
		//          = rho(k)      / aux_double
		m_alpha_prev = m_rho_prev / aux_double;

		// x(k + 1) = x(k) + alpha(k) * p(k)
		m_x = m_x_prev; m_x.add(m_alpha_prev,m_p_prev);

		// r(k + 1) = r(k) - alpha(k) * A * p(k)
		//          = r(k) - alpha(k) * q(k)
		m_r = m_r_prev; m_r.add(-m_alpha_prev,m_q_prev);

		// z(k + 1) = M_PC * r(k + 1)
		m_M_PC->vector_mult(m_z,m_r);

		// rho(k + 1) = r(k + 1) . z(k + 1)
		m_rho = m_r.dot(m_z);

		// beta(k + 1) = rho(k + 1) / rho(k)
		m_beta = m_rho / m_rho_prev;

		// p(k + 1) = M_proj * ( z(k + 1) + beta(k + 1) * p(k) ) ?
		// p(k + 1) = z(k + 1) + beta(k + 1) * p(k) ?
		m_aux = m_z;
		m_aux.add(m_beta,m_p_prev);

		if(m_bUseNullSpaceProjector)
		{
			m_M_null_proj->vector_mult(m_p,m_aux);
		}
		else
		{
			m_p = m_aux;
		}

		std::cout << "|" << std::endl;
		std::cout << "|     Iteration no. " << kkk << std::endl;
		std::cout << "|        rho(k)        :" << m_rho_prev << std::endl;
		std::cout << "|        alpha(k)      :" << m_alpha_prev << std::endl;
		std::cout << "|        rho(k + 1)    :" << m_rho << std::endl;
		std::cout << "|        beta(k + 1)   :" << m_beta << std::endl;
		std::cout << "|" << std::endl;

		// Advance iteration
		++kkk;

		// Check convergence
		bConverged = test_convergence(kkk,m_rho,m_rho_zero);
		bDiverged = test_divergence(kkk,m_rho,m_rho_zero);

		if(bConverged || bDiverged)
		{
			bKeepIterating = false;
		}
		else
		{
			m_p_prev = m_p;
			m_rho_prev = m_rho;
			m_x_prev = m_x;
			m_r_prev = m_r;
		}
	}

	*m_sol = m_x;

	// Set solution!
	if(bConverged)
	{
		std::cout << "| Converged after " << kkk << " iterations, residual = " << m_rho << std::endl;
	}
	if(bDiverged)
	{
		std::cout << "| DIVERGED after " << kkk << " iterations, residual = " << m_rho << std::endl;
	}
};

libMesh::PetscVector<libMesh::Number>& base_CG_solver::get_solution()
{
	return *m_sol;
};

bool base_CG_solver::test_convergence(unsigned int iter, double res_norm, double init_res_norm)
{
	if(res_norm < m_CG_conv_eps_abs || 		         // Absolute convergence
	   res_norm < m_CG_conv_eps_rel * init_res_norm) // Relative convergence
	{
		return true;
	}

	return false;
};

bool base_CG_solver::test_divergence(unsigned int iter, double res_norm, double init_res_norm)
{
	if(iter > m_CG_conv_max_n || 					 // Iteration divergence
	   res_norm > m_CG_div_tol * init_res_norm)      // Residual divergence
	{
		return true;
	}

	return false;
};

} /* namespace CArl */
