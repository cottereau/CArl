/*
 * base_CG_solver.h
 *
 *  Created on: Dec 21, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_BASE_CG_SOLVER_H_
#define COMMON_LIBMESH_CODE_BASE_CG_SOLVER_H_

#include "carl_headers.h"

#include "PETSC_matrix_operations.h"
#include "generic_solver_interface.h"

namespace carl
{

class base_CG_solver
{
private:

	// Typedef for the system operator function pointer
	typedef void (base_CG_solver::*sys_fptr)(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);
	typedef void (base_CG_solver::*precond_fptr)(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// libMesh communicator
	const libMesh::Parallel::Communicator * m_comm;

	// Convergence parameters
	double m_CG_conv_eps_abs;
	double m_CG_conv_eps_rel;
	int    m_CG_conv_max_n;
	double m_CG_div_tol;

	// System
	libMesh::PetscVector<libMesh::Number> * m_rhs;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_sol;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_initial_sol;
	libMesh::PetscMatrix<libMesh::Number> * m_M_PC;

	// Function pointers
	sys_fptr apply_system_matrix;
	precond_fptr apply_preconditioner;

	// System and coupling sizes
	unsigned int m_sys_N;
	unsigned int m_sys_local_N;

	unsigned int m_coupling_M;
	unsigned int m_coupling_local_M;

	unsigned int m_coupling_N;
	unsigned int m_coupling_local_N;

	// Preconditioner unique pointer
	std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> > m_coupling_precond_solver;

	// Internal matrix and solver pointers and co.
	libMesh::PetscMatrix<libMesh::Number> * m_sys_mat;
	generic_solver_interface * m_solver_A;
	generic_solver_interface * m_solver_B;

	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_temp_sol_A;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_temp_sol_B;

	double m_search_k;

	// Null space projector
	Mat * m_M_null_proj;

	// Set solution vectors
	void set_sol_vectors();

	// Convergence / divergence tests
	bool test_convergence(unsigned int iter, double res_norm, double init_res_norm);

	bool test_divergence(unsigned int iter, double res_norm, double init_res_norm);

	// Internal system operator functions
	void apply_M(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_LATIN_operator(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_CG_operator(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// Internal
	void apply_precond_matrix(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_coupled_sys_precon(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_inverse_coupling_precond(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// Flags
	bool m_bSystemOperatorSet;
	bool m_brhsSet;
	bool m_bUsePreconditioner;
	bool m_bBuiltExplicitPreconditioner;
	bool m_bPreconditionerSet;
	bool m_bUseNullSpaceProjector;

	// Preconditioner flag
	BaseCGPrecondType m_precond_type;

public:
	base_CG_solver(	const libMesh::Parallel::Communicator& comm ) :
		m_comm { &comm },
		m_CG_conv_eps_abs { 1E-5 },
		m_CG_conv_eps_rel { 1E-20 },
		m_CG_conv_max_n { 1000 },
		m_CG_div_tol { 10000 },
		m_rhs { NULL },
		m_M_PC { NULL },
		apply_system_matrix { NULL },
		apply_preconditioner { NULL },
		m_sys_N { 0 },
		m_sys_local_N { 0 },
		m_coupling_M { 0 },
		m_coupling_local_M { 0 },
		m_coupling_N { 0 },
		m_coupling_local_N { 0 },
		m_sys_mat { NULL },
		m_solver_A { NULL },
		m_solver_B { NULL },
		m_search_k { 0 },
		m_M_null_proj { NULL },
		m_bSystemOperatorSet { false },
		m_brhsSet { false },
		m_bUsePreconditioner { false },
		m_bBuiltExplicitPreconditioner { false },
		m_bPreconditionerSet { false },
		m_bUseNullSpaceProjector { false },
		m_precond_type { BaseCGPrecondType::NO_PRECONDITIONER }
	{
	};

	~base_CG_solver()
	{
		if(m_bBuiltExplicitPreconditioner)
		{
			delete m_M_PC;
		}
	};

	// Set up parameters
	void set_convergence_limits(	double conv_eps_abs_in = 1E-5,
										double conv_eps_rel_in = 1E-20,
										int conv_max_n_in = 1000,
										double div_tol_in = 10000)
	{
		m_CG_conv_eps_abs = conv_eps_abs_in;
		m_CG_conv_eps_rel = conv_eps_rel_in;
		m_CG_conv_max_n = conv_max_n_in;
		m_CG_div_tol = div_tol_in;
	}

	// Preconditioners
	void set_precond_matrix(libMesh::PetscMatrix<libMesh::Number>& m_in);

	void set_inverse_precond(libMesh::PetscMatrix<libMesh::Number>& m_in);

	void set_preconditioner_type(BaseCGPrecondType type_input = BaseCGPrecondType::NO_PRECONDITIONER);

	// Solver setup methods
	void set_solver_CG_projector(Mat& proj_in);

	void set_solver_matrix(libMesh::PetscMatrix<libMesh::Number>& sys_mat_in);

	void set_initial_sol(libMesh::PetscVector<libMesh::Number>& init_sol_in);

	void set_solver_LATIN(generic_solver_interface& solver_correction, libMesh::PetscVector<libMesh::Number>& sys_mat_in, double search_k_in);

	void set_solver_CG(generic_solver_interface& solver_in_A, generic_solver_interface& solver_in_B);

	void set_system_rhs(libMesh::PetscVector<libMesh::Number>& rhs_in);

	// Solve and output methods
	void solve();

	libMesh::PetscVector<libMesh::Number>& get_solution();

	void get_residual_vector(libMesh::PetscVector<libMesh::Number>& vec_out);

};

} /* namespace CArl */
#endif /* COMMON_LIBMESH_CODE_BASE_CG_SOLVER_H_ */

