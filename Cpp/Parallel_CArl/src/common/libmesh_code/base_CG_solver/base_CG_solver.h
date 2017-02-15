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
	typedef void (base_CG_solver::*set_nul_proj_fptr)(libMesh::PetscMatrix<libMesh::Number>& M_sys, libMesh::PetscMatrix<libMesh::Number>& C_sys);
	typedef void (base_CG_solver::*proj_residual_fptr)(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);
	typedef void (base_CG_solver::*proj_force_fptr)(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);
	typedef void (base_CG_solver::*corr_sol_fptr)(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);
	typedef void (base_CG_solver::*get_corr_fptr)(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// libMesh communicator
	const libMesh::Parallel::Communicator * m_comm;

	// Convergence parameters
	double m_CG_conv_eps_abs;
	double m_CG_conv_eps_rel;
	int    m_CG_conv_max_n;
	double m_CG_div_tol;
	double m_CG_conv_nullspace_corr_rel;

	std::vector<double> m_CG_Index;
	std::vector<double> m_CG_correction_norm;

	int    m_CG_conv_n;

	// System
	libMesh::PetscVector<libMesh::Number> * m_rhs;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_sol;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_initial_sol;
	libMesh::PetscMatrix<libMesh::Number> * m_M_PC;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_M_PC_jacobi;

	// Function pointers
	sys_fptr apply_system_matrix;
	precond_fptr apply_preconditioner;
	set_nul_proj_fptr set_nullspace_matrices;
	proj_force_fptr apply_force_projection;
	proj_residual_fptr apply_residual_projection;
	corr_sol_fptr correct_solution;
	get_corr_fptr get_correction;

	// System and coupling sizes
	unsigned int m_sys_N;
	unsigned int m_sys_local_N;

	unsigned int m_coupling_M;
	unsigned int m_coupling_local_M;

	unsigned int m_coupling_N;
	unsigned int m_coupling_local_N;

	// Preconditioner unique pointer
	std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> > m_coupling_precond_solver;

	// Auxiliar vectors
	libMesh::PetscVector<libMesh::Number> m_aux;
	libMesh::PetscVector<libMesh::Number> m_aux_bis;

	// Internal matrix and solver pointers and co.
	libMesh::PetscMatrix<libMesh::Number> * m_sys_mat;
	generic_solver_interface * m_solver_A;
	generic_solver_interface * m_solver_B;

	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_temp_sol_A;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_temp_sol_B;

	double m_search_k;

	MatNullSpace nullsp_sys;

	// Null space projectors - built
	Mat m_null_F;
	Mat m_null_PI;

	Mat m_null_R;
	Mat m_null_RT;

	Mat m_null_sol_correction;

	Mat RI_mat, RI_T_mat;

	// Null space projectors - runtime
	libMesh::PetscMatrix<libMesh::Number> * m_M_sys;
	libMesh::PetscMatrix<libMesh::Number> * m_C_sys;

	Mat	RITRI_mat, inv_RITRI_mat;

	PetscInt    null_nb_vecs;
	const Vec*	null_vecs;
	Vec*	    null_coupled_vecs;

	Vec			aux_null_vec_input;
	Vec			aux_null_vec_output;

	// Perf log
	libMesh::PerfLog m_perf_log;

	// Set solution vectors
	void set_sol_vectors();

	// Convergence / divergence tests
	bool test_convergence(unsigned int iter, double res_norm, double init_res_norm, double rel_sol_eps);

	bool test_correction_convergence(unsigned int iter, double corr_norm_prev, double corr_norm);

	bool test_divergence(unsigned int iter, double res_norm, double init_res_norm);

	// Internal system operator functions
	void apply_M(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_LATIN_operator(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_CG_operator(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// Internal
	void apply_precond_matrix(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_coupled_sys_precon(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_inverse_coupling_precond(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	void apply_jacobi_precond(libMesh::PetscVector<libMesh::Number>& v_in, libMesh::PetscVector<libMesh::Number>& v_out);

	// Nullspace methods

	// -> Brute force
	void build_CG_null_space_projection_matrices(libMesh::PetscMatrix<libMesh::Number>& M_sys, libMesh::PetscMatrix<libMesh::Number>& C_sys);

	void add_CG_built_nullspace_correction(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void apply_CG_built_nullspace_residual_projection(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void apply_CG_built_nullspace_force_projection(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	// -> Runtime
	void build_CG_runtime_null_space_projection_matrices(libMesh::PetscMatrix<libMesh::Number>& M_sys, libMesh::PetscMatrix<libMesh::Number>& C_sys);

	void get_CG_runtime_nullspace_correction(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void add_CG_runtime_nullspace_correction(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void apply_CG_runtime_nullspace_residual_projection(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void apply_CG_runtime_nullspace_force_projection(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void apply_precond_and_projection(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	// Flags
	bool m_bSystemOperatorSet;
	bool m_brhsSet;

	bool m_bUsePreconditioner;
	bool m_bBuiltExplicitPreconditioner;
	bool m_bPreconditionerSet;

	bool m_bUseNullSpaceProjector;
	bool m_bCreatedRigidBodyProjectors_built;
	bool m_bCreatedRigidBodyProjectors_runtime;

	bool m_bReorthogonalizeCorrections;

	// Preconditioner flag
	BaseCGPrecondType m_precond_type;

public:
	base_CG_solver(	const libMesh::Parallel::Communicator& comm ) :
		m_comm { &comm },
		m_CG_conv_eps_abs { 1E-5 },
		m_CG_conv_eps_rel { 1E-20 },
		m_CG_conv_max_n { 1000 },
		m_CG_div_tol { 10000 },
		m_CG_conv_nullspace_corr_rel { 1E-6 },
		m_CG_conv_n { -1 },
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
		m_aux(*m_comm),
		m_aux_bis(*m_comm),
		m_sys_mat { NULL },
		m_solver_A { NULL },
		m_solver_B { NULL },
		m_search_k { 0 },
		m_M_sys { NULL },
		m_C_sys { NULL },
		m_perf_log("Base CG solver"),
		m_bSystemOperatorSet { false },
		m_brhsSet { false },
		m_bUsePreconditioner { false },
		m_bBuiltExplicitPreconditioner { false },
		m_bPreconditionerSet { false },
		m_bUseNullSpaceProjector { false },
		m_bCreatedRigidBodyProjectors_built { false },
		m_bCreatedRigidBodyProjectors_runtime { false },
		m_bReorthogonalizeCorrections { true },
		m_precond_type { BaseCGPrecondType::NO_PRECONDITIONER }
	{
		set_nullspace_matrices = &base_CG_solver::build_CG_runtime_null_space_projection_matrices;
		apply_force_projection = &base_CG_solver::apply_CG_runtime_nullspace_force_projection;
		apply_residual_projection = &base_CG_solver::apply_CG_runtime_nullspace_residual_projection;
		correct_solution = &base_CG_solver::add_CG_runtime_nullspace_correction;
		get_correction = &base_CG_solver::get_CG_runtime_nullspace_correction;
	};

	~base_CG_solver()
	{
		if(m_bBuiltExplicitPreconditioner)
		{
			delete m_M_PC;
		}
		if(m_bCreatedRigidBodyProjectors_built)
		{
			MatDestroy(&m_null_F);
			MatDestroy(&m_null_PI);
			MatDestroy(&m_null_R);
			MatDestroy(&m_null_RT);
			MatDestroy(&m_null_sol_correction);

			MatDestroy(&RI_mat);
			MatDestroy(&RI_T_mat);
			MatDestroy(&RITRI_mat);
			MatDestroy(&inv_RITRI_mat);
		}
		if(m_bCreatedRigidBodyProjectors_runtime)
		{
			VecDestroyVecs(null_nb_vecs,&null_coupled_vecs);
			MatDestroy(&RITRI_mat);
			MatDestroy(&inv_RITRI_mat);
			VecDestroy(&aux_null_vec_input);
			VecDestroy(&aux_null_vec_output);
		}
	};

	// Set up parameters
	void set_convergence_limits(	double conv_eps_abs_in = 1E-20,
										double conv_eps_rel_in = 1E-5,
										int conv_max_n_in = 1000,
										double div_tol_in = 10000,
										double conv_nullspace_corr_rel_in = 1E-6)
	{
		m_CG_conv_eps_abs = conv_eps_abs_in;
		m_CG_conv_eps_rel = conv_eps_rel_in;
		m_CG_conv_max_n = conv_max_n_in;
		m_CG_div_tol = div_tol_in;
		m_CG_conv_nullspace_corr_rel = conv_nullspace_corr_rel_in;
	}

	// Preconditioners
	void set_precond_matrix(libMesh::PetscMatrix<libMesh::Number>& m_in);

	void set_inverse_precond(libMesh::PetscMatrix<libMesh::Number>& m_in);

	void set_jacobi_precond(libMesh::PetscMatrix<libMesh::Number>& m_in);

	void set_preconditioner_type(BaseCGPrecondType type_input = BaseCGPrecondType::NO_PRECONDITIONER);

	// Nullspace methods
	void set_CG_null_space_projection_matrices(libMesh::PetscMatrix<libMesh::Number>& M_sys, libMesh::PetscMatrix<libMesh::Number>& C_sys);

	void add_CG_nullspace_correction(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void apply_CG_nullspace_residual_projection(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void apply_CG_nullspace_force_projection(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	// Solver setup methods
	void set_solver_matrix(libMesh::PetscMatrix<libMesh::Number>& sys_mat_in);

	void set_initial_sol(libMesh::PetscVector<libMesh::Number>& init_sol_in);

	void set_solver_LATIN(generic_solver_interface& solver_correction, libMesh::PetscMatrix<libMesh::Number>& sys_mat_in, double search_k_in);

	void set_solver_CG(generic_solver_interface& solver_in_A, generic_solver_interface& solver_in_B);

	void set_system_rhs(libMesh::PetscVector<libMesh::Number>& rhs_in);

	// Solve and output methods
	void solve();

	libMesh::PetscVector<libMesh::Number>& get_solution();

	void get_residual_vector(libMesh::PetscVector<libMesh::Number>& vec_out);
	
	void get_residual_vector(libMesh::PetscVector<libMesh::Number>& vec_in, libMesh::PetscVector<libMesh::Number>& vec_out);

	void get_convergence_data(std::vector<double>& CG_Index_output);

	void get_correction_convergence_data(std::vector<double>& CG_correction_norm_output);

	void get_perf_log_timing(double& solve_time, double& precond_time, double& proj_time);
};

} /* namespace CArl */
#endif /* COMMON_LIBMESH_CODE_BASE_CG_SOLVER_H_ */

