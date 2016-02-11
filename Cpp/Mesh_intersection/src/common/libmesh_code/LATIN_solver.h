/*
 * LATIN_solver.h
 *
 *  Created on: Jan 29, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_LATIN_SOLVER_H_
#define COMMON_LIBMESH_CODE_LATIN_SOLVER_H_

#include "common_header.h"
#include "common_header_libmesh.h"
#include "common_functions.h"

#include "PETSC_matrix_operations.h"

const bool MASTER_bPerfLog_LATIN_solver_matrix_assemble = true;
const bool MASTER_bPerfLog_LATIN_solver_solve = true;

namespace carl
{

class PETSC_LATIN_solver
{
protected:

	// Communicator
	const libMesh::Parallel::Communicator * m_comm;

	// Matrices
	libMesh::PetscMatrix<libMesh::Number> * m_C_RA;
	libMesh::PetscMatrix<libMesh::Number> * m_C_RB;
	libMesh::PetscMatrix<libMesh::Number> * m_C_RR;
	libMesh::PetscVector<libMesh::Number> m_invC_RR_vec;

	libMesh::PetscMatrix<libMesh::Number> * m_M_A;
	libMesh::PetscMatrix<libMesh::Number> * m_M_B;

	libMesh::PetscMatrix<libMesh::Number> * m_H_A;
	libMesh::PetscMatrix<libMesh::Number> * m_H_B;
	Mat m_PETSC_H_A;
	Mat m_PETSC_H_B;

	libMesh::PetscMatrix<libMesh::Number> * m_P_A;
	libMesh::PetscMatrix<libMesh::Number> * m_P_B;
	Mat m_PETSC_P_A;
	Mat m_PETSC_P_B;

	// Forces
	libMesh::PetscVector<libMesh::Number> * m_F_A;
	libMesh::PetscVector<libMesh::Number> * m_F_B;

	// Solution
	libMesh::PetscVector<libMesh::Number> m_sol_A;
	libMesh::PetscVector<libMesh::Number> m_sol_B;

	// Constants
	double m_k_dA;
	double m_k_dB;

	double m_k_cA;
	double m_k_cB;

	// Numerical params
	std::vector<double> m_LATIN_Index;
	double m_LATIN_relax;
	double m_LATIN_conv_eps;
	int m_LATIN_conv_max_n;
	int m_LATIN_conv_n;

	double m_KSP_A_eps;
	int m_KSP_A_iter_max;

	double m_KSP_B_eps;
	int m_KSP_B_iter_max;

	// Some utility parameters
	bool m_bUseLumping;
	bool m_bMatricesSetUp;
	bool m_bForcesSetUp;
	bool m_bParamsSetUp;
	bool m_bCheckDimensions;
	bool m_bDeallocateMatrices;

	// System_names
	std::string m_ksp_name_A;
	std::string m_ksp_name_B;

private:
	PETSC_LATIN_solver();

public:

	// Constructors
	PETSC_LATIN_solver(const libMesh::Parallel::Communicator& comm) :
							m_comm { &comm },

							m_LATIN_relax { 0.8 },
							m_LATIN_conv_eps { 1E-2 },
							m_LATIN_conv_max_n { 10000 },
							m_LATIN_conv_n { 0 },

							m_KSP_A_eps { 1E-8 },
							m_KSP_A_iter_max { 10000 },

							m_KSP_B_eps { 1E-8 },
							m_KSP_B_iter_max { 10000 },

							m_bUseLumping { true },
							m_bMatricesSetUp { false },
							m_bForcesSetUp { false },
							m_bParamsSetUp { false },
							m_bCheckDimensions { false },
							m_bDeallocateMatrices { false },

							m_C_RA { NULL },
							m_C_RB { NULL },
							m_C_RR { NULL },
							m_invC_RR_vec { libMesh::PetscVector<libMesh::Number>(comm) },

							m_M_A { NULL },
							m_M_B { NULL },

							m_H_A { NULL },
							m_H_B { NULL },

							m_P_A { NULL },
							m_P_B { NULL },

							m_F_A { NULL },
							m_F_B { NULL },

							m_sol_A { libMesh::PetscVector<libMesh::Number>(comm) },
							m_sol_B { libMesh::PetscVector<libMesh::Number>(comm) },

							m_ksp_name_A { "macro_sys" },
							m_ksp_name_B { "micro_sys" }
	{
		m_LATIN_Index.resize(m_LATIN_conv_max_n);
	};

	PETSC_LATIN_solver(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB, const libMesh::Parallel::Communicator& comm )  :
			m_comm { &comm },

			m_LATIN_relax { 0.8 },
			m_LATIN_conv_eps { 1E-2 },
			m_LATIN_conv_max_n { 10000 },
			m_LATIN_conv_n { 0 },

			m_KSP_A_eps { 1E-8 },
			m_KSP_A_iter_max { 10000 },

			m_KSP_B_eps { 1E-8 },
			m_KSP_B_iter_max { 10000 },

			m_bUseLumping { true },
			m_bMatricesSetUp { false },
			m_bForcesSetUp { false },
			m_bParamsSetUp { true },
			m_bCheckDimensions { false },
			m_bDeallocateMatrices { false },

			m_C_RA { NULL },
			m_C_RB { NULL },
			m_C_RR { NULL },
			m_invC_RR_vec { libMesh::PetscVector<libMesh::Number>(comm) },

			m_M_A { NULL },
			m_M_B { NULL },

			m_H_A { NULL },
			m_H_B { NULL },

			m_P_A { NULL },
			m_P_B { NULL },

			m_F_A { NULL },
			m_F_B { NULL },

			m_sol_A { libMesh::PetscVector<libMesh::Number>(comm) },
			m_sol_B { libMesh::PetscVector<libMesh::Number>(comm) },

			m_ksp_name_A { "macro_sys_" },
			m_ksp_name_B { "micro_sys_" },

			m_k_dA { i_k_dA },
			m_k_dB { i_k_dB },
			m_k_cA { i_k_cA },
			m_k_cB { i_k_cB }
	{
		m_LATIN_Index.resize(m_LATIN_conv_max_n);
	};

	// Destructor
	~PETSC_LATIN_solver()
	{
		if(m_bDeallocateMatrices)
		{
			delete m_P_A;
			m_P_A = NULL;

			delete m_P_B;
			m_P_B = NULL;

			delete m_H_A;
			m_H_A = NULL;

			delete m_H_B;
			m_H_B = NULL;
		}
	};

	// Methods
	void use_exact_inverse();

	void use_lumped_inverse();

	void set_params(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB);

	void set_sys_names(const std::string& name_A, const std::string& name_B);

	void set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
							libMesh::PetscMatrix<libMesh::Number>& M_B,
							libMesh::PetscMatrix<libMesh::Number>& C_RA,
							libMesh::PetscMatrix<libMesh::Number>& C_RB,
							libMesh::PetscMatrix<libMesh::Number>& C_RR,
							double product_prealloc_P_A = 1,
							double product_prealloc_P_B = 1,
							double product_prealloc_H_A = 1500,
							double product_prealloc_H_B = 1500);

	void set_forces(	libMesh::PetscVector<libMesh::Number>& F_A,
						libMesh::PetscVector<libMesh::Number>& F_B);

	void solve();

	void check_dimensions();

	libMesh::PetscVector<libMesh::Number>& get_solution_BIG();

	libMesh::PetscVector<libMesh::Number>& get_solution_micro();

//	void use_exact_inverse()
//	{
//		m_bUseLumping = false;
//	};
//
//	void use_lumped_inverse()
//	{
//		m_bUseLumping = true;
//	};
//
//	void set_params(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB)
//	{
//		m_k_dA = i_k_dA;
//		m_k_dB = i_k_dB;
//		m_k_cA = i_k_cA;
//		m_k_cB = i_k_cB;
//		m_bParamsSetUp = true;
//	};
//
//	void set_sys_names(const std::string& name_A, const std::string& name_B)
//	{
//		m_ksp_name_A = name_A + "_";
//		m_ksp_name_B = name_B + "_";
//	};
//
//	void set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
//						libMesh::PetscMatrix<libMesh::Number>& M_B,
//						libMesh::PetscMatrix<libMesh::Number>& C_RA,
//						libMesh::PetscMatrix<libMesh::Number>& C_RB,
//						libMesh::PetscMatrix<libMesh::Number>& C_RR,
//						double product_prealloc_P_A = 1,
//						double product_prealloc_P_B = 1,
//						double product_prealloc_H_A = 1500,
//						double product_prealloc_H_B = 1500)
//	{
//		libMesh::PerfLog perf_log("Matrix setup",MASTER_bPerfLog_LATIN_solver_matrix_assemble);
//
//		// -> Will need k_d/cI
//		libmesh_assert_msg( m_bParamsSetUp , "LATIN parameters not set up!");
//
//		m_M_A = &M_A;
//		m_M_B = &M_B;
//
//		m_C_RA = &C_RA;
//		m_C_RB = &C_RB;
//		m_C_RR = &C_RR;
//
//
//		// -> Invert the C_RR matrix
//		if(m_bUseLumping)
//		{
//			perf_log.push("Lumping");
//			lump_matrix_and_invert(* m_C_RR,m_invC_RR_vec);
//			perf_log.pop("Lumping");
//		}
//		else
//		{
//			// TODO : implement inverse matrix
//			libmesh_error_msg("   set_matrices : Exact inverse matrix not implemented yet!!!");
//		}
//
//		/*
//		 * 	Since the definition of the PETSc matrix of a libMesh::PetscMatrix
//		 *	object is done at the declaration, we must use a dynamically
//		 *	allocated libMesh::PetscMatrix, and remember to de-allocate it
//		 *	during the destruction.
//		 */
//
//		m_bDeallocateMatrices = true;
//
//		// -> Calculate P_I = invC_RR * C_I
//		perf_log.push("P_I");
//		MatConvert(m_C_RA->mat(),MATSAME,MAT_INITIAL_MATRIX,&m_PETSC_P_A);
//		MatConvert(m_C_RB->mat(),MATSAME,MAT_INITIAL_MATRIX,&m_PETSC_P_B);
//
//		MatDiagonalScale(m_PETSC_P_A,m_invC_RR_vec.vec(),PETSC_NULL);
//		MatDiagonalScale(m_PETSC_P_B,m_invC_RR_vec.vec(),PETSC_NULL);
////		MatMatMult(m_invC_RR.mat(), m_C_RA->mat(), MAT_INITIAL_MATRIX, product_prealloc_P_A, &m_PETSC_P_A);
////		MatMatMult(m_invC_RR.mat(), m_C_RB->mat(), MAT_INITIAL_MATRIX, product_prealloc_P_B, &m_PETSC_P_B);
//
//		m_P_A = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_P_A, *m_comm);
//		m_P_B = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_P_B, *m_comm);
//		perf_log.pop("P_I");
//
//		std::cout << "| P_A " << std::endl;
//		print_matrix_dim(*m_P_A);
//
//		std::cout << "| P_B " << std::endl;
//		print_matrix_dim(*m_P_B);
//
//		// -> Calculate H_I = k_dI * C_I^t * P_I + M_I
//
//		// C_I^t * P_I
//		perf_log.push("Product","H_I");
//		MatTransposeMatMult(m_C_RA->mat(), m_PETSC_P_A, MAT_INITIAL_MATRIX, product_prealloc_H_A, &m_PETSC_H_A);
//		MatTransposeMatMult(m_C_RB->mat(), m_PETSC_P_B, MAT_INITIAL_MATRIX, product_prealloc_H_B, &m_PETSC_H_B);
//		perf_log.pop("Product","H_I");
//
//		// k_dI * ( C_I^t * P_I ) + M_I
//		perf_log.push("Sum","H_I");
//		MatAYPX(m_PETSC_H_A, m_k_dA, m_M_A->mat(), DIFFERENT_NONZERO_PATTERN);
//		MatAYPX(m_PETSC_H_B, m_k_dB, m_M_B->mat(), DIFFERENT_NONZERO_PATTERN);
//		perf_log.pop("Sum","H_I");
//
//		perf_log.push("Build","H_I");
//		m_H_A = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_H_A, *m_comm);
//		m_H_B = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_H_B, *m_comm);
//		perf_log.pop("Build","H_I");
//
//		std::cout << "| H_A " << std::endl;
//		print_matrix_dim(*m_H_A);
//
//		std::cout << "| H_B " << std::endl;
//		print_matrix_dim(*m_H_B);
//
//		m_bMatricesSetUp = true;
//	};
//
//	void set_forces(	libMesh::PetscVector<libMesh::Number>& F_A,
//						libMesh::PetscVector<libMesh::Number>& F_B)
//	{
//		m_F_A = &F_A;
//		m_F_B = &F_B;
//
//		m_bForcesSetUp = true;
//	};
//
//	void solve()
//	{
//		libMesh::PerfLog perf_log("Solve",MASTER_bPerfLog_LATIN_solver_solve);
//
//		// -> Test if the parameters are set up
//		libmesh_assert_msg( m_bParamsSetUp , "   solve : LATIN parameters not set up!");
//		libmesh_assert_msg( m_bMatricesSetUp , "   solve : Matrices not set up!");
//		libmesh_assert_msg( m_bForcesSetUp , "   solve : Forces not set up!");
//
//		// -> Matrix dimensions
//
//		// Check them beforehand
//		this->check_dimensions();
//
//		unsigned int 	dim_A, dim_A_local,
//						dim_B, dim_B_local,
//						dim_R, dim_R_local;
//
//		int silly_local = 0;
//
//		dim_A = m_H_A->m();
//		dim_B = m_H_B->m();
//		dim_R = m_P_A->m();
//
//		MatGetLocalSize(m_H_A->mat(),&silly_local,NULL);
//		dim_A_local = silly_local;
//		MatGetLocalSize(m_H_B->mat(),&silly_local,NULL);
//		dim_B_local = silly_local;
//		MatGetLocalSize(m_P_A->mat(),&silly_local,NULL);
//		dim_R_local = silly_local;
//
//		// -> Declarations
//		double iter_eps = 1 + m_LATIN_conv_eps;
//		int    iter_nb = 0;
//
//		double norm_diff = 0;
//		double norm_sum = 0;
//
//		perf_log.push("Vector declarations","Initialization");
//		m_sol_A.init(dim_A,dim_A_local);
//		m_sol_B.init(dim_B,dim_B_local);
//
//		// Temporary solutions
//		libMesh::PetscVector<libMesh::Number> w_A(*m_comm, dim_R, dim_R_local);
//		libMesh::PetscVector<libMesh::Number> w_B(*m_comm, dim_R, dim_R_local);
//		w_A.zero();
//		w_B.zero();
//
//		libMesh::PetscVector<libMesh::Number> u_A(*m_comm, dim_A, dim_A_local);
//		libMesh::PetscVector<libMesh::Number> u_B(*m_comm, dim_B, dim_B_local);
//		u_A.zero();
//		u_B.zero();
//
//		// Effective forces
//		libMesh::PetscVector<libMesh::Number> f_eff_A(*m_comm, dim_A, dim_A_local);
//		libMesh::PetscVector<libMesh::Number> f_eff_B(*m_comm, dim_B, dim_B_local);
//		f_eff_A.zero();
//		f_eff_B.zero();
//
//		// Auxiliary vectors
//		libMesh::PetscVector<libMesh::Number> aux_A(*m_comm, dim_R, dim_R_local);
//		libMesh::PetscVector<libMesh::Number> aux_B(*m_comm, dim_R, dim_R_local);
//		aux_A.zero();
//		aux_B.zero();
//
//		libMesh::PetscVector<libMesh::Number> phi_cA(*m_comm, dim_R, dim_R_local);
//		libMesh::PetscVector<libMesh::Number> phi_cB(*m_comm, dim_R, dim_R_local);
//		phi_cA.zero();
//		phi_cB.zero();
//
//		libMesh::PetscVector<libMesh::Number> phi_dA(*m_comm, dim_R, dim_R_local);
//		libMesh::PetscVector<libMesh::Number> phi_dB(*m_comm, dim_R, dim_R_local);
//		phi_dA.zero();
//		phi_dB.zero();
//
//		libMesh::PetscVector<libMesh::Number> phi_diff_A(*m_comm, dim_R, dim_R_local);
//		libMesh::PetscVector<libMesh::Number> phi_diff_B(*m_comm, dim_R, dim_R_local);
//		phi_diff_A.zero();
//		phi_diff_B.zero();
//
//		libMesh::PetscVector<libMesh::Number> phi_sum_A(*m_comm, dim_R, dim_R_local);
//		libMesh::PetscVector<libMesh::Number> phi_sum_B(*m_comm, dim_R, dim_R_local);
//		phi_sum_A.zero();
//		phi_sum_B.zero();
//
//		perf_log.pop("Vector declarations","Initialization");
//
//		// Solvers
//		perf_log.push("KSP solvers setup","Initialization");
//		libMesh::PetscLinearSolver<libMesh::Number> KSP_Solver_H_A(*m_comm);
//		libMesh::PetscLinearSolver<libMesh::Number> KSP_Solver_H_B(*m_comm);
//		KSP_Solver_H_A.reuse_preconditioner(true);
//		KSP_Solver_H_B.reuse_preconditioner(true);
//
//		KSP_Solver_H_A.init(m_H_A, m_ksp_name_A.c_str());
//		KSP_Solver_H_B.init(m_H_B, m_ksp_name_B.c_str());
//		perf_log.pop("KSP solvers setup","Initialization");
//
//		// -> Initialize the vectors
//		// u_0,I = H_I^-1 * F_I (KSP SOLVER!)
//
//		perf_log.push("KSP solver","Initialization");
//		std::cout << " init" << std::endl;
//		KSP_Solver_H_A.solve ( *m_H_A, m_sol_A, *m_F_A, m_KSP_A_eps, m_KSP_A_iter_max);
//		std::cout << " >  ";
//		KSP_Solver_H_A.print_converged_reason();
//		KSP_Solver_H_B.solve ( *m_H_B, m_sol_B, *m_F_B, m_KSP_B_eps, m_KSP_B_iter_max);
//		std::cout << " >  ";
//		KSP_Solver_H_B.print_converged_reason();
//		perf_log.pop("KSP solver","Initialization");
//
//		std::cout << std::endl;
//		// w_0,I = P_I * u_0,I
//		m_P_A->vector_mult(w_A,m_sol_A);
//		m_P_B->vector_mult(w_B,m_sol_B);
//
//		// phi_d0,I = - k_dI * w_0,I
//		phi_dA.add(-m_k_dA,w_A);
//		phi_dB.add(-m_k_dB,w_B);
//
//
//		while (iter_eps > m_LATIN_conv_eps && iter_nb < m_LATIN_conv_max_n)
//		{
//			std::cout << " iter no. " << iter_nb << std::endl;
//			// -> Coupled step
//			perf_log.push("Coupled iterations");
//
//			// aux_i,I = k_cI * w_i-1,I - phi_i-1,I
//			aux_A = phi_dA;
//			aux_A.scale(-1);
//			aux_A.add(m_k_cA,w_A);
//
//			aux_B = phi_dB;
//			aux_B.scale(-1);
//			aux_B.add(m_k_cB,w_B);
//
//			// w_i,I = ( aux_i,A + aux_i,B ) / ( k_cA + k_cB )
//			w_A = aux_A;
//			w_A.add(aux_B);
//			w_A.scale( 1./(m_k_cA + m_k_cB) );
//			w_B = w_A;
//
//			// phi_ci,I = k_cI * w_i,I - aux_i,I
//			phi_cA = aux_A;
//			phi_cA.scale(-1);
//			phi_cA.add(m_k_cA,w_A);
//
//			phi_cB = aux_B;
//			phi_cB.scale(-1);
//			phi_cB.add(m_k_cB,w_B);
//
//			perf_log.pop("Coupled - iterations");
//
//			// -> Decoupled step
//
//			// aux_i,I = - k_dI * w_i,I - phi_i,I
//			aux_A = phi_cA;
//			aux_A.scale(-1);
//			aux_A.add( - m_k_dA,w_A);
//
//			aux_B = phi_cB;
//			aux_B.scale(-1);
//			aux_B.add( - m_k_dB,w_B);
//
//			perf_log.push("Effective force","Decoupled - iterations");
//			// f_eff_i,I = F_A - C_I^t aux_i,I
//			MatMultTranspose(m_C_RA->mat(), aux_A.vec(), f_eff_A.vec());
//			f_eff_A.scale(-1);
//			f_eff_A.add(*m_F_A);
//			MatMultTranspose(m_C_RB->mat(), aux_B.vec(), f_eff_B.vec());
//			f_eff_B.scale(-1);
//			f_eff_B.add(*m_F_B);
//			perf_log.pop("Effective force","Decoupled - iterations");
//
//			perf_log.push("KSP solver","Decoupled - iterations");
//			// u_i,I = H_I^-1 * f_eff_i,I (KSP SOLVER!)
//			KSP_Solver_H_A.solve ( *m_H_A, u_A, f_eff_A, m_KSP_A_eps, m_KSP_A_iter_max);
//			std::cout << " >  ";
//			KSP_Solver_H_A.print_converged_reason();
//			KSP_Solver_H_B.solve ( *m_H_B, u_B, f_eff_B, m_KSP_B_eps, m_KSP_B_iter_max);
//			std::cout << " >  ";
//			KSP_Solver_H_B.print_converged_reason();
//			perf_log.pop("KSP solver","Decoupled - iterations");
//
//			// w_i,I = P_I * u_i,I
//			perf_log.push("Projection","Decoupled - iterations");
//			m_P_A->vector_mult(w_A,u_A);
//			m_P_B->vector_mult(w_B,u_B);
//			perf_log.pop("Projection","Decoupled - iterations");
//
//			// phi_i,I = - k_dI * w_i,I - aux_i,I
//			phi_dA = aux_A;
//			phi_dA.scale(-1);
//			phi_dA.add( - m_k_dA,w_A);
//
//			phi_dB = aux_B;
//			phi_dB.scale(-1);
//			phi_dB.add( - m_k_dB,w_B);
//
//			// Relaxation : sol_i,I = LATINRelax*u_i,I + (1 - LATINRelax)*sol_i-1,I;
//			m_sol_A.scale(1 - m_LATIN_relax);
//			m_sol_A.add(m_LATIN_relax,u_A);
//
//			m_sol_B.scale(1 - m_LATIN_relax);
//			m_sol_B.add(m_LATIN_relax,u_B);
//
//			// -> Test convergence over the phi's
//			phi_diff_A = phi_dA;
//			phi_diff_A.add(-1,phi_cA);
//
//			phi_diff_B = phi_dB;
//			phi_diff_B.add(-1,phi_cB);
//
//			phi_sum_A = phi_dA;
//			phi_sum_A.add(phi_cA);
//
//			phi_sum_B = phi_dB;
//			phi_sum_B.add(phi_cB);
//
//			norm_diff = phi_diff_A.l2_norm() + phi_diff_B.l2_norm();
//			norm_sum = phi_sum_A.l2_norm() + phi_sum_B.l2_norm();
//
//			if(norm_sum < 1E-18)
//			{
//				break;
//			}
//
//			iter_eps = 2*norm_diff / norm_sum;
//			m_LATIN_Index[iter_nb] = iter_eps;
//			std::cout << " > eps : " << iter_eps << std::endl << std::endl;
//			++iter_nb;
//		}
//
//		std::cout << "| LATIN solver: " << std::endl;
//		std::cout << "|     nb. of iterations : " << iter_nb;
//		if(iter_nb == m_LATIN_conv_max_n)
//		{
//			std::cout << " (MAX!)";
//		}
//		std::cout << std::endl;
//		std::cout << "|     eps               : " << iter_eps << std::endl << std::endl;
//	}
//
//	void check_dimensions()
//	{
//		int H_A_mmm, H_A_nnn, H_A_local_mmm, H_A_local_nnn,
//			H_B_mmm, H_B_nnn, H_B_local_mmm, H_B_local_nnn,
//
//			C_A_mmm, C_A_nnn, C_A_local_mmm, C_A_local_nnn,
//			C_B_mmm, C_B_nnn, C_B_local_mmm, C_B_local_nnn,
//
//			P_A_mmm, P_A_nnn, P_A_local_mmm, P_A_local_nnn,
//			P_B_mmm, P_B_nnn, P_B_local_mmm, P_B_local_nnn,
//
//			F_A_size, F_A_local_size,
//			F_B_size, F_B_local_size;
//
//
//		H_A_mmm = m_H_A->m(); H_A_nnn = m_H_A->n();
//		MatGetLocalSize(m_H_A->mat(),&H_A_local_mmm,&H_A_local_nnn);
//
//		H_B_mmm = m_H_B->m(); H_B_nnn = m_H_B->n();
//		MatGetLocalSize(m_H_B->mat(),&H_B_local_mmm,&H_B_local_nnn);
//
//		C_A_mmm = m_C_RA->m(); C_A_nnn = m_C_RA->n();
//		MatGetLocalSize(m_C_RA->mat(),&C_A_local_mmm,&C_A_local_nnn);
//
//		C_B_mmm = m_C_RB->m(); C_B_nnn = m_C_RB->n();
//		MatGetLocalSize(m_C_RB->mat(),&C_B_local_mmm,&C_B_local_nnn);
//
//		P_A_mmm = m_P_A->m(); P_A_nnn = m_P_A->n();
//		MatGetLocalSize(m_P_A->mat(),&P_A_local_mmm,&P_A_local_nnn);
//
//		P_B_mmm = m_P_B->m(); P_B_nnn = m_P_B->n();
//		MatGetLocalSize(m_P_B->mat(),&P_B_local_mmm,&P_B_local_nnn);
//
//		F_A_size = m_F_A->size(); F_A_size = m_F_A->local_size();
//		F_B_size = m_F_B->size(); F_B_size = m_F_B->local_size();
//
//		// Test if the matrices are squared
//		libmesh_assert_msg( H_A_mmm == H_A_nnn , "   check_dimensions : H_A.m() != H_A.n() !");
//		libmesh_assert_msg( H_B_mmm == H_B_nnn , "   check_dimensions : H_B.m() != H_B.n() !");
//
//		// Test the projections
//		libmesh_assert_msg( H_A_nnn == P_A_nnn , "   check_dimensions : H_A.n() != P_A.n() !");
//		libmesh_assert_msg( H_B_nnn == P_B_nnn , "   check_dimensions : H_B.n() != P_B.n() !");
//		libmesh_assert_msg( P_A_mmm == P_B_mmm , "   check_dimensions : P_A.m() != P_B.m() !");
//
//		// Test the couplings
//		libmesh_assert_msg( H_A_nnn == C_A_nnn , "   check_dimensions : H_A.n() != C_A.n() !");
//		libmesh_assert_msg( H_B_nnn == C_B_nnn , "   check_dimensions : H_B.n() != C_B.n() !");
//		libmesh_assert_msg( C_A_mmm == C_B_mmm , "   check_dimensions : C_A.m() != C_B.m() !");
//
//		// Test the forces
//		libmesh_assert_msg( H_A_nnn == F_A_size , "   check_dimensions : H_A.n() != F_A.size() !");
//		libmesh_assert_msg( H_B_nnn == F_B_size , "   check_dimensions : H_B.n() != F_B.size() !");
//
//
//		// -> Now local !
//
//		// Test if the matrices are squared
//		libmesh_assert_msg( H_A_local_mmm == H_A_local_nnn , "   check_dimensions : H_A.m() != H_A.n() (local) !");
//		libmesh_assert_msg( H_B_local_mmm == H_B_local_nnn , "   check_dimensions : H_B.m() != H_B.n() (local) !");
//
//		// Test the projections
//		libmesh_assert_msg( H_A_local_nnn == P_A_local_nnn , "   check_dimensions : H_A.n() != P_A.n() (local) !");
//		libmesh_assert_msg( H_B_local_nnn == P_B_local_nnn , "   check_dimensions : H_B.n() != P_B.n() (local) !");
//		libmesh_assert_msg( P_A_local_mmm == P_B_local_mmm , "   check_dimensions : P_A.m() != P_B.m() (local) !");
//
//		// Test the couplings
//		libmesh_assert_msg( H_A_local_nnn == C_A_local_nnn , "   check_dimensions : H_A.n() != C_A.n() (local) !");
//		libmesh_assert_msg( H_B_local_nnn == C_B_local_nnn , "   check_dimensions : H_B.n() != C_B.n() (local) !");
//		libmesh_assert_msg( C_A_local_mmm == C_B_local_mmm , "   check_dimensions : C_A.m() != C_B.m() (local) !");
//
//		// Test the forces
//		libmesh_assert_msg( H_A_local_nnn == F_A_local_size , "   check_dimensions : H_A.n() != F_A.size() (local) !");
//		libmesh_assert_msg( H_B_local_nnn == F_B_local_size , "   check_dimensions : H_B.n() != F_B.size() (local) !");
//
//		m_bCheckDimensions = true;
//	}
//
//	libMesh::PetscVector<libMesh::Number>& get_solution_BIG()
//	{
//		return m_sol_A;
//	}
//
//	libMesh::PetscVector<libMesh::Number>& get_solution_micro()
//	{
//		return m_sol_B;
//	}
};

}

#endif /* COMMON_LIBMESH_CODE_LATIN_SOLVER_H_ */
