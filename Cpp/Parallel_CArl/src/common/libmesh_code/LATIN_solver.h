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

#include "assemble_functions_nonlinear_elasticity_3D.h"

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
	libMesh::PetscVector<libMesh::Number> * m_invC_RR_vec;

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

	libMesh::PetscMatrix<libMesh::Number> * m_Extra_M_B;
	Mat m_PETSC_Extra_M_B;

	// Forces
	libMesh::PetscVector<libMesh::Number> * m_F_A;
	libMesh::PetscVector<libMesh::Number> * m_F_B;

	// Solution
	libMesh::PetscVector<libMesh::Number> * m_sol_A;
	libMesh::PetscVector<libMesh::Number> * m_sol_B;

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
	bool m_bDeallocateNonlinear;
	bool m_bSolved;

	// Restart parameters
	bool m_bUseRestart;
	bool m_bPrintRestart;

	std::string m_conv_filename;
	std::string m_phi_A_filename;
	std::string m_phi_B_filename;
	std::string m_sol_A_filename;
	std::string m_sol_B_filename;

	// System_names
	std::string m_ksp_name_A;
	std::string m_ksp_name_B;

	std::vector<std::string>   ksp_solver_table;



private:
	PETSC_LATIN_solver();

public:

	// Constructors
	PETSC_LATIN_solver(const libMesh::Parallel::Communicator& comm) :
							m_comm { &comm },

							m_C_RA { NULL },
							m_C_RB { NULL },
							m_C_RR { NULL },
						//	m_invC_RR_vec { libMesh::PetscVector<libMesh::Number>(comm) },
							m_invC_RR_vec { NULL },
							m_M_A { NULL },
							m_M_B { NULL },

							m_H_A { NULL },
							m_H_B { NULL },

							m_P_A { NULL },
							m_P_B { NULL },

							m_F_A { NULL },
							m_F_B { NULL },

							// m_sol_A { libMesh::PetscVector<libMesh::Number>(comm) },
							// m_sol_B { libMesh::PetscVector<libMesh::Number>(comm) },

							m_sol_A { NULL },
							m_sol_B { NULL },
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
							m_bDeallocateNonlinear { false },
							m_bSolved {false},
							m_bUseRestart {false},

							m_ksp_name_A { "macro_sys" },
							m_ksp_name_B { "micro_sys" }
	{
		m_LATIN_Index.resize(m_LATIN_conv_max_n);

		m_invC_RR_vec = new libMesh::PetscVector<libMesh::Number>(comm);

		m_sol_A = new libMesh::PetscVector<libMesh::Number>(comm);
		m_sol_B = new libMesh::PetscVector<libMesh::Number>(comm);
		ksp_solver_table = {"CG","CGN","CGS","CR","QMR","TCQMR","TFQMR",
				   "BICG","BICGSTAB","MINRES","GMRES","LSQR","JACOBI","SOR_FORWARD",
				   "SOR_BACKWARD","SSOR","RICHARDSON","CHEBYSHEV","SPARSELU","INVALID_SOLVER"};
	};

	PETSC_LATIN_solver(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB, const libMesh::Parallel::Communicator& comm )  :
			m_comm { &comm },

			m_C_RA { NULL },
			m_C_RB { NULL },
			m_C_RR { NULL },
//			m_invC_RR_vec { libMesh::PetscVector<libMesh::Number>(comm) },

			m_invC_RR_vec { NULL },
			m_M_A { NULL },
			m_M_B { NULL },

			m_H_A { NULL },
			m_H_B { NULL },

			m_P_A { NULL },
			m_P_B { NULL },

			m_F_A { NULL },
			m_F_B { NULL },

//			m_sol_A { libMesh::PetscVector<libMesh::Number>(comm) },
//			m_sol_B { libMesh::PetscVector<libMesh::Number>(comm) },

			m_sol_A { NULL },
			m_sol_B { NULL },
			m_k_dA { i_k_dA },
			m_k_dB { i_k_dB },
			m_k_cA { i_k_cA },
			m_k_cB { i_k_cB },

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
			m_bSolved {false},
			m_bUseRestart {false},

			m_ksp_name_A { "macro_sys_" },
			m_ksp_name_B { "micro_sys_" }
	{
		m_LATIN_Index.resize(m_LATIN_conv_max_n);

		m_invC_RR_vec = new libMesh::PetscVector<libMesh::Number>(comm);

		m_sol_A = new libMesh::PetscVector<libMesh::Number>(comm);
		m_sol_B = new libMesh::PetscVector<libMesh::Number>(comm);

		ksp_solver_table = {"CG","CGN","CGS","CR","QMR","TCQMR","TFQMR",
				   "BICG","BICGSTAB","MINRES","GMRES","LSQR","JACOBI","SOR_FORWARD",
				   "SOR_BACKWARD","SSOR","RICHARDSON","CHEBYSHEV","SPARSELU","INVALID_SOLVER"};
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

			delete m_sol_A;
			delete m_sol_B;
			delete m_invC_RR_vec;
			m_sol_A = NULL;
			m_sol_B = NULL;
			m_invC_RR_vec = NULL;
		}
		if(m_bDeallocateNonlinear)
		{
			delete m_Extra_M_B;
			m_Extra_M_B = NULL;
		}
	};

	// Methods
	void use_exact_inverse();

	void use_lumped_inverse();

	void set_params(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB);

	void set_sys_names(const std::string& name_A, const std::string& name_B);

	void set_restart( 	bool bUseRestart,
						bool bPrintRestart,
						const std::string& restart_base_filename);

	void set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
							libMesh::PetscMatrix<libMesh::Number>& M_B,
							libMesh::PetscMatrix<libMesh::Number>& C_RA,
							libMesh::PetscMatrix<libMesh::Number>& C_RB,
							libMesh::PetscMatrix<libMesh::Number>& C_RR,
							double product_prealloc_P_A = 1,
							double product_prealloc_P_B = 1,
							double product_prealloc_H_A = 1500,
							double product_prealloc_H_B = 1500);

	void set_matrices_nonlinear(	libMesh::PetscMatrix<libMesh::Number>& M_A,
							libMesh::PetscMatrix<libMesh::Number>& C_RA,
							libMesh::PetscMatrix<libMesh::Number>& C_RB,
							libMesh::PetscMatrix<libMesh::Number>& C_RR,
							double product_prealloc_P_A = 1,
							double product_prealloc_P_B = 1,
							double product_prealloc_H_A = 1500,
							double product_prealloc_H_B = 1500);

	void set_forces(	libMesh::PetscVector<libMesh::Number>& F_A,
						libMesh::PetscVector<libMesh::Number>& F_B);

	void set_forces_nonlinear(	libMesh::PetscVector<libMesh::Number>& F_A);

	void set_convergence_limits(double eps, int convIter);

	void set_relaxation(double relax);

	void solve();

	void solve_nonlinear(libMesh::EquationSystems& EqSys_micro, const std::string type_name_micro);

	void check_dimensions();

	libMesh::NumericVector<libMesh::Number>& get_solution_BIG();

	libMesh::NumericVector<libMesh::Number>& get_solution_micro();

	void print_convergence(std::ostream& convergenceOut);


};

}

#endif /* COMMON_LIBMESH_CODE_LATIN_SOLVER_H_ */
