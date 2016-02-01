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
	libMesh::PetscMatrix<libMesh::Number> m_invC_RR;

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

	// Some utility parameters
	bool m_bUseLumping;
	bool m_bMatricesSetUp;
	bool m_bForcesSetUp;
	bool m_bParamsSetUp;

	bool m_bDeallocateMatrices;

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

							m_bUseLumping { true },
							m_bMatricesSetUp { false },
							m_bForcesSetUp { false },
							m_bParamsSetUp { false },
							m_bDeallocateMatrices { false },

							m_C_RA { NULL },
							m_C_RB { NULL },
							m_C_RR { NULL },
							m_invC_RR { libMesh::PetscMatrix<libMesh::Number>(comm) },

							m_M_A { NULL },
							m_M_B { NULL },

							m_H_A { NULL },
							m_H_B { NULL },

							m_P_A { NULL },
							m_P_B { NULL },

							m_F_A { NULL },
							m_F_B { NULL },

							m_sol_A { libMesh::PetscVector<libMesh::Number>(comm) },
							m_sol_B { libMesh::PetscVector<libMesh::Number>(comm) }
	{
		m_LATIN_Index.resize(m_LATIN_conv_max_n);

	};

	PETSC_LATIN_solver(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB, const libMesh::Parallel::Communicator& comm )  :
			m_comm { &comm },

			m_LATIN_relax { 0.8 },
			m_LATIN_conv_eps { 1E-2 },
			m_LATIN_conv_max_n { 10000 },
			m_LATIN_conv_n { 0 },

			m_bUseLumping { true },
			m_bMatricesSetUp { false },
			m_bForcesSetUp { false },
			m_bParamsSetUp { true },
			m_bDeallocateMatrices { false },

			m_C_RA { NULL },
			m_C_RB { NULL },
			m_C_RR { NULL },
			m_invC_RR { libMesh::PetscMatrix<libMesh::Number>(comm) },

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

	void use_exact_inverse()
	{
		m_bUseLumping = false;
	};

	void use_lumped_inverse()
	{
		m_bUseLumping = true;
	};

	void set_params(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB)
	{
		m_k_dA = i_k_dA;
		m_k_dB = i_k_dB;
		m_k_cA = i_k_cA;
		m_k_cB = i_k_cB;
		m_bParamsSetUp = true;
	};

	void set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
						libMesh::PetscMatrix<libMesh::Number>& M_B,
						libMesh::PetscMatrix<libMesh::Number>& C_RA,
						libMesh::PetscMatrix<libMesh::Number>& C_RB,
						libMesh::PetscMatrix<libMesh::Number>& C_RR)
	{
		// -> Will need k_d/cI
		libmesh_assert_msg( m_bParamsSetUp , "LATIN parameters not set up!");

		m_M_A = &M_A;
		m_M_B = &M_B;

		m_C_RA = &C_RA;
		m_C_RB = &C_RB;
		m_C_RR = &C_RR;

		// -> Invert the C_RR matrix
		if(m_bUseLumping)
		{
			lump_matrix_and_invert(* m_C_RR,m_invC_RR);
		}
		else
		{
			// TODO : implement inverse matrix
			libmesh_error_msg(" Exact inverse matrix not implemented yet!!!");
		}

		/*
		 * 	Since the definition of the PETSc matrix of a libMesh::PetscMatrix
		 *	object is done at the declaration, we must use a dynamically
		 *	allocated libMesh::PetscMatrix, and remember to de-allocate it
		 *	during the destruction.
		 */

		m_bDeallocateMatrices = true;

		// -> Calculate P_I = invC_RR * C_I
		MatMatMult(m_invC_RR.mat(), m_C_RA->mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &m_PETSC_P_A);
		MatMatMult(m_invC_RR.mat(), m_C_RB->mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &m_PETSC_P_B);

		m_P_A = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_P_A, *m_comm);
		m_P_B = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_P_B, *m_comm);

		std::cout << "| P_A " << std::endl;
		print_matrix(*m_P_A);

		std::cout << "| P_B " << std::endl;
		print_matrix(*m_P_B);

		// -> Calculate H_I = k_dI * C_I^t * P_I + M_I

		// C_I^t * P_I
		MatTransposeMatMult(m_C_RA->mat(), m_PETSC_P_A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &m_PETSC_H_A);
		MatTransposeMatMult(m_C_RB->mat(), m_PETSC_P_B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &m_PETSC_H_B);

		// k_dI * ( C_I^t * P_I ) + M_I
		MatAYPX(m_PETSC_H_A, m_k_dA, m_M_A->mat(), DIFFERENT_NONZERO_PATTERN);
		MatAYPX(m_PETSC_H_B, m_k_dB, m_M_B->mat(), DIFFERENT_NONZERO_PATTERN);

		m_H_A = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_H_A, *m_comm);
		m_H_B = new libMesh::PetscMatrix<libMesh::Number>(m_PETSC_H_B, *m_comm);

		std::cout << "| H_A " << std::endl;
		print_matrix(*m_H_A);

		std::cout << "| H_B " << std::endl;
		print_matrix(*m_H_B);

		m_bMatricesSetUp = true;
	};

	void set_forces(	libMesh::PetscVector<libMesh::Number>& F_A,
						libMesh::PetscVector<libMesh::Number>& F_B)
	{
		m_F_A = &F_A;
		m_F_B = &F_B;

		m_bForcesSetUp = true;
	};

	void solve()
	{
		// Test if the parameters are set up
		libmesh_assert_msg( m_bParamsSetUp , "LATIN parameters not set up!");
		libmesh_assert_msg( m_bMatricesSetUp , "Matrices not set up!");
		libmesh_assert_msg( m_bForcesSetUp , "Forces not set up!");

		// TODO : (finally) implement the solver ...
	}

	libMesh::PetscVector<libMesh::Number>& get_solution_micro()
	{
		return m_sol_A;
	}

	libMesh::PetscVector<libMesh::Number>& get_solution_BIG()
	{
		return m_sol_B;
	}
};

}

#endif /* COMMON_LIBMESH_CODE_LATIN_SOLVER_H_ */
