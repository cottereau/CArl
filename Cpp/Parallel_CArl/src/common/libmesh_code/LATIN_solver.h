/*
 * LATIN_solver.h
 *
 *  Created on: Jan 29, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_LATIN_SOLVER_H_
#define COMMON_LIBMESH_CODE_LATIN_SOLVER_H_

#include "carl_headers.h"
#include "coupled_solver.h"
#include "assemble_functions_nonlinear_elasticity_3D.h"

#include "PETSC_matrix_operations.h"

const bool MASTER_bPerfLog_LATIN_solver_matrix_assemble = true;
const bool MASTER_bPerfLog_LATIN_solver_solve = true;

namespace carl
{

class PETSC_LATIN_solver : public coupled_solver
{
protected:

	// Matrices
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_invC_RR_vec;

	std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> > m_H_A;
	std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> > m_H_B;
	Mat m_PETSC_H_A;
	Mat m_PETSC_H_B;

	std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> > m_P_A;
	std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> > m_P_B;
	Mat m_PETSC_P_A;
	Mat m_PETSC_P_B;

	std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> > m_Extra_M_B;
	Mat m_PETSC_Extra_M_B;

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

	// Restart parameters
	std::string m_phi_A_filename;
	std::string m_phi_B_filename;

	// Partitioning debug parameters
	std::string m_info_matrix_H_A_filename;
	std::string m_info_matrix_H_B_filename;

	std::string m_matrix_H_A_filename;
	std::string m_matrix_H_B_filename;

private:
	PETSC_LATIN_solver();

public:

	// Constructors
	PETSC_LATIN_solver(	const libMesh::Parallel::Communicator& comm,
						carl::CoupledSolverType solver_type = carl::LATIN_MODIFIED_STIFFNESS) :
							coupled_solver (comm,solver_type),

							m_LATIN_relax { 0.8 },
							m_LATIN_conv_eps { 1E-2 },
							m_LATIN_conv_max_n { 10000 },
							m_LATIN_conv_n { 0 },

							m_bUseLumping { true }
	{
		m_LATIN_Index.resize(m_LATIN_conv_max_n);
	};

	PETSC_LATIN_solver(	double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB,
						const libMesh::Parallel::Communicator& comm,
						carl::CoupledSolverType solver_type = carl::LATIN_MODIFIED_STIFFNESS)  :
							coupled_solver (comm,solver_type),

							m_k_dA { i_k_dA },
							m_k_dB { i_k_dB },
							m_k_cA { i_k_cA },
							m_k_cB { i_k_cB },

							m_LATIN_relax { 0.8 },
							m_LATIN_conv_eps { 1E-2 },
							m_LATIN_conv_max_n { 10000 },
							m_LATIN_conv_n { 0 },

							m_bUseLumping { true }
	{
		m_LATIN_Index.resize(m_LATIN_conv_max_n);
	};

	// Methods - reimplemented from "coupled_solver.h"
	void set_restart( 	bool bUseRestart,
							bool bPrintRestart,
							const std::string& restart_base_filename);

	void set_info(	bool bSavePartitionInfo,
						const std::string& info_base_filename);

	void solve();

	// Methods
	void use_exact_inverse();

	void use_lumped_inverse();

	void set_params(double i_k_dA, double i_k_dB, double i_k_cA, double i_k_cB);

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

	void set_forces_nonlinear(	libMesh::PetscVector<libMesh::Number>& F_A);

	void set_convergence_limits(double eps, int convIter);

	void set_relaxation(double relax);

	void solve_modified_stiffness();

	void solve_original_stiffness();

	void solve_nonlinear(libMesh::EquationSystems& EqSys_micro, const std::string type_name_micro);

	void print_convergence(std::ostream& convergenceOut);
};

}

#endif /* COMMON_LIBMESH_CODE_LATIN_SOLVER_H_ */
