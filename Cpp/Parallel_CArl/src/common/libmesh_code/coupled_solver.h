/*
 * coupled_solver.h
 *
 *  Created on: Nov 3, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_COUPLED_SOLVER_H_
#define COMMON_LIBMESH_CODE_COUPLED_SOLVER_H_

#include "carl_headers.h"

#include "PETSC_matrix_operations.h"
#include "generic_solver_interface.h"
#include "base_CG_solver.h"

namespace carl
{

class coupled_solver
{
protected:

	// Communicator
	const libMesh::Parallel::Communicator * m_comm;

	// Matrices
	libMesh::PetscMatrix<libMesh::Number> * m_C_RA;
	libMesh::PetscMatrix<libMesh::Number> * m_C_RB;
	libMesh::PetscMatrix<libMesh::Number> * m_C_RR;

	libMesh::PetscMatrix<libMesh::Number> * m_M_A;
	libMesh::PetscMatrix<libMesh::Number> * m_M_B;

	// Forces
	libMesh::PetscVector<libMesh::Number> * m_F_A;
	libMesh::PetscVector<libMesh::Number> * m_F_B;

	// Null space projectors
	Mat m_null_F;
	Mat m_null_PI;

	Mat m_null_R;
	Mat m_null_RT;

	Mat 		RI_mat, RI_T_mat;

	// Solution
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_sol_A;
	std::unique_ptr<libMesh::PetscVector<libMesh::Number> > m_sol_B;

	double m_KSP_A_eps;
	int m_KSP_A_iter_max;

	double m_KSP_B_eps;
	int m_KSP_B_iter_max;

	// Some utility parameters
	bool m_bMatricesSetUp;
	bool m_bCheckDimensions;
	bool m_bForcesSetUp;
	bool m_bParamsSetUp;
	bool m_bSolved;
	bool m_bCreatedRigidBodyProjectors;

	// Restart parameters
	bool m_bUseRestart;
	bool m_bPrintRestart;

	std::string m_conv_filename;
	std::string m_sol_A_filename;
	std::string m_sol_B_filename;

	// Partitioning debug parameters
	bool m_bSavePartitionInfo;
	std::string m_info_matrix_M_A_filename;
	std::string m_info_matrix_M_B_filename;

	std::string m_matrix_M_A_filename;
	std::string m_matrix_M_B_filename;

	std::string m_matrix_C_A_filename;
	std::string m_matrix_C_B_filename;

	// System_names
	std::string m_ksp_name_A;
	std::string m_ksp_name_B;

	std::vector<std::string>   ksp_solver_table;
	carl::CoupledSolverType    m_solver_type;

private:
	coupled_solver();

public:
	// Constructors
	coupled_solver(	const libMesh::Parallel::Communicator& comm,
						carl::CoupledSolverType solver_type = carl::LATIN_MODIFIED_STIFFNESS) :
							m_comm { &comm },

							m_C_RA { NULL },
							m_C_RB { NULL },
							m_C_RR { NULL },

							m_M_A { NULL },
							m_M_B { NULL },

							m_F_A { NULL },
							m_F_B { NULL },

							m_KSP_A_eps { 1E-5 },
							m_KSP_A_iter_max { 10000 },

							m_KSP_B_eps { 1E-5 },
							m_KSP_B_iter_max { 10000 },

							m_bMatricesSetUp { false },
							m_bCheckDimensions { false },
							m_bForcesSetUp { false },
							m_bParamsSetUp { false },
							m_bSavePartitionInfo { false },

							m_bSolved {false},
							m_bCreatedRigidBodyProjectors { false },
							m_bUseRestart {false},
							m_bPrintRestart {false},

							m_ksp_name_A { "macro_sys" },
							m_ksp_name_B { "micro_sys" },
							m_solver_type { solver_type }
	{
		m_sol_A = std::unique_ptr<libMesh::PetscVector<libMesh::Number> >(new libMesh::PetscVector<libMesh::Number>(comm));
		m_sol_B = std::unique_ptr<libMesh::PetscVector<libMesh::Number> >(new libMesh::PetscVector<libMesh::Number>(comm));
		ksp_solver_table = {"CG","CGN","CGS","CR","QMR","TCQMR","TFQMR",
				   "BICG","BICGSTAB","MINRES","GMRES","LSQR","JACOBI","SOR_FORWARD",
				   "SOR_BACKWARD","SSOR","RICHARDSON","CHEBYSHEV","SPARSELU","INVALID_SOLVER"};
	};

	~coupled_solver()
	{
		if(m_bCreatedRigidBodyProjectors)
		{
			MatDestroy(&m_null_F);
			MatDestroy(&m_null_PI);
			MatDestroy(&m_null_R);
			MatDestroy(&m_null_RT);

			MatDestroy(&RI_mat);
			MatDestroy(&RI_T_mat);
		}
	}
	// Methods
	void set_sys_names(const std::string& name_A, const std::string& name_B);

	void set_forces(	libMesh::PetscVector<libMesh::Number>& F_A,
						libMesh::PetscVector<libMesh::Number>& F_B);

	void set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
						libMesh::PetscMatrix<libMesh::Number>& M_B,
						libMesh::PetscMatrix<libMesh::Number>& C_RA,
						libMesh::PetscMatrix<libMesh::Number>& C_RB,
						libMesh::PetscMatrix<libMesh::Number>& C_RR);

	libMesh::NumericVector<libMesh::Number>& get_solution_BIG();

	libMesh::NumericVector<libMesh::Number>& get_solution_micro();

	void set_restart( 	bool bUseRestart,
						bool bPrintRestart,
						const std::string& restart_base_filename);

	void set_info(	bool bSavePartitionInfo,
						const std::string& info_base_filename);

	void check_dimensions();

	void build_null_space_projection_matrices(libMesh::PetscMatrix<libMesh::Number>& M_sys, libMesh::PetscMatrix<libMesh::Number>& C_sys);

	// Virtual methods
	virtual void solve() = 0;
};

};

#endif /* COMMON_LIBMESH_CODE_COUPLED_SOLVER_H_ */
