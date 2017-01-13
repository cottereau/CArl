/*
 * CG_solver.h
 *
 *  Created on: Jan 29, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_CG_COUPLED_SOLVER_H_
#define COMMON_LIBMESH_CODE_CG_COUPLED_SOLVER_H_

#include "carl_headers.h"
#include "coupled_solver.h"
#include "KSP_linear_solver.h"
#include "assemble_functions_nonlinear_elasticity_3D.h"
#include "base_CG_solver.h"

#include "PETSC_matrix_operations.h"

const bool MASTER_bPerfLog_CG_solver_matrix_assemble = true;
const bool MASTER_bPerfLog_CG_solver_solve = true;

namespace carl
{

class PETSC_CG_coupled_solver : public coupled_solver
{
protected:

	// Preallocator
	std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> > m_PC;

	// Convergence parameters
	double m_CG_conv_eps_abs;
	double m_CG_conv_eps_rel;
	int    m_CG_conv_max_n;
	double m_CG_div_tol;

	std::vector<double> m_CG_Index;
	int    m_CG_conv_n;

	// Flags
	bool m_bCoordsSetup;
	BaseCGPrecondType m_precond_type;

	// Partitioning debug parameters
	std::string m_info_matrix_PC_filename;
	std::string m_matrix_PC_filename;

	// Set up the system solvers
	generic_solver_interface * m_sys_A_solver;
	generic_solver_interface * m_sys_B_solver;
	base_CG_solver m_coupling_solver;

	// Restart parameters
	std::string m_u0_A_filename;
	std::string m_u0_B_filename;
	std::string m_p_i_filename;
	std::string m_r_i_filename;
	std::string m_lambda_i_filename;
	std::string m_rho_filename;

private:
	PETSC_CG_coupled_solver();

public:

	// Constructors
	PETSC_CG_coupled_solver(	const libMesh::Parallel::Communicator& comm) :
							coupled_solver (comm,carl::CG),

//							m_coord_vect_A { NULL },
//							m_coord_vect_B { NULL },
							m_CG_conv_eps_abs { 1E-5 },
							m_CG_conv_eps_rel { 1E-20 },
							m_CG_conv_max_n { 1000 },
							m_CG_div_tol { 10000 },
							m_bCoordsSetup { false },
							m_precond_type { BaseCGPrecondType::NO_PRECONDITIONER },
							m_sys_A_solver { NULL },
							m_sys_B_solver { NULL },
							m_coupling_solver(comm)
	{
		m_CG_Index.resize(m_CG_conv_max_n);
		m_bParamsSetUp = true;
	};

	// Set the coupled system solver
	void set_solvers(generic_solver_interface * solver_A, generic_solver_interface * solver_B);

	// Methods - reimplemented from "coupled_solver.h"
	void set_restart( 		bool bUseRestart,
							bool bPrintRestart,
							const std::string& restart_base_filename);

	void set_info(	bool bSavePartitionInfo,
						const std::string& info_base_filename);

	void solve();

	void set_matrices(	libMesh::PetscMatrix<libMesh::Number>& M_A,
							libMesh::PetscMatrix<libMesh::Number>& M_B,
							libMesh::PetscMatrix<libMesh::Number>& C_RA,
							libMesh::PetscMatrix<libMesh::Number>& C_RB,
							libMesh::PetscMatrix<libMesh::Number>& C_RR);

	// Methods
	void set_convergence_limits(double eps_abs, double eps_rel, int convIter = 1E4, double div_tol = 1E4);

	void set_preconditioner_type(BaseCGPrecondType type_input = BaseCGPrecondType::NO_PRECONDITIONER );

	void build_preconditioner();

	void print_convergence(std::ostream& convergenceOut);

	void set_null_space_projector();

	void set_preconditioner_matrix(libMesh::PetscMatrix<libMesh::Number>& M_precond);
};

}

#endif /* COMMON_LIBMESH_CODE_CG_COUPLED_SOLVER_H_ */
