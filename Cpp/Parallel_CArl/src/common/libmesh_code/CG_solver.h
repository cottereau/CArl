/*
 * CG_solver.h
 *
 *  Created on: Jan 29, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_CG_SOLVER_H_
#define COMMON_LIBMESH_CODE_CG_SOLVER_H_

#include "carl_headers.h"
#include "coupled_solver.h"
#include "KSP_linear_solver.h"
#include "assemble_functions_nonlinear_elasticity_3D.h"

#include "PETSC_matrix_operations.h"

const bool MASTER_bPerfLog_CG_solver_matrix_assemble = true;
const bool MASTER_bPerfLog_CG_solver_solve = true;

namespace carl
{

class PETSC_CG_solver : public coupled_solver
{
protected:

	// Preallocator
	std::unique_ptr<libMesh::PetscMatrix<libMesh::Number> > m_PC;

//	// Coordinates vectors
//	libMesh::PetscVector<libMesh::Number> * m_coord_vect_A;
//	libMesh::PetscVector<libMesh::Number> * m_coord_vect_B;

	// Convergence parameters
	double m_CG_conv_eps_abs;
	double m_CG_conv_eps_rel;
	int    m_CG_conv_max_n;
	double m_CG_div_tol;

	std::vector<double> m_CG_Index;
	int    m_CG_conv_n;

	// Flags
	bool m_bUsePreconditioner;
	bool m_bRecalculatePreconditioner;
	bool m_bCoordsSetup;

	// Partitioning debug parameters
	std::string m_info_matrix_PC_filename;
	std::string m_matrix_PC_filename;

	// Restart parameters
	std::string m_u0_A_filename;
	std::string m_u0_B_filename;
	std::string m_p_i_filename;
	std::string m_r_i_filename;
	std::string m_lambda_i_filename;
	std::string m_rho_filename;

private:
	PETSC_CG_solver();

public:

	// Constructors
	PETSC_CG_solver(	const libMesh::Parallel::Communicator& comm) :
							coupled_solver (comm,carl::CG),

//							m_coord_vect_A { NULL },
//							m_coord_vect_B { NULL },
							m_CG_conv_eps_abs { 1E-5 },
							m_CG_conv_eps_rel { 1E-4 },
							m_CG_conv_max_n { 10000 },
							m_CG_div_tol { 1E4 },
							m_bUsePreconditioner { false },
							m_bCoordsSetup { false }
	{
		m_CG_Index.resize(m_CG_conv_max_n);
		m_bParamsSetUp = true;
	};

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

	void use_preconditioner(bool flag = true);

	void build_preconditioner()
	{
	};

	void print_convergence(std::ostream& convergenceOut);
//
//	void set_coordinates(	libMesh::PetscVector<libMesh::Number>& coord_vect_A,
//						libMesh::PetscVector<libMesh::Number>& coord_vect_B);
};

}

#endif /* COMMON_LIBMESH_CODE_CG_SOLVER_H_ */
