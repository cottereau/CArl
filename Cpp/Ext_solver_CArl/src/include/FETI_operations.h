/*
 * FETI_operations.h
 *
 *  Created on: Apr 23, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef FETI_OPERATIONS_H_
#define FETI_OPERATIONS_H_

#include "carl_headers.h"
#include "PETSC_matrix_operations.h"

namespace carl
{

/**	\brief Class containing the operations needed for the FETI solver.
 *
 *		This class is used by the several `CArl_FETI` programs to do operations of the FETI solver,
 *  including matrix and vector I/O and iteration operations. Due to the need to read vectors and 
 *  matrices and the usage of several PETSc operations for which the libMesh interface was not 
 *  implemented, direct PETSc `Vec`'s and `Mat`'s are used instead of their libMesh interfaces.
 *  libMesh's parallel communicators are still used, though.
 *
 */
class	FETI_Operations
{
protected:

	//  --- Miscellaneous declarations
	libMesh::Parallel::Communicator& m_comm;	///< Communicator
	bool        m_bScratchFolderSet;			///< Have the scratch folder been set?
	bool        m_bCouplingFolderSet;			///< Have the Coupling matrices folder been set?
	std::string m_scratch_folder_path;			///< Scratch folder path
	std::string m_coupling_path_base;			///< Coupling matrices path
	IterationStatus	m_FETI_solver_status;		///< Current FETI / CG solver status
	int         m_kkk;							///< Current iteration index
	bool		m_bIterationSet;				///< Have the current iteration index been set?

	//  --- Coupling matrix and preconditioner declarations

	// Coupling matrices
	Mat m_C_R_micro;	///< Mediator - micro system coupling matrix
	Mat m_C_R_BIG;		///< Mediator - macro system coupling matrix
	Mat m_C_RR;			///< Mediator - mediator system coupling matrix

	// Matrix dimensions
	/*
	 *					 M	   x N
	 * m_C_R_micro     : n_med x n_sys_micro
	 * m_C_R_BIG       : n_med x n_sys_BIG
	 * m_C_RR          : n_med x n_med
	 *
	 */
	PetscInt	m_C_R_micro_M, m_C_R_micro_N, m_C_R_micro_M_local, m_C_R_micro_N_local;
	PetscInt	m_C_R_BIG_M, m_C_R_BIG_N, m_C_R_BIG_M_local, m_C_R_BIG_N_local;
	PetscInt	m_C_RR_M, m_C_RR_M_local;

	// Objects used by the preconditioner
	KSP 		m_coupling_precond_solver; 	///< Preconditioner system solver
	Vec			m_coupling_jacobi_precond_vec;	///< Preconditioner Jacobi vector

	// Boolean flags
	bool 	    m_bC_R_micro_MatrixSet;	///< Mediator - micro coupling matrix has been set? 
	bool 	    m_bC_R_BIG_MatrixSet;	///< Mediator - macro coupling matrix has been set? 
	bool 	    m_bC_RR_MatrixSet;		///< Mediator - mediator coupling matrix has been set? 
	bool 	    m_bCouplingMatricesSet;	///< All coupling matrices has been set? 

	BaseCGPrecondType	m_precond_type;	///< Preconditioner type. *Default*: `BaseCGPrecondType::NO_PRECONDITIONER`.

	//  --- Null space / rigid body modes declarations
	RBModesSystem m_RB_modes_system;	///< Which model is associated to the RB modes? *Default*: `RBModesSystem::MICRO`.

	bool		m_bCreatedPrecondSolver;	///< Have the preconditioner system solver been set?
	bool		m_bCreatedPrecondJacobiVec;	///< Have the preconditioner system solver been set?

	// Vectors
	/*
	 *		Note: using Vec arrays does not pose a memory problem, since a PETSc `Vec` is essentially a 
	 *  pointer in the stack. The `Vec`'s contents arre allocated in the heap when calling functions such as
	 *  `VecCreate` and such.
	 *  
	 */ 
	PetscInt    m_null_nb_vecs;	///< Number of null space vectors
	enum 		{ maxVecLimit = 6 };	///< Hack to allow setting the Vec arrays size
	Vec        	m_null_vecs[maxVecLimit];			///< Null space vectors
	Vec		    m_null_coupled_vecs[maxVecLimit]; 	///< Null space vectors times the coupling matrix
	Vec			m_RB_mode_correction;

	// Null vector dimensions
	PetscInt	m_null_vecs_N, m_null_vecs_N_local;

	// Matrices
	Mat m_RITRI_mat;		///< Matrix \f${R^I}^t \cdot R^I = R^t C^t \cdot C R \f$
	Mat m_inv_RITRI_mat;	///< Matrix \f$ inv({R^I}^t \cdot R^I)\f$, used in some projectors

	// Boolean flags
	bool		m_bUsingNullVecs;				///< Do we need to use the null space vectors?
	bool 	    m_bNullVecsSet;					///< Have the null space vectors been set?
	bool 	    m_bNullVecsDimensionsSet;		///< Have the null space vectors' dimensions been set?
	bool 	    m_binvRITRIMatSet;				///< Have the inv(RI^T * RI) matrix been set?
	bool		m_bRBModeSet;					///< Have the RB modes correction vector been set?

	//  --- FETI operations declarations
	// Vectors
	Vec		m_u_0_micro;				///< Solution of the decoupled micro model
	Vec		m_u_0_BIG;					///< Solution of the decoupled macro model

	Vec		m_ext_solver_sol_micro;		///< Micro model external solver solution
	Vec		m_ext_solver_sol_BIG;		///< Macro model external solver solution

	Vec		m_current_residual;			///< Current residual vector, r(k)
	Vec		m_current_z;				///< Current z(k)
	Vec		m_current_p;				///< Current p(k)

	Vec		m_current_phi;				///< Current Lagrange multipliers
	Vec		m_current_rb_correction;	///< Current rigid body modes correction

	// Boolean flags
	bool		m_bSet_u_0;				///< Have the u_0 vectors been set?
	bool		m_bSet_ext_solver_sol;	///< Have the external solvers solutions been set vectors been set?
	bool		m_bSet_current_residual; ///< Have the current_r vector been set?
	bool		m_bSet_current_z;		 ///< Have the current_z vector been set?

	//  --- Protected methods
	// Preconditioner methods
	void set_inverse_precond_solver();	///< Set up the full inversed coupling matrix preconditioner

	void set_jacobi_precond_vector();	///< Set up the Jacobi coupling matrix preconditioner vector

	void read_jacobi_precond_vector();	///< Read the Jacobi coupling matrix preconditioner vector

	void apply_inverse_coupling_precond(Vec vec_in, Vec vec_out);	///< Apply the full inversed coupling matrix preconditioner

	void apply_jacobi_coupling_precond(Vec vec_in, Vec vec_out);	///< Apply the Jacobi coupling matrix preconditioner vector

	void apply_precond(Vec vec_in, Vec vec_out);	///< Common interface to the preconditionners 

	// Projection methods
	void apply_RB_projection(Vec vec_in, Vec vec_out);	///< Apply the rigid body modes projection operation

	/// Calculate and export the external solver RHS's
	void export_ext_solver_rhs(Vec vec_in);

	/// Default constructor
	FETI_Operations();

public:
	/// Constructor with scratch folder path, coupling matrices base filename, and libMesh communicator
	FETI_Operations(libMesh::Parallel::Communicator& comm, const std::string& scratch_folder_path , const std::string& coupling_path_base) :
		m_comm { comm },
		m_bScratchFolderSet { true },
		m_bCouplingFolderSet { true },
		m_scratch_folder_path { scratch_folder_path },
		m_coupling_path_base { coupling_path_base },
		m_FETI_solver_status { IterationStatus::ITERATING },
		m_kkk { -1 },
		m_bIterationSet { false },
		m_C_R_micro_M { -1},
		m_C_R_micro_N { -1},
		m_C_R_micro_M_local { -1},
		m_C_R_micro_N_local { -1},
		m_C_R_BIG_M { -1},
		m_C_R_BIG_N { -1},
		m_C_R_BIG_M_local { -1},
		m_C_R_BIG_N_local { -1},
		m_C_RR_M { -1},
		m_C_RR_M_local { -1},
		m_bC_R_micro_MatrixSet  { false },
		m_bC_R_BIG_MatrixSet  { false },
		m_bC_RR_MatrixSet  { false },
		m_bCouplingMatricesSet  { false },
		m_precond_type { BaseCGPrecondType::NO_PRECONDITIONER },
		m_RB_modes_system { RBModesSystem::MICRO },
		m_bCreatedPrecondSolver { false },
		m_bCreatedPrecondJacobiVec { false },
		m_null_nb_vecs { -1 },
		m_null_vecs_N { -1 },
		m_null_vecs_N_local { -1 },
		m_bUsingNullVecs { false },
		m_bNullVecsSet { false },
		m_bNullVecsDimensionsSet { false },
		m_binvRITRIMatSet { false },
		m_bRBModeSet { false },
		m_bSet_u_0 { false },
		m_bSet_ext_solver_sol { false },
		m_bSet_current_residual { false }
	{
	};

	~FETI_Operations()
	{
		if(m_bC_R_BIG_MatrixSet)
			MatDestroy(&m_C_R_BIG);
		if(m_bC_R_micro_MatrixSet)
			MatDestroy(&m_C_R_micro);
		if(m_bC_RR_MatrixSet)
			MatDestroy(&m_C_RR);
		if(m_bNullVecsSet)
		{
			for(int iii = 0; iii < m_null_nb_vecs; ++iii)
				VecDestroy(&m_null_vecs[iii]);
		}
		if(m_binvRITRIMatSet)
		{
			MatDestroy(&m_inv_RITRI_mat);
		}
		if(m_bSet_u_0)
		{
			VecDestroy(&m_u_0_BIG);
			VecDestroy(&m_u_0_micro);	
		}
		if(m_bSet_ext_solver_sol)
		{
			VecDestroy(&m_ext_solver_sol_BIG);
			VecDestroy(&m_ext_solver_sol_micro);	
		}
		if(m_bSet_current_residual)
		{
			VecDestroy(&m_current_residual);
		}
		if(m_bCreatedPrecondSolver)
		{
			KSPDestroy(&m_coupling_precond_solver);
		}
		if(m_bCreatedPrecondJacobiVec)
		{
			VecDestroy(&m_coupling_jacobi_precond_vec);
		}
		if(m_bRBModeSet)
		{
			VecDestroy(&m_RB_mode_correction);
		}
	};

	/// Set the current iteration index
	void set_iteration(int kkk);

	//  --- Coupling matrix and preconditioner methods
	/// Read the mediator - micro system coupling matrix
	void set_coupling_matrix_R_micro();

	/// Read the mediator - macro system coupling matrix
	void set_coupling_matrix_R_BIG();

	/// Read the mediator - mediator system coupling matrix
	void set_coupling_matrix_RR();

	/// Read all the coupling matrices 
	void read_coupling_matrices();

	/// Set the preconditioner type. The boolean `bInitialSet` determinates if, in the case
	/// of a COUPLING_JACOBI preconditioner, if it will calculate the needed vector (`true`) 
	/// or if it will read it (`false`).
	void set_preconditioner(BaseCGPrecondType CG_precond_type, bool bInitialSet = true);

	//  --- Null space / rigid body modes methods
	/// Set up the 
	void using_rb_modes(bool bUseRigidBodyModes);

	/// Set and print the null space vectors and matrices
	void set_null_space(const std::string& input_filename_base, int nb_of_vecs);

	/// Read the null space vectors
	void read_null_space_vecs(const std::string& RB_vectors_base, int nb_of_rb_vectors);

	/// Read the \f$ inv({R^I}^t \cdot R^I)\f$ matrix
	void read_null_space_inv_RITRI_mat();

	/// Calculate the inital solution, \f$\phi_0\f$
	void calculate_null_space_phi_0(const std::string& force_path);

	//  --- Read methods
	/// Read the decoupled solutions, \f$ u_0,i\f$
	void read_decoupled_solutions();

	/// Read the latest external solver output
	void read_ext_solver_output();

	//  --- FETI steps methods
	/// Calculate the inital residual vector, `r_0`
	void calculate_initial_r();

	/// Calculate the initial `z_0` = `p_0` vectors
	void calculate_initial_z_and_p();

	/// Calculate the current `z` vector
	void calculate_z();
	
	/// Calculate the rigid body modes corrections
	void calculate_rb_correction();

	//  --- Write methods
	
	/// Calculate and export the external solver RHS's for the next iteration
	void export_ext_solver_rhs_iteration();

	/// Export the initial iteration vectors, `r_0`, `z_0` and `p_0`
	void export_inital_vecs();

	/// Export the iteration scalar data, `rho_0` and `| RB_corr_0 |`
	void export_scalar_data();
};
}

#endif /* FETI_OPERATIONS_H_ */