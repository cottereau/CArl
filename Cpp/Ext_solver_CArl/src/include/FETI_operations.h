/*
 * FETI_operations.h
 *
 *  Created on: Apr 23, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef FETI_OPERATIONS_H_
#define FETI_OPERATIONS_H_

#include "carl_headers.h"
#include "common_enums.h"
#include "PETSC_matrix_operations.h"

namespace carl
{

/**	\brief Class containing the operations needed for the FETI solver.
 *
 *	This class is used by the several `CArl_FETI` programs to do operations of the FETI solver,
 *  including matrix and vector I/O and iteration operations. Due to the need to read vectors and 
 *  matrices and the usage of several PETSc operations for which the libMesh interface was not 
 *  implemented, direct PETSc `Vec`'s and `Mat`'s are used instead of their libMesh interfaces.
 *  libMesh's parallel communicators are still used, though.
 *
 *	To avoid any confusion involving the iteration index `m_kkk`, it is only changed when reading
 *  the previous iterations's scalar data, and is kept constant during an iteration - following the
 *  of Alg. 1 in the article.
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

	//  --- FETI operations declarations
	// Vectors
	Vec		m_u_0_micro;				///< Solution of the decoupled micro model
	Vec		m_u_0_BIG;					///< Solution of the decoupled macro model

	Vec		m_ext_solver_sol_micro;		///< Micro model external solver solution
	Vec		m_ext_solver_sol_BIG;		///< Macro model external solver solution

	Vec		m_current_residual;			///< Current residual vector, `r(kkk+1)` (`r(0)` for the initialization)
	Vec		m_current_z;				///< Current `z(kkk+1)` (`z(0)` for the initialization)
	Vec		m_current_p;				///< Current `p(kkk+1)` (`p(0)` for the initialization)

	Vec		m_current_phi;				///< Current Lagrange multipliers / solution, `phi(kkk+1)` (`phi(0)` for the initialization)
	Vec		m_current_rb_correction;	///< Current rigid body modes correction, `RB_corr(kkk+1)` (`RB_corr(0)` for the initialization)

	Vec     m_previous_residual;		///< Previous residual vector, `r(kkk)`
	Vec		m_previous_phi;				///< Previous Lagrange multipliers / solution, `phi(kkk)`
	Vec *   m_previous_p_ptr;			///< Pointer to the previous `p` vectors, `p(0 ... kkk)`
	Vec *   m_previous_q_ptr;			///< Pointer to the previous `q` vectors, `q(0 ... kkk)`

	// Real values
	double m_gamma;						///< Double containing the value of `rho(kkk) / p(kkk) . q(kkk)`
	double m_rho_0;						///< Double containing the initial `rho(0)`
	double m_current_rho;				///< Double containing the current `rho(kkk+1)`
	double m_current_RB_mode_corr;		///< Double containing the current `RB_mode_corr(kkk+1)`
	double m_previous_rho;				///< Double containing the previous `rho(kkk)`
	double m_previous_RB_mode_corr;		///< Double containing the previous `RB_mode_corr(kkk)`

	// Boolean flags
	bool	m_bSet_u_0;				 		///< Have the u_0 vectors been set?
	bool	m_bSet_ext_solver_sol;	 		///< Have the external solvers solutions been set vectors been set?

	bool	m_bSet_current_residual; 		///< Have the current `r(kkk+1)` vector been set?
	bool	m_bSet_current_z;		 		///< Have the current `z(kkk+1)` vector been set?
	bool    m_bSet_current_p;				///< Have the current `p(kkk+1)` vector been set?
	bool	m_bSet_current_phi;		 		///< Have the current Lagrange multipliers / solution `phi(kkk+1)` been set?
	bool	m_bSet_current_RB_correction;	///< Have the RB modes correction  `RB_corr(kkk+1)` vector been set?

	bool	m_bSet_previous_residual;	///< Have the previous `r(kkk)` vector been set?
	bool	m_bSet_previous_phi;		///< Have the previous Lagrange multipliers / solution `phi(kkk)` been set?
	bool	m_bSet_previous_p_ptr;		///< Have the previous `p` vectors, `p(0 ... kkk)`, been set?
	bool	m_bSet_previous_q_ptr;		///< Have the previous `q` vectors, `q(0 ... kkk)`, been set?

	//  --- Vectors containing the scalar data
	std::vector<double> m_p_dot_q;		///< Vector containing the `p.q` scalar products

	// Boolean flags
	bool 	m_bReadPreviousScalar;		///< Read the previous iterations scalar data?
	bool	m_bCalculatedScalar;		///< Calculated the current iterations scalar data?

	// --- Convergence parameters
	double m_abs_residual_conv;			///< Absolute residual convergence
	double m_rel_residual_conv;			///< Relative residual convergence (relative to initial value)
	double m_rb_modes_conv;				///< Relative RB correction convergence (relative to previous value)
	double m_rel_residual_div;			///< Relative residual divergence (relative to initial value)
 	int m_max_iter_div;					///< Number of iterations divergence

	// Boolean flags
	bool 	m_bConvResidualAbs;			///< The residual converged? (absolute)
	bool 	m_bConvResidualRel;			///< The residual converged? (relative to initial value)
	bool 	m_bConvRBCorrRel;			///< The RB correction converged? (relative to previous value)
	bool 	m_bDivResidualRel;			///< The residual diverged? (relative to initial value)
	bool 	m_bDivResidualNeg;			///< The residual is negative? (more usefull for debugging, really)
	bool 	m_bDivIter;					///< The number of iterations diverged?
	bool    m_bConv;					///< Did the solver converge?
	bool    m_bDiv;						///< Did the solver diverge?

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

	/// PETSc Vec and Mat deallocation, called by the destructor
	void clear_PETSc();

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
		m_kkk { 0 },
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
		m_gamma { -1 },
		m_rho_0 { -1 },
		m_current_rho { -1 },
		m_current_RB_mode_corr { -1 },
		m_previous_rho { -1 },
		m_previous_RB_mode_corr { -1 },
		m_bSet_u_0 { false },
		m_bSet_ext_solver_sol { false },
		m_bSet_current_residual { false },
		m_bSet_current_z { false },
		m_bSet_current_p { false },
		m_bSet_current_phi { false },
		m_bSet_current_RB_correction { false },
		m_bSet_previous_residual { false },
		m_bSet_previous_phi { false },
		m_bSet_previous_p_ptr { false },
		m_bSet_previous_q_ptr { false },
		m_bReadPreviousScalar { false },
		m_bCalculatedScalar { false },
		m_abs_residual_conv { -1 },
		m_rel_residual_conv { -1 },
		m_rb_modes_conv { -1 },
		m_rel_residual_div { -1 },
		m_max_iter_div { -1 },
		m_bConvResidualAbs { false },
		m_bConvResidualRel { false },
		m_bConvRBCorrRel { false },
		m_bDivResidualRel { false },
		m_bDivResidualNeg { false },
		m_bDivIter { false },
		m_bConv { false },
		m_bDiv { false }
	{
	};

	/// Destructor, deallocates the PETSc 
	~FETI_Operations()
	{
		this->clear_PETSc();
	};

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

	/// Read the previous Lagrage multiplier / solution
	void read_previous_phi();

	/// Read the previous residual
	void read_previous_r();

	/// Read all the previous `p` vectors
	void read_all_previous_p();

	/// Read all the previous `q` vectors
	void read_all_previous_q();

	/// Read the scalar data, `rho(0 ... kkk)`, `| RB_corr(0 ... kkk) |` and `p(0 ... kkk - 1).q(0 ... kkk - 1)`
	void read_scalar_data();

	/// Read the vector data - essentially calls the "read_previous" and "read_all_previous" methods
	void read_vector_data();

	//  --- FETI steps methods
	/// Calculate the initial `p(0)` vector
	void calculate_initial_p();

	/// Calculate the inital residual vector, `r(0)`
	void calculate_initial_r();

	/// Calculate the current `p(kkk+1)` vector
	void calculate_p();

	/// Calculate the current `q(kkk)` vector
	void calculate_q();

	/// Calculate the current `r(kkk+1)` residual vector
	void calculate_r();

	/// Calculate the current `z(kkk+1)` vector
	void calculate_z();

	/// Calculate the current `phi(kkk+1)` solution vector
	void calculate_phi();
	
	/// Calculate the rigid body modes corrections
	void calculate_rb_correction();

	/// Calculate the scalar quantities
	void calculate_scalar_data();

	/// Check the convergence
	IterationStatus check_convergence(double rel_residual_conv, double abs_residual_conv, int max_iter_div, double rel_residual_div, double rb_modes_conv = -1);

	/// Increase the iteration counter
	void increase_iter_counter();

	//  --- Write methods
	/// Calculate and export the external solver RHS's for the next iteration
	void export_ext_solver_rhs_iteration();

	/// Calculate and export the external solver RHS's for the first iteration
	void export_ext_solver_rhs_initial();

	/// Calculate and export the external solver RHS's for the decoupled problem
	void export_ext_solver_rhs_decoupled();

	/// Export the current `p(kkk+1)` vector
	void export_p();

	/// Export the current `q(kkk)` vector
	void export_q();

	/// Export the current `r(kkk+1)` residual vector
	void export_r();

	/// Export the current Lagrange multiplier / solution
	void export_phi();

	/// Export the initial iteration vectors, `r(0)` and `p(0)`
	void export_inital_vecs();

	/// Export the initial iteration scalar data, `rho(0)` and `| RB_corr(0) |`
	void export_initial_scalar_data();

	/// Export the iteration scalar data, `rho(kkk+1)`, `| RB_corr(kkk+1) |` and `p(kkk).q(kkk)`
	void export_scalar_data();

	// Export the iteration vectors
	void export_iter_vecs();

	// Print on `std::cout` the current values of the convergence parameters - and if we converged
	void print_previous_iters_conv(int nb_of_iters = 5);
};
}

#endif /* FETI_OPERATIONS_H_ */