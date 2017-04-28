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
class	FETI_Operations
{
protected:

	//  --- Miscellaneous declarations
	libMesh::Parallel::Communicator& m_comm;	///< Communicator
	bool        m_bScratchFolderSet;			///< Have the scratch folder been set?
	std::string m_scratch_folder_path;			///< Scratch folder path
	IterationStatus	m_FETI_solver_status;		///< Current FETI / CG solver status

	//  --- Coupling matrix and preconditioner declarations

	// Coupling matrices
	Mat m_C_R_micro;	///< Mediator - micro system coupling matrix
	Mat m_C_R_BIG;		///< Mediator - macro system coupling matrix
	Mat m_C_RR;			///< Mediator - mediator system coupling matrix

	// Matrix dimensions
	PetscInt	m_C_R_micro_M, m_C_R_micro_N, m_C_R_micro_M_local, m_C_R_micro_N_local;
	PetscInt	m_C_R_BIG_M, m_C_R_BIG_N, m_C_R_BIG_M_local, m_C_R_BIG_N_local;
	PetscInt	m_C_RR_M, m_C_RR_M_local;

	// Boolean flags
	bool 	    m_bC_R_micro_MatrixSet;	///< Mediator - micro coupling matrix has been set? 
	bool 	    m_bC_R_BIG_MatrixSet;	///< Mediator - macro coupling matrix has been set? 
	bool 	    m_bC_RR_MatrixSet;		///< Mediator - mediator coupling matrix has been set? 
	bool 	    m_bCouplingMatricesSet;	///< All coupling matrices has been set? 

	BaseCGPrecondType	m_precond_type;	///< Preconditioner type. *Default*: `BaseCGPrecondType::NO_PRECONDITIONER`.

	//  --- Null space / rigid body modes declarations

	RBModesSystem m_RB_modes_system;	///< Which model is associated to the RB modes? *Default*: `RBModesSystem::MICRO`.

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
	Vec		m_u_0_micro;				///< Solution of the decoupled micro model
	Vec		m_u_0_BIG;					///< Solution of the decoupled macro model

	Vec		m_ext_solver_sol_micro;		///< Micro model external solver solution
	Vec		m_ext_solver_sol_BIG;		///< Macro model external solver solution

	Vec		m_current_residual;			///< Current residual vector, r(k)
	Vec		m_current_z;				///< Current z(k)
	Vec		m_current_p;				///< Current p(k)

	Vec		m_current_phi;				///< Current Lagrange multipliers
	Vec		m_current_rb_correction;	///< Current rigid body modes correction

	/// Default constructor, set as protected
	FETI_Operations();


public:
	FETI_Operations(libMesh::Parallel::Communicator& comm, const std::string& scratch_folder_path ) :
		m_comm { comm },
		m_bScratchFolderSet { true },
		m_scratch_folder_path { scratch_folder_path },
		m_FETI_solver_status { IterationStatus::ITERATING },
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
		m_null_nb_vecs { -1 },
		m_null_vecs_N { -1 },
		m_null_vecs_N_local { -1 },
		m_bUsingNullVecs { false },
		m_bNullVecsSet { false },
		m_bNullVecsDimensionsSet { false },
		m_binvRITRIMatSet { false }
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
	};

	//  --- Coupling matrix and preconditioner methods
	/// Read the mediator - micro system coupling matrix
	void set_coupling_matrix_R_micro(const std::string& filename);

	/// Read the mediator - macro system coupling matrix
	void set_coupling_matrix_R_BIG(const std::string& filename);

	/// Read the mediator - mediator system coupling matrix
	void set_coupling_matrix_RR(const std::string& filename);

	/// Read all the coupling matrices 
	void read_coupling_matrices(const std::string& filename_base);

	/// Set the preconditioner type
	void set_preconditioner(BaseCGPrecondType CG_precond_type);

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

	/// Calculate the rigid body modes corrections
	void calculate_rb_correction();

	//  --- Write methods
	/// Export the initial iteration vectors, `r_0`, `z_0` and `p_0`
	void export_inital_vecs();

	/// Calculate and export the external solver RHS's
	void export_ext_solver_rhs();

	/// Export the iteration scalar data, `rho_0` and `| RB_corr_0 |`
	void export_scalar_data();
};
}

#endif /* FETI_OPERATIONS_H_ */