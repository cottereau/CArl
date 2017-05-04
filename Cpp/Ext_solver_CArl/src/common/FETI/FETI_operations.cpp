#include "FETI_operations.h"

namespace carl
{
//  --- Protected methods
void FETI_Operations::set_inverse_precond_solver()
{
	homemade_assert_msg(m_bC_RR_MatrixSet,"Preconditioner matrix not set yet!");

	KSPCreate(m_comm.get(), &m_coupling_precond_solver);
	KSPSetOperators(m_coupling_precond_solver, m_C_RR, m_C_RR);
	KSPSetFromOptions(m_coupling_precond_solver);

	m_bCreatedPrecondSolver = true;
}

void FETI_Operations::set_jacobi_precond_vector()
{
	homemade_assert_msg(m_bC_RR_MatrixSet,"Preconditioner matrix not set yet!");

	// Create and set the vector
	VecCreate(m_comm.get(),&m_coupling_jacobi_precond_vec);
	VecSetSizes(m_coupling_jacobi_precond_vec,m_C_RR_M_local,m_C_RR_M);
	VecSetFromOptions(m_coupling_jacobi_precond_vec);

	// Get the diagonal
	MatGetDiagonal(m_C_RR,m_coupling_jacobi_precond_vec);

	// Calculate the reciprocal
	VecReciprocal(m_coupling_jacobi_precond_vec);

	// Export it
	write_PETSC_vector(m_coupling_jacobi_precond_vec,m_scratch_folder_path + "/precond_Jacobi_vector.petscvec",m_comm.rank(),m_comm.get());
	write_PETSC_vector_MATLAB(m_coupling_jacobi_precond_vec,m_scratch_folder_path + "/precond_Jacobi_vector.m",m_comm.get());

	// Set flag
	m_bCreatedPrecondJacobiVec = true;
}

void FETI_Operations::read_jacobi_precond_vector()
{
	// Create and set the vector
	VecCreate(m_comm.get(),&m_coupling_jacobi_precond_vec);
	VecSetSizes(m_coupling_jacobi_precond_vec,m_C_RR_M_local,m_C_RR_M);

	// Read it
	read_PETSC_vector(m_coupling_jacobi_precond_vec,m_scratch_folder_path + "/precond_Jacobi_vector.petscvec", m_comm.get());

	// Set flag
	m_bCreatedPrecondJacobiVec = true;
}

void FETI_Operations::apply_inverse_coupling_precond(Vec vec_in, Vec vec_out)
{
	homemade_assert_msg(m_bCreatedPrecondSolver,"Preconditioner system not set yet!");
	KSPSolve(m_coupling_precond_solver, vec_in, vec_out);
}

void FETI_Operations::apply_jacobi_coupling_precond(Vec vec_in, Vec vec_out)
{	
	homemade_assert_msg(m_bCreatedPrecondJacobiVec,"Preconditioner vector not set yet!");
	VecPointwiseMult(vec_out, m_coupling_jacobi_precond_vec, vec_in);
}

void FETI_Operations::apply_precond(Vec vec_in, Vec vec_out)
{	
	switch (m_precond_type)
	{
		case BaseCGPrecondType::NO_PRECONDITIONER : 
				// Shouldn't call this function in this case
				homemade_error_msg("No preconditioner to be applied!");
				break;

		case BaseCGPrecondType::COUPLING_OPERATOR :	
				this->apply_inverse_coupling_precond(vec_in, vec_out);
				break;

		case BaseCGPrecondType::COUPLING_JACOBI :
				this->apply_jacobi_coupling_precond(vec_in, vec_out);
				break;
		default :
				// Undefined preconditioner
				homemade_error_msg("Undefined preconditioner");
	}
}

void FETI_Operations::apply_RB_projection(Vec vec_in, Vec vec_out)
{
	homemade_assert_msg(m_bNullVecsSet,"Null space vectors not set yet!");
	homemade_assert_msg(m_binvRITRIMatSet,"Null space matrices not set yet!");

	// vec_out = [ I - RC * (inv_RITRI_mat) * RC^t ] * vec_in

	// Declaration of Vecs with size 'm_null_nb_vecs'
	Vec dummy_seq_vec;
	Vec dummy_seq_vec_bis;
	VecCreateSeq(PETSC_COMM_SELF,m_null_nb_vecs,&dummy_seq_vec);
	VecZeroEntries(dummy_seq_vec);
	VecDuplicate(dummy_seq_vec,&dummy_seq_vec_bis);

	// dummy_seq_vec = RC^t * vec_in
	// -> All the communications are done here!
	PetscScalar *dummy_seq_array;
	VecGetArray(dummy_seq_vec,&dummy_seq_array);
	VecMDot(vec_in,m_null_nb_vecs,m_null_coupled_vecs,dummy_seq_array);
	VecRestoreArray(dummy_seq_vec,&dummy_seq_array);

	// dummy_seq_vec_bis = - inv_RITRI_mat * dummy_seq_vec
	// -> Completely local operation!
	MatMult(m_inv_RITRI_mat,dummy_seq_vec,dummy_seq_vec_bis);
	VecScale(dummy_seq_vec_bis,-1);

	// vec_out = vec_in + sum ( dummy_seq_vec_bis[i] * vec_RC[i])
	// -> This should have no communications at all!
	VecCopy(vec_in,vec_out);
	
	VecGetArray(dummy_seq_vec_bis,&dummy_seq_array);
	VecMAXPY(vec_out,m_null_nb_vecs,dummy_seq_array,m_null_coupled_vecs);
	VecRestoreArray(dummy_seq_vec_bis,&dummy_seq_array);

	// Cleanup
	VecDestroy(&dummy_seq_vec);
	VecDestroy(&dummy_seq_vec_bis);
}

void FETI_Operations::export_ext_solver_rhs(Vec vec_in)
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

	Vec vec_C_micro_t_p_kkk_PETSc;
	VecCreate(m_comm.get(),&vec_C_micro_t_p_kkk_PETSc);
	VecSetSizes(vec_C_micro_t_p_kkk_PETSc,m_C_R_micro_N_local,m_C_R_micro_N);
	VecSetFromOptions(vec_C_micro_t_p_kkk_PETSc);

	Vec vec_C_BIG_t_p_kkk_PETSc;
	VecCreate(m_comm.get(),&vec_C_BIG_t_p_kkk_PETSc);
	VecSetSizes(vec_C_BIG_t_p_kkk_PETSc,m_C_R_BIG_N_local,m_C_R_BIG_N);
	VecSetFromOptions(vec_C_BIG_t_p_kkk_PETSc);

	MatMultTranspose(m_C_R_micro,vec_in,vec_C_micro_t_p_kkk_PETSc);
	MatMultTranspose(m_C_R_BIG,vec_in,vec_C_BIG_t_p_kkk_PETSc);

	write_PETSC_vector(vec_C_BIG_t_p_kkk_PETSc,m_scratch_folder_path + "/ext_solver_A_rhs.petscvec",m_comm.rank(),m_comm.get());
	write_PETSC_vector_MATLAB(vec_C_BIG_t_p_kkk_PETSc,m_scratch_folder_path + "/ext_solver_A_rhs.m",m_comm.get());

	write_PETSC_vector(vec_C_micro_t_p_kkk_PETSc,m_scratch_folder_path + "/ext_solver_B_rhs.petscvec",m_comm.rank(),m_comm.get());
	write_PETSC_vector_MATLAB(vec_C_micro_t_p_kkk_PETSc,m_scratch_folder_path + "/ext_solver_B_rhs.m",m_comm.get());


	VecDestroy(&vec_C_micro_t_p_kkk_PETSc);
	VecDestroy(&vec_C_BIG_t_p_kkk_PETSc);
}

void FETI_Operations::clear_PETSc()
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
	if(m_bSet_current_phi)
	{
		VecDestroy(&m_current_phi);
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
	if(m_bSet_current_RB_correction)
	{
		VecDestroy(&m_current_rb_correction);
	}
	if(m_bSet_previous_residual)
	{
		VecDestroy(&m_previous_residual);
	}
	if(m_bSet_previous_phi)
	{
		VecDestroy(&m_previous_phi);
	}
	if(m_bSet_previous_p_ptr)
	{
		VecDestroyVecs(m_kkk+1,&m_previous_p_ptr);
		delete[] m_previous_p_ptr;
	}
	if(m_bSet_previous_q_ptr)
	{
		VecDestroyVecs(m_kkk+1,&m_previous_q_ptr);
		delete[] m_previous_q_ptr;
	}
}

//  --- Coupling matrix and preconditioner methods
void FETI_Operations::set_coupling_matrix_R_micro()
{
	homemade_assert_msg(m_bCouplingFolderSet,"Common coupling matrix path not set yet!");
	MatCreate(m_comm.get(),&m_C_R_micro);
	read_PETSC_matrix(m_C_R_micro,m_coupling_path_base + "_micro.petscmat",m_comm.get());
	MatGetLocalSize(m_C_R_micro,&m_C_R_micro_M_local,&m_C_R_micro_N_local);
	MatGetSize(m_C_R_micro,&m_C_R_micro_M,&m_C_R_micro_N);
	m_C_RR_M = m_C_R_micro_M; m_C_RR_M_local = m_C_R_micro_M_local;

	if(m_RB_modes_system == RBModesSystem::MICRO)
	{
		m_null_vecs_N = m_C_R_micro_N;
		m_null_vecs_N_local = m_C_R_micro_N_local;
		m_bNullVecsDimensionsSet = true;
	}

	m_bC_R_micro_MatrixSet = true;
	m_bCouplingMatricesSet = m_bC_R_BIG_MatrixSet && m_bC_R_micro_MatrixSet && m_bC_RR_MatrixSet;
}

void FETI_Operations::set_coupling_matrix_R_BIG()
{
	homemade_assert_msg(m_bCouplingFolderSet,"Common coupling matrix path not set yet!");
	MatCreate(m_comm.get(),&m_C_R_BIG);
	read_PETSC_matrix(m_C_R_BIG,m_coupling_path_base + "_macro.petscmat",m_comm.get());
	MatGetLocalSize(m_C_R_BIG,&m_C_R_BIG_M_local,&m_C_R_BIG_N_local);
	MatGetSize(m_C_R_BIG,&m_C_R_BIG_M,&m_C_R_BIG_N);
	m_C_RR_M = m_C_R_BIG_M; m_C_RR_M_local = m_C_R_BIG_M_local;

	if(m_RB_modes_system == RBModesSystem::MACRO)
	{
		m_null_vecs_N = m_C_R_BIG_N;
		m_null_vecs_N_local = m_C_R_BIG_N_local;
		m_bNullVecsDimensionsSet = true;
	}

	m_bC_R_BIG_MatrixSet = true;
	m_bCouplingMatricesSet = m_bC_R_BIG_MatrixSet && m_bC_R_micro_MatrixSet && m_bC_RR_MatrixSet;
}

void FETI_Operations::set_coupling_matrix_RR()
{
	homemade_assert_msg(m_bCouplingFolderSet,"Common coupling matrix path not set yet!");
	MatCreate(m_comm.get(),&m_C_RR);
	read_PETSC_matrix(m_C_RR,m_coupling_path_base + "_mediator.petscmat",m_comm.get());
	MatGetLocalSize(m_C_RR,&m_C_RR_M_local,NULL);
	MatGetSize(m_C_RR,&m_C_RR_M,NULL);
	m_bC_RR_MatrixSet = true;
	m_bCouplingMatricesSet = m_bC_R_BIG_MatrixSet && m_bC_R_micro_MatrixSet && m_bC_RR_MatrixSet;
}

void FETI_Operations::read_coupling_matrices()
{
	homemade_assert_msg(m_bCouplingFolderSet,"Common coupling matrix path not set yet!");
	this->set_coupling_matrix_R_micro();
	this->set_coupling_matrix_R_BIG();
	this->set_coupling_matrix_RR();
}

void FETI_Operations::set_preconditioner(BaseCGPrecondType CG_precond_type, bool bInitialSet)
{
	m_precond_type = CG_precond_type;

	switch (m_precond_type)
	{
		case BaseCGPrecondType::NO_PRECONDITIONER : 
				// Well ... nothing to do ...
				break;

		case BaseCGPrecondType::COUPLING_OPERATOR :	
				// Read the mediator - mediator coupling matrix and set the solver
				this->set_coupling_matrix_RR();
				this->set_inverse_precond_solver();
				m_bC_RR_MatrixSet = true;
				break;

		case BaseCGPrecondType::COUPLING_JACOBI :
				if(bInitialSet)
				{
					// Read the mediator - mediator coupling matrix and build the Jacobi coupling preconditioner vector
					this->set_coupling_matrix_RR();
					this->set_jacobi_precond_vector();
					m_bC_RR_MatrixSet = true;
				} else {
					// Just read the Jacobi coupling preconditioner vector
					this->read_jacobi_precond_vector();
					m_bC_RR_MatrixSet = true;
				}
				break;
	}
}

//  --- Null space / rigid body modes methods
void FETI_Operations::using_rb_modes(bool bUseRigidBodyModes)
{
	m_bUsingNullVecs = bUseRigidBodyModes;
}

void FETI_Operations::set_null_space(const std::string& input_filename_base, int nb_of_vecs)
{
	homemade_assert_msg(m_bNullVecsDimensionsSet,"Null vectors sizes not set yet!");
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

	// Set the null vec arrays
	m_null_nb_vecs = nb_of_vecs;

	// Set the (dummy) coupling matrix pointer
	Mat * coupling_matrix;

	switch (m_RB_modes_system)
	{
		case RBModesSystem::MICRO : 
				coupling_matrix = &m_C_R_micro;
				break;

		case RBModesSystem::MACRO :	
				// TODO: Remove the error message below after the more general code was implemented
				homemade_error_msg("Option not implemented yet!");
				coupling_matrix = &m_C_R_BIG;
				break;
	}

	// Matrix dimensions
	//                       M     x N
	// coupling_matrix     : n_med x n_sys
	// m_null_vecs         : n_sys x nb_of_vecs	( R )  -> nb_of_vecs vectors of dim n_sys
	// m_null_coupled_vecs : n_med x nb_of_vecs	( RC ) -> nb_of_vecs vectors of dim n_coupl

	// Set the first nullspace vectors
	std::string input_filename = input_filename_base + "_0_n_" + std::to_string(m_null_nb_vecs) + ".petscvec";
	VecCreate(m_comm.get(),&m_null_vecs[0]);
	VecSetSizes(m_null_vecs[0],m_null_vecs_N_local,m_null_vecs_N);
	read_PETSC_vector(m_null_vecs[0],input_filename, m_comm.get());
	
	std::string output_filename = m_scratch_folder_path + "/rb_coupl_vector_0_n_" + std::to_string(m_null_nb_vecs) + ".petscvec";
	VecCreate(m_comm.get(),&m_null_coupled_vecs[0]);
	VecSetSizes(m_null_coupled_vecs[0],m_C_RR_M_local,m_C_RR_M);
	VecSetFromOptions(m_null_coupled_vecs[0]);
	MatMult(*coupling_matrix,m_null_vecs[0],m_null_coupled_vecs[0]);
	write_PETSC_vector(m_null_coupled_vecs[0],output_filename,m_comm.rank(),m_comm.get());

	// DEBUG print
	// write_PETSC_vector_MATLAB(m_null_coupled_vecs[0],m_scratch_folder_path + "/rb_coupl_vector_0_n_" + std::to_string(m_null_nb_vecs) + ".m",m_comm.get());
	// write_PETSC_vector_MATLAB(m_null_vecs[0],m_scratch_folder_path + "/rb_vector_0_n_" + std::to_string(m_null_nb_vecs) + ".m",m_comm.get());

	// Read and calculate the rest of the nullspace vectors
	for(int iii = 1; iii < m_null_nb_vecs; ++iii)
	{
		input_filename = input_filename_base + "_" + std::to_string(iii) + "_n_" + std::to_string(m_null_nb_vecs) + ".petscvec";
		VecDuplicate(m_null_vecs[0],&m_null_vecs[iii]);
		read_PETSC_vector(m_null_vecs[iii],input_filename, m_comm.get());

		std::string output_filename = m_scratch_folder_path + "/rb_coupl_vector_" + std::to_string(iii) + "_n_" + std::to_string(m_null_nb_vecs) + ".petscvec";
		VecDuplicate(m_null_coupled_vecs[0],&m_null_coupled_vecs[iii]);
		MatMult(*coupling_matrix,m_null_vecs[iii],m_null_coupled_vecs[iii]);

		write_PETSC_vector(m_null_coupled_vecs[iii],output_filename,m_comm.rank(),m_comm.get());

		// DEBUG print
		// write_PETSC_vector_MATLAB(m_null_coupled_vecs[iii],m_scratch_folder_path + "/rb_coupl_vector_" + std::to_string(iii) + "_n_" + std::to_string(m_null_nb_vecs) + ".m",m_comm.get());
		// write_PETSC_vector_MATLAB(m_null_vecs[iii],m_scratch_folder_path + "/rb_vector_" + std::to_string(iii) + "_n_" + std::to_string(m_null_nb_vecs) + ".m",m_comm.get());
	}

	// Build the LOCAL dense matrix
	std::vector<PetscScalar> dummy_vec_val(m_null_nb_vecs,0);
	std::vector<PetscInt>    dummy_vec_row(m_null_nb_vecs,0);

	for(PetscInt iii = 0; iii < m_null_nb_vecs; ++iii)
	{
		dummy_vec_row[iii] = iii;
	}

	MatCreateSeqDense(PETSC_COMM_SELF,m_null_nb_vecs,m_null_nb_vecs,NULL,&m_RITRI_mat);

	for(PetscInt iii = 0; iii < m_null_nb_vecs; ++iii)
	{
		VecMDot(m_null_coupled_vecs[iii],m_null_nb_vecs,m_null_coupled_vecs,dummy_vec_val.data());
		MatSetValues(m_RITRI_mat,m_null_nb_vecs,dummy_vec_row.data(),1,&iii,dummy_vec_val.data(),INSERT_VALUES);
	}
	MatAssemblyBegin(m_RITRI_mat,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(m_RITRI_mat,MAT_FINAL_ASSEMBLY);

	PETSC_invert_dense_matrix(m_RITRI_mat,m_inv_RITRI_mat);

	if(m_comm.rank() == 0)
	{
		write_PETSC_matrix(m_inv_RITRI_mat,m_scratch_folder_path + "/rb_inv_RITRI.petscmat",0,PETSC_COMM_SELF);

		// DEBUG print
		// write_PETSC_matrix_MATLAB(m_inv_RITRI_mat,m_scratch_folder_path + "/rb_inv_RITRI.m",PETSC_COMM_SELF);
	}

	// Set up flag
	m_bNullVecsSet = true;
	m_binvRITRIMatSet = true;

	// Cleanup
	MatDestroy(&m_RITRI_mat);
}

void FETI_Operations::read_null_space_vecs(const std::string& RB_vectors_base, int nb_of_rb_vectors)
{

	/*
	 *		Note: for now, this function reads both sets of null space vectors (the "original" ones and
	 *  the ones multiplied by the coupling matrix), but I'm not sure if this less efficient than re-calculating
	 *  these vectors each time. [Thiago]
	 */

	homemade_assert_msg(m_bNullVecsDimensionsSet,"Null vectors sizes not set yet!");
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

	// Set the null vec arrays
	m_null_nb_vecs = nb_of_rb_vectors;

	// Set the first nullspace vectors
	std::string rb_vec_filename = RB_vectors_base + "_0_n_" + std::to_string(m_null_nb_vecs) + ".petscvec";
	VecCreate(m_comm.get(),&m_null_vecs[0]);
	VecSetSizes(m_null_vecs[0],m_null_vecs_N_local,m_null_vecs_N);
	read_PETSC_vector(m_null_vecs[0],rb_vec_filename, m_comm.get());
	
	std::string rb_coupl_vec_filename = m_scratch_folder_path + "/rb_coupl_vector_0_n_" + std::to_string(m_null_nb_vecs) + ".petscvec";
	VecCreate(m_comm.get(),&m_null_coupled_vecs[0]);
	VecSetSizes(m_null_coupled_vecs[0],m_C_RR_M_local,m_C_RR_M);
	read_PETSC_vector(m_null_coupled_vecs[0],rb_coupl_vec_filename, m_comm.get());

	// Read and calculate the rest of the nullspace vectors
	for(int iii = 1; iii < m_null_nb_vecs; ++iii)
	{
		rb_vec_filename = RB_vectors_base + "_" + std::to_string(iii) + "_n_" + std::to_string(m_null_nb_vecs) + ".petscvec";
		VecDuplicate(m_null_vecs[0],&m_null_vecs[iii]);
		read_PETSC_vector(m_null_vecs[iii],rb_vec_filename, m_comm.get());

		rb_coupl_vec_filename = m_scratch_folder_path + "/rb_coupl_vector_" + std::to_string(iii) + "_n_" + std::to_string(m_null_nb_vecs) + ".petscvec";
		VecDuplicate(m_null_coupled_vecs[0],&m_null_coupled_vecs[iii]);
		read_PETSC_vector(m_null_coupled_vecs[iii],rb_coupl_vec_filename, m_comm.get());
	}

	// Set up flag
	m_bNullVecsSet = true;
}

void FETI_Operations::read_null_space_inv_RITRI_mat()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	
	// A dummy vector which will be used to sync the matrices
	std::vector<PetscScalar> dummy_matrix_contents(m_null_nb_vecs*m_null_nb_vecs, 0);

	MatCreateSeqDense(PETSC_COMM_SELF,m_null_nb_vecs,m_null_nb_vecs,NULL,&m_inv_RITRI_mat);

	// Only really read the matrix in the first processor
	if(m_comm.rank() == 0)
	{
		read_PETSC_matrix(m_inv_RITRI_mat,m_scratch_folder_path + "/rb_inv_RITRI.petscmat",PETSC_COMM_SELF);

		// Get the data
		PetscScalar *dummy_array;
		MatDenseGetArray(m_inv_RITRI_mat,&dummy_array);

		for(int iii = 0; iii < m_null_nb_vecs*m_null_nb_vecs; ++iii)
		{
			dummy_matrix_contents[iii] = dummy_array[iii];
		}

		MatDenseRestoreArray(m_inv_RITRI_mat,&dummy_array);
	}

	// Sync the matrix!
	m_comm.barrier();
	m_comm.broadcast(dummy_matrix_contents);

	if(m_comm.rank() != 0)
	{
		std::vector<PetscInt> dummy_indexes(m_null_nb_vecs,0);

		for(int iii = 0; iii < m_null_nb_vecs; ++iii)
		{
			dummy_indexes[iii] = iii;
		}

		MatSetValues(m_inv_RITRI_mat,m_null_nb_vecs,&dummy_indexes[0],m_null_nb_vecs,&dummy_indexes[0],&dummy_matrix_contents[0],INSERT_VALUES);
		MatAssemblyBegin(m_inv_RITRI_mat,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(m_inv_RITRI_mat,MAT_FINAL_ASSEMBLY);
	}

	// Set up flag
	m_binvRITRIMatSet = true;
}

void FETI_Operations::calculate_null_space_phi_0(const std::string& force_path)
{
	// 	phi0 = - RC * (inv_RITRI_mat) * R^t* Force
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bNullVecsSet,"Null space vectors not set yet!");
	homemade_assert_msg(m_binvRITRIMatSet,"Null space matrices not set yet!");
	homemade_assert_msg(m_bC_R_micro_MatrixSet,"Micro coupling matrix not set yet!");
	homemade_assert_msg(m_bC_R_BIG_MatrixSet,"Macro coupling matrix not set yet!");

	// --- Declare / create the vectors
	// phi(0)
	VecCreate(m_comm.get(),&m_current_phi);
	VecSetSizes(m_current_phi,m_C_RR_M_local,m_C_RR_M);
	VecSetFromOptions(m_current_phi);

	// Force vector
	Vec vec_force_PETSc;
	VecCreate(m_comm.get(),&vec_force_PETSc);
	read_PETSC_vector(vec_force_PETSc,force_path,m_comm.get());

	// Create the local auxiliary vectors
	Vec aux_null_vec_input;
	Vec aux_null_vec_output;

	VecCreateSeq(PETSC_COMM_SELF,m_null_nb_vecs,&aux_null_vec_input);
	VecCreateSeq(PETSC_COMM_SELF,m_null_nb_vecs,&aux_null_vec_output);

	VecZeroEntries(aux_null_vec_input);
	VecZeroEntries(aux_null_vec_output);

	// aux_vec_input = R^t * vec_force
	// -> All the communications are done here!
	// Cannot use VecMDot because "m_null_vecs" is a constant pointer ...
	PetscScalar *dummy_array_input;
	VecGetArray(aux_null_vec_input,&dummy_array_input);
	for(int iii = 0; iii < m_null_nb_vecs; ++iii)
	{
		VecDot(vec_force_PETSc,m_null_vecs[iii],&dummy_array_input[iii]);
	}
	VecRestoreArray(aux_null_vec_input,&dummy_array_input);

	// aux_vec_output = inv_RITRI_mat * aux_vec_input
	// -> Completely local operation!
	MatMult(m_inv_RITRI_mat,aux_null_vec_input,aux_null_vec_output);

	// vec_out = sum ( aux_null_vec_output[i] * vec_RC[i])
	// -> This should have no communications at all!
	VecZeroEntries(m_current_phi);

	PetscScalar *dummy_array_output;
	VecGetArray(aux_null_vec_output,&dummy_array_output);
	VecMAXPY(m_current_phi,m_null_nb_vecs,dummy_array_output,m_null_coupled_vecs);
	VecRestoreArray(aux_null_vec_output,&dummy_array_output);

	VecScale(m_current_phi,-1);

	// Cleanup
	VecDestroy(&vec_force_PETSc);
	VecDestroy(&aux_null_vec_input);
	VecDestroy(&aux_null_vec_output);

	// Set flags
	m_bSet_current_phi = true;
}

//  --- FETI read methods
void FETI_Operations::read_decoupled_solutions()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bC_R_micro_MatrixSet,"Micro system dimensions not set yet!");
	homemade_assert_msg(m_bC_R_BIG_MatrixSet,"Macro system dimensions not set yet!");

	// Create the vectors
	VecCreate(m_comm.get(),&m_u_0_BIG);
	VecSetSizes(m_u_0_BIG,m_C_R_BIG_N_local,m_C_R_BIG_N);
	read_PETSC_vector(m_u_0_BIG,m_scratch_folder_path + "/ext_solver_u0_A_sys_sol_vec.petscvec", m_comm.get());

	VecCreate(m_comm.get(),&m_u_0_micro);
	VecSetSizes(m_u_0_micro,m_C_R_micro_N_local,m_C_R_micro_N);
	read_PETSC_vector(m_u_0_micro,m_scratch_folder_path + "/ext_solver_u0_B_sys_sol_vec.petscvec", m_comm.get());

	// Set the flag
	m_bSet_u_0 = true;
}

void FETI_Operations::read_ext_solver_output()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bC_R_micro_MatrixSet,"Micro system dimensions not set yet!");
	homemade_assert_msg(m_bC_R_BIG_MatrixSet,"Macro system dimensions not set yet!");

	// Create the vectors
	VecCreate(m_comm.get(),&m_ext_solver_sol_BIG);
	VecSetSizes(m_ext_solver_sol_BIG,m_C_R_BIG_N_local,m_C_R_BIG_N);
	read_PETSC_vector(m_ext_solver_sol_BIG,m_scratch_folder_path + "/ext_solver_A_sys_sol_vec.petscvec", m_comm.get());

	VecCreate(m_comm.get(),&m_ext_solver_sol_micro);
	VecSetSizes(m_ext_solver_sol_micro,m_C_R_micro_N_local,m_C_R_micro_N);
	read_PETSC_vector(m_ext_solver_sol_micro,m_scratch_folder_path + "/ext_solver_B_sys_sol_vec.petscvec", m_comm.get());

	// Set the flag
	m_bSet_ext_solver_sol = true;
}

void FETI_Operations::read_previous_phi()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

	// Create the vectors
	VecCreate(m_comm.get(),&m_previous_phi);
	VecSetSizes(m_previous_phi,m_C_RR_M_local,m_C_RR_M);
	read_PETSC_vector(m_previous_phi,m_scratch_folder_path + "/FETI_iter__phi__current.petscvec", m_comm.get());

	// Set flag
	m_bSet_previous_phi = true;
}

void FETI_Operations::read_previous_r()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

	// Create the vectors
	VecCreate(m_comm.get(),&m_previous_residual);
	VecSetSizes(m_previous_residual,m_C_RR_M_local,m_C_RR_M);
	read_PETSC_vector(m_previous_residual,m_scratch_folder_path + "/FETI_iter__r__current.petscvec", m_comm.get());

	// Set flag
	m_bSet_previous_residual = true;
}

void FETI_Operations::read_all_previous_p()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSet_previous_residual,"Previous residual vector not set yet!");

	// Create the vectors
	m_previous_p_ptr = new Vec[m_kkk+1];
	VecDuplicateVecs(m_previous_residual,m_kkk+1,&m_previous_p_ptr);

	for(int iii = 0; iii < m_kkk+1; ++iii)
	{
		read_PETSC_vector(m_previous_p_ptr[iii],m_scratch_folder_path + "/FETI_iter__p__" + std::to_string(iii) + ".petscvec", m_comm.get());
	}

	// Set flag
	m_bSet_previous_p_ptr = true;
}

void FETI_Operations::read_all_previous_q()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSet_previous_residual,"Previous residual vector not set yet!");

	// Create q(0 ... kkk), read q(0 ... kkk - 1)
	m_previous_q_ptr = new Vec[m_kkk+1];
	VecDuplicateVecs(m_previous_residual,m_kkk+1,&m_previous_q_ptr);
	if(m_kkk > 0)
	{
		for(int iii = 0; iii < m_kkk; ++iii)
		{
			read_PETSC_vector(m_previous_q_ptr[iii],m_scratch_folder_path + "/FETI_iter__q__" + std::to_string(iii) + ".petscvec", m_comm.get());
		}
	}

	// Set flags
	m_bSet_previous_q_ptr = true;
}

void FETI_Operations::read_scalar_data()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

	// Import the scalar data
	// File format:  kkk   rho(0)   rho(kkk)   [ RB_corr(kkk) ]

	// ONLY read in proc 0!
	if(m_comm.rank() == 0)
	{
		std::ifstream scalar_data;

		scalar_data.open(m_scratch_folder_path + "/FETI_iter_scalar_data.dat");
		scalar_data >> m_kkk >> m_rho_0 >> m_previous_rho;

		if(m_bUsingNullVecs)
		{
			scalar_data >> m_previous_RB_mode_corr;
		}
		scalar_data.close();
	}

	// Broadcast the values
	m_comm.broadcast(m_kkk);
	m_comm.broadcast(m_rho_0);
	m_comm.broadcast(m_previous_rho);
	if(m_bUsingNullVecs)
	{
		m_comm.broadcast(m_previous_RB_mode_corr);
	}

	// Import `p(0 ... kkk - 1).q(0 ... kkk - 1)`
	m_p_dot_q.resize(m_kkk+1,0);
	if(m_kkk != 0)
	{
		// ONLY read in proc 0!
		if(m_comm.rank() == 0)
		{
			std::ifstream scalar_data;
			scalar_data.open(m_scratch_folder_path + "/FETI_iter_p_dot_q.dat");

			for(int iii = 0; iii < m_kkk; ++iii)
			{
				scalar_data >> m_p_dot_q[iii];
			}
			scalar_data.close();
		}

		// Broadcast the values
		m_comm.broadcast(m_p_dot_q);
	}

	// Set flag
	m_bReadPreviousScalar = true;
}

void FETI_Operations::read_vector_data()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bReadPreviousScalar,"Scalar data not set yet!");

	/* Vectors to read: 'r(kkk)'
	 *					'phi(kkk)'
	 *					'p(0 ... kkk)'
	 * 					'q(0 ... kkk - 1)'
	 */

	// Create and read the vectors
	// r(kkk)
	this->read_previous_r();

	// phi(kkk)
	this->read_previous_phi();

	// p(0 ... kkk)
	this->read_all_previous_p();

	// q(0 ... kkk)
	// Create q(0 ... kkk), read q(0 ... kkk - 1)
	this->read_all_previous_q();
}

//  --- FETI steps methods
void FETI_Operations::calculate_initial_p()
{
	// `p(0)` is identical to `z(0)`, so just calculate the latter
	this->calculate_z();
}

void FETI_Operations::calculate_initial_r()
{
	/*
	 *		If m_bUsingNullVecs == true :
	 *			r(0) = ( C_1 * u_0,1 - C_2 * u_0,2 ) - C_1 * x_0,1 - C_2 * x_0,2
	 *
	 *		eles
	 *			r(0) = ( C_1 * u_0,1 - C_2 * u_0,2 )
	 */

	// Common asserts
	homemade_assert_msg(m_bC_R_micro_MatrixSet,"Micro coupling matrix not set yet!");
	homemade_assert_msg(m_bC_R_BIG_MatrixSet,"Macro coupling matrix not set yet!");
	homemade_assert_msg(m_bC_R_micro_MatrixSet,"Micro system dimensions not set yet!");
	homemade_assert_msg(m_bC_R_BIG_MatrixSet,"Macro system dimensions not set yet!");
	homemade_assert_msg(m_bSet_u_0,"Decoupled solution not set yet!");

	// Assert for only m_bUsingNullVecs == true
	homemade_assert_msg(m_bSet_ext_solver_sol && m_bUsingNullVecs,"Ext. solver solutions not set yet!");

	// Create the vector
	VecCreate(m_comm.get(),&m_current_residual);
	VecSetSizes(m_current_residual,m_C_RR_M_local,m_C_RR_M);
	VecSetFromOptions(m_current_residual);

	// - C_2 * u_0,2
	MatMult(m_C_R_micro,m_u_0_micro,m_current_residual);
	VecScale(m_current_residual,-1);

	// C_1 * u_0,1
	// --> r(0) = ( C_1 * u_0,1 - C_2 * u_0,2 )
	MatMultAdd(m_C_R_BIG,m_u_0_BIG,m_current_residual,m_current_residual);
	
	if(m_bUsingNullVecs)
	{
		Vec dummy_vec;
		VecDuplicate(m_current_residual,&dummy_vec);

		// C_1 * x_0,1
		MatMult(m_C_R_BIG,m_ext_solver_sol_BIG,dummy_vec);

		// C_2 * x_0,2
		MatMultAdd(m_C_R_micro,m_ext_solver_sol_micro,dummy_vec,dummy_vec);

		// --> r(0) = ( C_1 * u_0,1 - C_2 * u_0,2 ) - C_1 * x_0,1 - C_2 * x_0,2
		VecAXPY(m_current_residual,-1,dummy_vec);
		VecDestroy(&dummy_vec);
	}

	// Set flag
	m_bSet_current_residual = true;
}

void FETI_Operations::calculate_p()
{
	homemade_assert_msg(m_bSet_current_z,"Current 'z' not calculated yet!");
	homemade_assert_msg(m_bSet_previous_p_ptr,"'p' vectors not set yet!");
	homemade_assert_msg(m_bSet_previous_q_ptr,"'q' vectors not set yet!");

	// p(kkk+1) = z(kkk+1) + \sum ( beta(iii) * p(iii) ), iii = 0 ... kkk
	// beta(iii) = - z(kkk+1).q(iii) / p(iii).q(iii)
	// For simplicity, the signal was included in 'beta' - WATCH OUT FOR THE SIGNAL!
	std::vector<PetscScalar> beta(m_kkk+1,0);
	VecMDot(m_current_z,m_kkk+1,m_previous_q_ptr,beta.data());
	for(int iii = 0; iii < m_kkk+1; ++iii)
	{
		beta[iii] = - beta[iii] / m_p_dot_q[iii];
	}

	// p(kkk+1) = z(kkk+1)
	VecDuplicate(m_current_z,&m_current_p);
	VecCopy(m_current_z,m_current_p);

	// p(kkk+1) = z(kkk+1) + \sum ( beta(iii) * p(iii) )
	VecMAXPY(m_current_p,m_kkk+1,beta.data(),m_previous_q_ptr);

	m_bSet_current_p = true;
}

void FETI_Operations::calculate_q()
{
	homemade_assert_msg(m_bC_R_micro_MatrixSet,"Micro coupling matrix not set yet!");
	homemade_assert_msg(m_bC_R_BIG_MatrixSet,"Macro coupling matrix not set yet!");
	homemade_assert_msg(m_bSet_ext_solver_sol,"Ext. solver solutions not set yet!");
	homemade_assert_msg(m_bSet_previous_q_ptr,"'q' vectors not set yet!");

	// q(kkk) = C_1 * x_1(kkk) + C_2 * x_2(kkk)
	MatMult(m_C_R_BIG,m_ext_solver_sol_BIG,m_previous_q_ptr[m_kkk]);
	MatMultAdd(m_C_R_micro,m_ext_solver_sol_micro,m_previous_q_ptr[m_kkk],m_previous_q_ptr[m_kkk]);
}

void FETI_Operations::calculate_r()
{
	homemade_assert_msg(m_bSet_previous_q_ptr,"'q' vectors not set yet!");
	homemade_assert_msg(m_bSet_previous_residual,"Previous 'r' not set yet!");
	homemade_assert_msg(m_bSet_current_phi,"Current 'phi' not set yet!");
	// Last assert is needed to guarantee that 'm_p_dot_q(kkk)' has been calculated

	// Calculate 'r(kkk + 1) = r(kkk) - gamma * q(kkk)'

	double dummy_double = - m_gamma;
	VecDuplicate(m_previous_residual,&m_current_residual);
	VecWAXPY(m_current_residual,dummy_double,m_previous_q_ptr[m_kkk],m_previous_residual);

	m_bSet_current_residual = true;
}

void FETI_Operations::calculate_z()
{
	/*	
	 *  Four possibilities:
	 *	 - m_bUsingNullVecs == false
	 *     m_precond_type == BaseCGPrecondType::NO_PRECONDITIONER
	 *     -> m_current_z = m_current_residual
	 *
	 *	 - m_bUsingNullVecs == false
	 *     m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER
	 *     -> m_current_z = M_PC^-1 * m_current_residual
	 *
	 *	 - m_bUsingNullVecs == true
	 *     m_precond_type == BaseCGPrecondType::NO_PRECONDITIONER
	 *     -> m_current_z = M_proj * m_current_residual
	 *
	 *	 - m_bUsingNullVecs == true
	 *     m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER
	 *     -> m_current_z = M_proj * M_PC^-1 * M_proj * m_current_residual
	 */

	homemade_assert_msg((m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER) && m_bC_RR_MatrixSet,"Preconditioner matrix not set yet!");
	homemade_assert_msg(m_bSet_current_residual,"Current residual not calculated yet!");

	if(m_bUsingNullVecs)
	{
		// Use projections!
		if(m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER)
		{
	 		// -> m_current_z = M_proj * M_PC^-1 * M_proj * m_current_residual
			VecDuplicate(m_current_residual,&m_current_z);

			Vec dummy_vec, dummy_vec_bis;
			VecDuplicate(m_current_residual,&dummy_vec);
			VecDuplicate(m_current_residual,&dummy_vec_bis);

			// dummy_vec = M_proj * m_current_residual
			this->apply_RB_projection(m_current_residual,dummy_vec);

			// dummy_vec_bis =  M_PC^-1 * dummy_vec
			this->apply_precond(dummy_vec,dummy_vec_bis);

			// m_current_z = M_proj * dummy_vec_bis
			this->apply_RB_projection(dummy_vec_bis,m_current_z);

			VecDestroy(&dummy_vec);
			VecDestroy(&dummy_vec_bis);
			// Set flag
			m_bSet_current_z = true;
		}
		else
		{
			// -> m_current_z = M_proj * m_current_residual
			VecDuplicate(m_current_residual,&m_current_z);

			this->apply_RB_projection(m_current_residual,m_current_z);
			m_bSet_current_z = true;
		}
	}
	else
	{
		// Do not use projections!
		if(m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER)
		{
			// -> m_current_z = M_PC^-1 * m_current_residual
			VecDuplicate(m_current_residual,&m_current_z);

			this->apply_precond(m_current_residual,m_current_z);
			m_bSet_current_z = true;
		}
		// else, m_current_z = m_current_residual -> DO NOTHING!
	}
}

void FETI_Operations::calculate_phi()
{
	homemade_assert_msg(m_bSet_previous_p_ptr,"'p' vectors not set yet!");
	homemade_assert_msg(m_bSet_previous_q_ptr,"'q' vectors not set yet!");
	homemade_assert_msg(m_bSet_previous_phi,"Previous 'phi' not set yet!");

	// Calculate 'phi(kkk + 1) = phi(kkk) + gamma * p(kkk)'

	// gamma = rho(kkk) / ( p(kkk).q(kkk) )
	VecDot(m_previous_p_ptr[m_kkk],m_previous_q_ptr[m_kkk],&m_p_dot_q[m_kkk]);
	m_gamma = m_previous_rho / m_p_dot_q[m_kkk];

	// phi(kkk + 1) = phi(kkk) + gamma * p(kkk)
	VecDuplicate(m_previous_phi,&m_current_phi);
	VecWAXPY(m_current_phi,m_gamma,m_previous_p_ptr[m_kkk],m_previous_phi);

	m_bSet_current_phi = true;
}

void FETI_Operations::calculate_rb_correction()
{
	homemade_assert_msg(m_bNullVecsSet,"Null space vectors not set yet!");
	homemade_assert_msg(m_bSet_current_residual,"Current residual not calculated yet!");

	// Declare and create vectors
	Vec dummy_seq_vec;
	Vec dummy_seq_vec_bis;
	VecCreateSeq(PETSC_COMM_SELF,m_null_nb_vecs,&dummy_seq_vec);
	VecZeroEntries(dummy_seq_vec);
	VecDuplicate(dummy_seq_vec,&dummy_seq_vec_bis);

	VecDuplicate(m_null_vecs[0],&m_current_rb_correction);
	VecZeroEntries(m_current_rb_correction);

	// m_current_rb_correction = R * (inv_RITRI_mat) * RC^t * m_current_residual

	// dummy_seq_vec = RC^t * m_current_residual
	// -> All the communications are done here!
	PetscScalar *dummy_seq_array;
	VecGetArray(dummy_seq_vec,&dummy_seq_array);
	VecMDot(m_current_residual,m_null_nb_vecs,m_null_coupled_vecs,dummy_seq_array);
	VecRestoreArray(dummy_seq_vec,&dummy_seq_array);

	// dummy_seq_vec_bis = inv_RITRI_mat * dummy_seq_vec
	// -> Completely local operation!
	MatMult(m_inv_RITRI_mat,dummy_seq_vec,dummy_seq_vec_bis);

	// m_current_rb_correction = sum ( dummy_seq_vec_bis[i] * m_null_vecs[i])
	// -> This should have no communications at all!
	VecGetArray(dummy_seq_vec_bis,&dummy_seq_array);
	VecMAXPY(m_current_rb_correction,m_null_nb_vecs,dummy_seq_array,m_null_vecs);
	VecRestoreArray(dummy_seq_vec_bis,&dummy_seq_array);

	// Cleanup
	VecDestroy(&dummy_seq_vec);
	VecDestroy(&dummy_seq_vec_bis);

	// Set flag
	m_bSet_current_RB_correction = true;
}

void FETI_Operations::calculate_scalar_data()
{
	homemade_assert_msg(m_bSet_current_residual,"Current residual not calculated yet!");
	homemade_assert_msg(m_bSet_current_z,"Current 'z' not calculated yet!");
	if(m_bUsingNullVecs)
	{
		homemade_assert_msg(m_bSet_current_RB_correction,"Current RB modes correction not calculated yet!");
	}

	// --- Calculate the values
	// rho(kkk+1)
	if(m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER || m_bUsingNullVecs)
	{
		VecDot(m_current_residual,m_current_z,&m_current_rho);
	} else {
		VecDot(m_current_residual,m_current_residual,&m_current_rho);
	}

	// RB_corr(kkk+1)
	if(m_bUsingNullVecs)
	{
		VecDot(m_current_rb_correction,m_current_rb_correction,&m_current_RB_mode_corr);
	}

	// gamma(kkk)
	VecDot(m_previous_p_ptr[m_kkk],m_previous_q_ptr[m_kkk],&m_p_dot_q[m_kkk]);

	// Set flags
	m_bCalculatedScalar = true;
}

IterationStatus FETI_Operations::check_convergence(double rel_residual_conv, double abs_residual_conv, int max_iter_div, double rel_residual_div, double rb_modes_conv)
{
	// NOT IMPLEMENTED YET!!!
	return IterationStatus::ITERATING;
}

//  --- Write methods
void FETI_Operations::export_ext_solver_rhs_iteration()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

	if(m_kkk == 0)
	{
		// p(0) = z(0)
		homemade_assert_msg(m_bSet_current_z,"Current 'p' not calculated yet!");
		this->export_ext_solver_rhs(m_current_z);
	} else {
		homemade_assert_msg(m_bSet_current_p,"Current 'p' not calculated yet!");
		this->export_ext_solver_rhs(m_current_p);
	}
}

void FETI_Operations::export_ext_solver_rhs_decoupled()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSet_current_phi,"Current 'phi' not calculated yet!");

	this->export_ext_solver_rhs(m_current_phi);
}

void FETI_Operations::export_p()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSet_current_p,"Current 'p' not calculated yet!");

	write_PETSC_vector(m_current_p,m_scratch_folder_path + "/FETI_iter__p__" + std::to_string(m_kkk+1) + ".petscvec",m_comm.rank(),m_comm.get());
	write_PETSC_vector_MATLAB(m_current_p,m_scratch_folder_path + "/FETI_iter__p__" + std::to_string(m_kkk+1) + ".m",m_comm.get());
}

void FETI_Operations::export_q()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSet_previous_q_ptr,"'q' vectors not calculated yet!");

	write_PETSC_vector(m_previous_q_ptr[m_kkk],m_scratch_folder_path + "/FETI_iter__q__" + std::to_string(m_kkk) + ".petscvec",m_comm.rank(),m_comm.get());
	write_PETSC_vector_MATLAB(m_previous_q_ptr[m_kkk],m_scratch_folder_path + "/FETI_iter__q__" + std::to_string(m_kkk) + ".m",m_comm.get());
}

void FETI_Operations::export_r()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSet_current_residual,"Current residual not calculated yet!");

	write_PETSC_vector(m_current_residual,m_scratch_folder_path + "/FETI_iter__r__current.petscvec",m_comm.rank(),m_comm.get());
	write_PETSC_vector_MATLAB(m_current_residual,m_scratch_folder_path + "/FETI_iter__r__current.m",m_comm.get());
}

void FETI_Operations::export_phi()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSet_current_phi,"Current 'phi' not calculated yet!");

	write_PETSC_vector(m_current_phi,m_scratch_folder_path + "/FETI_iter__phi__current.petscvec",m_comm.rank(),m_comm.get());
	write_PETSC_vector_MATLAB(m_current_phi,m_scratch_folder_path + "/FETI_iter__phi__current.m",m_comm.get());
}

void FETI_Operations::export_inital_vecs()
{
	// Export r(0) and p(0) ( p(0) is identical to z(0))
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bSet_current_residual,"Current residual not calculated yet!");

	// In all cases, print r(0)
	write_PETSC_vector(m_current_residual,m_scratch_folder_path + "/FETI_iter__r__current.petscvec",m_comm.rank(),m_comm.get());
	write_PETSC_vector_MATLAB(m_current_residual,m_scratch_folder_path + "/FETI_iter__r__current.m",m_comm.get());

	// p(0) is only identical to r(0) if neither a preconditioner or the RB modes are used
	if(m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER || m_bUsingNullVecs)
	{
		write_PETSC_vector(m_current_z,m_scratch_folder_path + "/FETI_iter__p__0.petscvec",m_comm.rank(),m_comm.get());
		write_PETSC_vector_MATLAB(m_current_z,m_scratch_folder_path + "/FETI_iter__p__0.m",m_comm.get());
	}
}

void FETI_Operations::export_initial_scalar_data()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");


	homemade_assert_msg(m_bSet_current_residual,"Current residual not calculated yet!");
	PetscScalar residual = 0;
	PetscScalar RB_correct = 0;

	if(m_bUsingNullVecs)
	{
		homemade_assert_msg(m_bSet_current_RB_correction,"Current RB modes correction not calculated yet!");
	}

	// --- Calculate the values
	// rho(0)
	if(m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER || m_bUsingNullVecs)
	{
		VecDot(m_current_residual,m_current_z,&residual);
	} else {
		VecDot(m_current_residual,m_current_residual,&residual);
	}

	// RB_corr(0)
	if(m_bUsingNullVecs)
	{
		VecDot(m_current_rb_correction,m_current_rb_correction,&RB_correct);
	}

	// ONLY write in proc 0!
	if(m_comm.rank() == 0)
	{
		std::ofstream scalar_data;

		// Export the scalar data
		scalar_data.open(m_scratch_folder_path + "/FETI_iter_scalar_data.dat");
		scalar_data << m_kkk << " " << residual << " " << residual;

		if(m_bUsingNullVecs)
		{
			scalar_data <<  " " << RB_correct;
		}

		scalar_data << std::endl;

		scalar_data.close();

		// Export the convergence data
		scalar_data.open(m_scratch_folder_path + "/FETI_convergence.dat");
		scalar_data << m_kkk << " " << residual;

		if(m_bUsingNullVecs)
		{
			scalar_data <<  " " << RB_correct;
		}

		scalar_data << std::endl;

		scalar_data.close();
	}
}

void FETI_Operations::export_scalar_data()
{
	homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
	homemade_assert_msg(m_bCalculatedScalar,"Scalar data not calculated yet!");

	// ONLY write in proc 0!
	if(m_comm.rank() == 0)
	{
		std::ofstream scalar_data;

		// Export the scalar data
		scalar_data.open(m_scratch_folder_path + "/FETI_iter_scalar_data.dat");
		scalar_data << m_kkk << " " << m_rho_0 << " " << m_current_rho;

		if(m_bUsingNullVecs)
		{
			scalar_data <<  " " << m_current_RB_mode_corr;
		}

		scalar_data << std::endl;

		scalar_data.close();

		// Export the convergence data
		scalar_data.open(m_scratch_folder_path + "/FETI_convergence.dat",std::ofstream::app);
		scalar_data << m_kkk + 1 << " " << m_current_rho;

		if(m_bUsingNullVecs)
		{
			scalar_data <<  " " << m_current_RB_mode_corr;
		}
		scalar_data << std::endl;

		scalar_data.close();

		// Export `p(kkk).q(kkk)`
		if(m_kkk == 0)
		{
			scalar_data.open(m_scratch_folder_path + "/FETI_iter_p_dot_q.dat");
		} else {

			scalar_data.open(m_scratch_folder_path + "/FETI_iter_p_dot_q.dat",std::ofstream::app);
		}

		scalar_data << m_p_dot_q[m_kkk] << std::endl;

		scalar_data.close();
	}
}

void FETI_Operations::export_iter_vecs()
{
	/* Export the iteration vectors
	 * Vectors to export: 'r(kkk+1)'
	 *                    'phi(kkk+1)'
	 *                    'p(kkk+1)'
	 *                    'q(kkk)'
	 */

	// r(kkk+1)
	this->export_r();

	// phi(kkk+1)
	this->export_phi();

	// p(kkk+1)
	this->export_p();

	// q(kkk)
	this->export_q();
}

void FETI_Operations::print_previous_iters_conv(int nb_of_iters)
{
	// NOT IMPLEMENTED YET!!!
}

}