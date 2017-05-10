#include "FETI_operations.h"

namespace carl
{
//  --- Null space / rigid body modes methods
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
	 *			r(0) = ( C_1 * u_1,0 - C_2 * u_2,0 ) - C_1 * x_1(0) - C_2 * x_2(0)
	 *
	 *		eles
	 *			r(0) = ( C_1 * u_1,0 - C_2 * u_2,0 )
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

		// C_1 * x_1(0)
		MatMult(m_C_R_BIG,m_ext_solver_sol_BIG,dummy_vec);

		// C_2 * x_2(0)
		MatMultAdd(m_C_R_micro,m_ext_solver_sol_micro,dummy_vec,dummy_vec);

		// --> r(0) = ( C_1 * u_1,0 - C_2 * u_2,0 ) - C_1 * x_1(0) - C_2 * x_2(0)
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
	VecMAXPY(m_current_p,m_kkk+1,beta.data(),m_previous_p_ptr);

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
	homemade_assert_msg(m_bSet_current_RB_correction && m_bUsingNullVecs,"Current RB modes correction not calculated yet!");

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
		VecNorm(m_current_rb_correction,NORM_2,&m_current_RB_mode_corr);
	}

	// gamma(kkk)
	VecDot(m_previous_p_ptr[m_kkk],m_previous_q_ptr[m_kkk],&m_p_dot_q[m_kkk]);

	// Set flags
	m_bCalculatedScalar = true;
}

void FETI_Operations::calculate_coupled_solution()
{
	homemade_assert_msg(m_bSet_ext_solver_sol,"Ext. solver solutions not set yet!");
	homemade_assert_msg(m_bSet_u_0,"Decoupled solution not set yet!");
	homemade_assert_msg(m_bSet_current_RB_correction && m_bUsingNullVecs,"RB modes correction not set yet!");

	// Set the solution vectors
	VecDuplicate(m_u_0_micro,&m_coupled_sol_micro);
	VecDuplicate(m_u_0_BIG,&m_coupled_sol_BIG);

	// u_1 = - x_1(FINAL) + u_0,1 
	VecWAXPY(m_coupled_sol_BIG,-1,m_ext_solver_sol_BIG,m_u_0_BIG);

	// u_2 =   x_2(FINAL) + u_0,2 
	VecWAXPY(m_coupled_sol_micro,1,m_ext_solver_sol_micro,m_u_0_micro);

	// Add rigid body modes correction, if needed
	if(m_bUsingNullVecs)
	{
		switch (m_RB_modes_system)
		{
			case RBModesSystem::MACRO :
							VecAXPY(m_coupled_sol_BIG,1,m_current_rb_correction);
							break;

			case RBModesSystem::MICRO :
							VecAXPY(m_coupled_sol_micro,1,m_current_rb_correction);
							break;
		}
	}
	
	// Set flags
	m_bCoupled_sols_set = true;
}

IterationStatus FETI_Operations::check_convergence(double rel_residual_conv, double abs_residual_conv, int max_iter_div, double rel_residual_div, double rb_modes_conv)
{
	homemade_assert_msg(m_bCalculatedScalar,"Scalar quantities not calculated yet!");

	m_abs_residual_conv = abs_residual_conv;
	m_rel_residual_conv = rel_residual_conv;
	m_rb_modes_conv = rb_modes_conv;
	m_rel_residual_div = rel_residual_div;
	m_max_iter_div = max_iter_div;

	// Check the iteration divergence
	if(m_kkk == m_max_iter_div) // Iteration divergence
	{
		m_bDivIter = true;
	}

	// Check the residual convergence / divergence
	if(m_current_rho > 0)
	{
		// Check the residual convergence
		if(m_current_rho < m_abs_residual_conv) // Absolute convergence
		{
			m_bConvResidualAbs = true;
		}

		if(m_current_rho < m_rel_residual_conv * m_rho_0) // Relative convergence
		{
			m_bConvResidualRel = true;
		}

		// Check the residual divergence
		if(m_current_rho > m_rel_residual_div * m_rho_0) // Relative divergence
		{
			m_bDivResidualRel = true;
		}
	} else {
		// There's something wrong ...
		m_bDivResidualNeg = true;
	} 

	// Check the rigid body modes correction convergence
	if(m_bUsingNullVecs)
	{
		if(std::abs(m_current_RB_mode_corr - m_previous_RB_mode_corr)/m_current_RB_mode_corr < m_rb_modes_conv) // Correction convergence
		{
			m_bConvRBCorrRel = true;
		}
	} else {
		// Short-circuit this boolean
		m_bConvRBCorrRel = true;
	}

	// Possible cases:
	// - Any diverged flag is true: diverged!
	// - Both the RB modes correction AND one of the residual convergence flags are true: converged!
	// - All flags false: continue iterating!

	m_bConv = ( m_bConvResidualAbs || m_bConvResidualRel ) && m_bConvRBCorrRel;
	m_bDiv = m_bDivResidualNeg || m_bDivResidualRel || m_bDivIter;

	IterationStatus output_status = IterationStatus::ITERATING;
	if(m_bDiv) {
		output_status = IterationStatus::DIVERGED;
	} else if(m_bConv) {
		output_status = IterationStatus::CONVERGED;
	} else {
		output_status = IterationStatus::ITERATING;
	}

	return output_status;
}

}