/*
 * \file FETI_operations_setup.cpp
 *
 *  Created on: Apr 23, 2017
 *      Author: Thiago Milanetto Schlittler
 * 
 * \brief **STAT/DYN-CG**   setup and internal calculations (including all `protected` methods)  for the FETI solver.
 */

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

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_coupling_jacobi_precond_vec,m_scratch_folder_path + "/precond_Jacobi_vector.m",m_comm.get());
#endif

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
  // -> Calculate dummy_seq_vec_bis on the first proc, and then broadcast the value
  
  /*    
   *    Originally, this operation was done locally, but due to a syncing issue,
   *    we have to do it this way to avoid a "Value must the same in all processors" error
   *    when calling VecMAXPY below.
   */ 
  PETSC_MatMultScale_Bcast(m_inv_RITRI_mat,dummy_seq_vec,dummy_seq_vec_bis,-1);

  m_comm.barrier();

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
  write_PETSC_vector(vec_C_micro_t_p_kkk_PETSc,m_scratch_folder_path + "/ext_solver_B_rhs.petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(vec_C_micro_t_p_kkk_PETSc,m_scratch_folder_path + "/ext_solver_B_rhs.m",m_comm.get());
  write_PETSC_vector_MATLAB(vec_C_BIG_t_p_kkk_PETSc,m_scratch_folder_path + "/ext_solver_A_rhs.m",m_comm.get());
#endif

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
  if(m_bCoupled_sols_set)
  {
    VecDestroy(&m_coupled_sol_micro);
    VecDestroy(&m_coupled_sol_BIG);
  }
}

//  --- Public methods
//  --- Coupling matrix and preconditioner methods
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
  // m_null_vecs         : n_sys x nb_of_vecs ( R )  -> nb_of_vecs vectors of dim n_sys
  // m_null_coupled_vecs : n_med x nb_of_vecs ( RC ) -> nb_of_vecs vectors of dim n_coupl

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
    #ifdef PRINT_MATLAB_DEBUG
      carl::write_PETSC_matrix_MATLAB(m_inv_RITRI_mat,m_scratch_folder_path + "/rb_inv_RITRI.m",PETSC_COMM_SELF);
      carl::write_PETSC_matrix_MATLAB(m_RITRI_mat,m_scratch_folder_path + "/rb_RITRI.m",PETSC_COMM_SELF);
    #endif
  }

  // Set up flag
  m_bNullVecsSet = true;
  m_binvRITRIMatSet = true;

  // Cleanup
  MatDestroy(&m_RITRI_mat);
}

}
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */
