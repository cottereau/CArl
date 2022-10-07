/*
 * \file FETI_operations_IO.cpp
 *
 *  Created on: Apr 23, 2017
 *      Author: Thiago Milanetto Schlittler
 * 
 * \brief **STAT/DYN-CG**  read/write methods for the FETI solver.
 */


#include "FETI_operations.h"

namespace carl
{
//  --- Coupling matrix and preconditioner methods
void FETI_Operations::set_coupling_matrix_R_micro()
{
  homemade_assert_msg(m_bCouplingFolderSet,"Common coupling matrix path not set yet!");

  // Create matrix and, if needed, get sizes
  MatCreate(m_comm.get(),&m_C_R_micro);

  // If the sizes were defined already, set them
  if(m_bmicro_sizes_set && m_bR_sizes_set)
  {
    MatSetSizes(m_C_R_micro,m_C_R_micro_M_local,m_C_R_micro_N_local,m_C_R_micro_M,m_C_R_micro_N);
  }

  // Read the matrix
  read_PETSC_matrix(m_C_R_micro,m_coupling_folder_path + "/coupling_matrix_micro.petscmat",m_comm.get());
  PetscInt N,M;
  MatGetSize(m_C_R_micro,&N,&M);

  // If the micro sizes were not defined yet, define them 
  if(!m_bmicro_sizes_set)
  {
    MatGetLocalSize(m_C_R_micro,&m_C_R_micro_M_local,&m_C_R_micro_N_local);
    MatGetSize(m_C_R_micro,&m_C_R_micro_M,&m_C_R_micro_N);
    m_bmicro_sizes_set = true;
  }

  // If the mediator sizes were not defined yet, define them
  if(!m_bR_sizes_set)
  {
    m_C_RR_M = m_C_R_micro_M; m_C_RR_M_local = m_C_R_micro_M_local;
    m_bR_sizes_set = true;
  }

  // Set null vector dimensions (if needed)
  if(m_bUsingNullVecs && m_RB_modes_system == RBModesSystem::MICRO)
  {
    m_null_vecs_N = m_C_R_micro_N;
    m_null_vecs_N_local = m_C_R_micro_N_local;
    m_bNullVecsDimensionsSet = true;
  }

  // Set flags
  m_bC_R_micro_MatrixSet = true;
  m_bCouplingMatricesSet = m_bC_R_BIG_MatrixSet && m_bC_R_micro_MatrixSet && m_bC_RR_MatrixSet;
}

void FETI_Operations::set_coupling_matrix_R_BIG()
{
  homemade_assert_msg(m_bCouplingFolderSet,"Common coupling matrix path not set yet!");

  // Create matrix and, if needed, get sizes
  MatCreate(m_comm.get(),&m_C_R_BIG);
  
  // If the sizes were defined already, set them
  if(m_bBIG_sizes_set && m_bR_sizes_set)
  {
    MatSetSizes(m_C_R_BIG,m_C_R_BIG_M_local,m_C_R_BIG_N_local,m_C_R_BIG_M,m_C_R_BIG_N);
  }

  // Read the matrix
  read_PETSC_matrix(m_C_R_BIG,m_coupling_folder_path + "/coupling_matrix_macro.petscmat",m_comm.get());

  // If the macro sizes were not defined yet, define them 
  if(!m_bBIG_sizes_set)
  {
    MatGetLocalSize(m_C_R_BIG,&m_C_R_BIG_M_local,&m_C_R_BIG_N_local);
    MatGetSize(m_C_R_BIG,&m_C_R_BIG_M,&m_C_R_BIG_N);
    m_bBIG_sizes_set = true;
  }

  // If the mediator sizes were not defined yet, define them
  if(!m_bR_sizes_set)
  {
    m_C_RR_M = m_C_R_BIG_M; m_C_RR_M_local = m_C_R_BIG_M_local;
    m_bR_sizes_set = true;
  }

  // Set null vector dimensions (if needed)
  if(m_bUsingNullVecs && m_RB_modes_system == RBModesSystem::MACRO)
  {
    m_null_vecs_N = m_C_R_BIG_N;
    m_null_vecs_N_local = m_C_R_BIG_N_local;
    m_bNullVecsDimensionsSet = true;
  }

  // Set flags
  m_bC_R_BIG_MatrixSet = true;
  m_bCouplingMatricesSet = m_bC_R_BIG_MatrixSet && m_bC_R_micro_MatrixSet && m_bC_RR_MatrixSet;
}

void FETI_Operations::set_coupling_matrix_RR()
{
  homemade_assert_msg(m_bCouplingFolderSet,"Common coupling matrix path not set yet!");

  // Create matrix and, if needed, get sizes
  MatCreate(m_comm.get(),&m_C_RR);

  // If the sizes were defined already, set them
  if(m_bR_sizes_set)
  {
    MatSetSizes(m_C_RR,m_C_RR_M_local,m_C_RR_M_local,m_C_RR_M,m_C_RR_M);
  }

  read_PETSC_matrix(m_C_RR,m_coupling_folder_path + "/coupling_matrix_mediator.petscmat",m_comm.get());

  // If the mediator sizes were not defined yet, define them
  if(!m_bR_sizes_set)
  {
    MatGetLocalSize(m_C_RR,&m_C_RR_M_local,NULL);
    MatGetSize(m_C_RR,&m_C_RR_M,NULL);
    m_bR_sizes_set = true;
  }

  // Set flags
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

//  --- Null space / rigid body modes methods
void FETI_Operations::read_null_space_vecs(const std::string& RB_vectors_base, int nb_of_rb_vectors)
{

  /*
   *    Note: for now, this function reads both sets of null space vectors (the "original" ones and
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

//  --- FETI read methods
void FETI_Operations::read_decoupled_solutions()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

  // Create the vectors and, if needed, set sizes
  VecCreate(m_comm.get(),&m_u_0_BIG);
  VecCreate(m_comm.get(),&m_u_0_micro);

  // If the sizes were defined already, set them
  if(m_bBIG_sizes_set)
  {
    VecSetSizes(m_u_0_BIG,m_C_R_BIG_N_local,m_C_R_BIG_N);
  }
  
  if(m_bmicro_sizes_set)
  {
    VecSetSizes(m_u_0_micro,m_C_R_micro_N_local,m_C_R_micro_N);
  }
  
  // Read them
  read_PETSC_vector(m_u_0_BIG,m_scratch_folder_path + "/ext_solver_u0_A_sys_sol_vec.petscvec", m_comm.get());
  read_PETSC_vector(m_u_0_micro,m_scratch_folder_path + "/ext_solver_u0_B_sys_sol_vec.petscvec", m_comm.get());

  // If the sizes were not defined yet, define them 
  if(!m_bBIG_sizes_set)
  {
    VecGetLocalSize(m_u_0_BIG,&m_C_R_BIG_N_local);
    VecGetSize(m_u_0_BIG,&m_C_R_BIG_N);
    m_bBIG_sizes_set = true;
  }

  if(!m_bmicro_sizes_set)
  {
    VecGetLocalSize(m_u_0_micro,&m_C_R_micro_N_local);
    VecGetSize(m_u_0_micro,&m_C_R_micro_N);
    m_bmicro_sizes_set = true;
  }

  // Set null vector dimensions (if needed)
  if(m_bUsingNullVecs && !m_bNullVecsDimensionsSet)
  {
    switch (m_RB_modes_system)
    {
      case RBModesSystem::MACRO :
              m_null_vecs_N = m_C_R_BIG_N;
              m_null_vecs_N_local = m_C_R_BIG_N_local;
              break;

      case RBModesSystem::MICRO :
              m_null_vecs_N = m_C_R_micro_N;
              m_null_vecs_N_local = m_C_R_micro_N_local;
              break;
    }
    m_bNullVecsDimensionsSet = true;
  }

  // Set the flag
  m_bSet_u_0 = true;
}

void FETI_Operations::read_ext_solver_output()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

  // Create the vectors
  VecCreate(m_comm.get(),&m_ext_solver_sol_BIG);
  VecCreate(m_comm.get(),&m_ext_solver_sol_micro);

  // If the sizes were defined already, set them
  if(m_bBIG_sizes_set)
  {
    VecSetSizes(m_ext_solver_sol_BIG,m_C_R_BIG_N_local,m_C_R_BIG_N);
  }
  
  if(m_bmicro_sizes_set)
  {
    VecSetSizes(m_ext_solver_sol_micro,m_C_R_micro_N_local,m_C_R_micro_N);
  }

  // Read them
  read_PETSC_vector(m_ext_solver_sol_BIG,m_scratch_folder_path + "/ext_solver_A_sys_sol_vec.petscvec", m_comm.get());
  read_PETSC_vector(m_ext_solver_sol_micro,m_scratch_folder_path + "/ext_solver_B_sys_sol_vec.petscvec", m_comm.get());

  // If the sizes were not defined yet, define them 
  if(!m_bBIG_sizes_set)
  {
    VecGetLocalSize(m_ext_solver_sol_BIG,&m_C_R_BIG_N_local);
    VecGetSize(m_ext_solver_sol_BIG,&m_C_R_BIG_N);
    m_bBIG_sizes_set = true;
  }

  if(!m_bmicro_sizes_set)
  {
    VecGetLocalSize(m_ext_solver_sol_micro,&m_C_R_micro_N_local);
    VecGetSize(m_ext_solver_sol_micro,&m_C_R_micro_N);
    m_bmicro_sizes_set = true;
  }

  // Set null vector dimensions (if needed)
  if(m_bUsingNullVecs && !m_bNullVecsDimensionsSet)
  {
    switch (m_RB_modes_system)
    {
      case RBModesSystem::MACRO :
              m_null_vecs_N = m_C_R_BIG_N;
              m_null_vecs_N_local = m_C_R_BIG_N_local;
              break;

      case RBModesSystem::MICRO :
              m_null_vecs_N = m_C_R_micro_N;
              m_null_vecs_N_local = m_C_R_micro_N_local;
              break;
    }
    m_bNullVecsDimensionsSet = true;
  }

  // Set the flag
  m_bSet_ext_solver_sol = true;
}

void FETI_Operations::read_rb_corr()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bNullVecsDimensionsSet,"Null vectors sizes not set yet!");

  // Create the vectors
  VecCreate(m_comm.get(),&m_current_rb_correction);
  VecSetSizes(m_current_rb_correction,m_null_vecs_N_local,m_null_vecs_N);
  read_PETSC_vector(m_current_rb_correction,m_scratch_folder_path + "/FETI_RB_correction.petscvec", m_comm.get());

  // Set flag
  m_bSet_current_RB_correction = true;
}

void FETI_Operations::read_previous_phi()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bR_sizes_set,"Mediator space dimensions not set yet!");

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
  homemade_assert_msg(m_bR_sizes_set,"Mediator space dimensions not set yet!");

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
  homemade_assert_msg(m_bR_sizes_set,"Mediator space dimensions not set yet!");

  // Create p(0 ... kkk),
  m_previous_p_ptr = new Vec[m_kkk+1];
  VecDuplicateVecs(m_previous_residual,m_kkk+1,&m_previous_p_ptr);

  // Read p(0 ... kkk)
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
  homemade_assert_msg(m_bR_sizes_set,"Mediator space dimensions not set yet!");

  // Create q(0 ... kkk)
  m_previous_q_ptr = new Vec[m_kkk+1];
  VecDuplicateVecs(m_previous_residual,m_kkk+1,&m_previous_q_ptr);

  // Read q(0 ... kkk - 1)
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
   *          'phi(kkk)'
   *          'p(0 ... kkk)'
   *          'q(0 ... kkk - 1)'
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

//  --- Write methods
void FETI_Operations::export_ext_solver_rhs_Ct_p()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

  homemade_assert_msg(m_bSet_current_p,"Current 'p' not calculated yet!");
  this->export_ext_solver_rhs(m_current_p);
}

void FETI_Operations::export_ext_solver_rhs_initial()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");

  homemade_assert_msg(m_bSet_current_z,"Current 'p' not calculated yet!");
  this->export_ext_solver_rhs(m_current_z);
}

void FETI_Operations::export_ext_solver_rhs_Ct_phi()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSet_current_phi,"Current 'phi' not calculated yet!");

  this->export_ext_solver_rhs(m_current_phi);
}

void FETI_Operations::export_rb_correction_vector()
{ 
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSet_current_RB_correction,"Current RB correction not calculated yet!");

  
  write_PETSC_vector(m_current_rb_correction,m_scratch_folder_path + "/FETI_RB_correction.petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_current_rb_correction,m_scratch_folder_path + "/FETI_RB_correction.m",m_comm.get());
#endif
}

void FETI_Operations::export_p()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSet_current_p,"Current 'p' not calculated yet!");

  write_PETSC_vector(m_current_p,m_scratch_folder_path + "/FETI_iter__p__" + std::to_string(m_kkk+1) + ".petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_current_p,m_scratch_folder_path + "/FETI_iter__p__" + std::to_string(m_kkk+1) + ".m",m_comm.get());
#endif
}

void FETI_Operations::export_q()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSet_previous_q_ptr,"'q' vectors not calculated yet!");

  write_PETSC_vector(m_previous_q_ptr[m_kkk],m_scratch_folder_path + "/FETI_iter__q__" + std::to_string(m_kkk) + ".petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_previous_q_ptr[m_kkk],m_scratch_folder_path + "/FETI_iter__q__" + std::to_string(m_kkk) + ".m",m_comm.get());
#endif
}

void FETI_Operations::export_r()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSet_current_residual,"Current residual not calculated yet!");

  write_PETSC_vector(m_current_residual,m_scratch_folder_path + "/FETI_iter__r__current.petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_current_residual,m_scratch_folder_path + "/FETI_iter__r__current.m",m_comm.get());
#endif
}

void FETI_Operations::export_phi()
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSet_current_phi,"Current 'phi' not calculated yet!");

  write_PETSC_vector(m_current_phi,m_scratch_folder_path + "/FETI_iter__phi__current.petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_current_phi,m_scratch_folder_path + "/FETI_iter__phi__current.m",m_comm.get());
#endif
}

void FETI_Operations::export_inital_vecs()
{
  // Export r(0) and p(0) ( p(0) is identical to z(0))
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bSet_current_residual,"Current residual not calculated yet!");

  // In all cases, print r(0)
  write_PETSC_vector(m_current_residual,m_scratch_folder_path + "/FETI_iter__r__current.petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_current_residual,m_scratch_folder_path + "/FETI_iter__r__current.m",m_comm.get());
#endif

  // p(0) is only identical to r(0) if neither a preconditioner or the RB modes are used
  if(m_precond_type != BaseCGPrecondType::NO_PRECONDITIONER || m_bUsingNullVecs)
  {
    write_PETSC_vector(m_current_z,m_scratch_folder_path + "/FETI_iter__p__0.petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
    write_PETSC_vector_MATLAB(m_current_z,m_scratch_folder_path + "/FETI_iter__p__0.m",m_comm.get());
#endif
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
    VecNorm(m_current_rb_correction,NORM_2,&RB_correct);
  }

  // ONLY write in proc 0!
  if(m_comm.rank() == 0)
  {
    std::ofstream scalar_data;

    // Export the scalar data
    scalar_data.open(m_scratch_folder_path + "/FETI_iter_scalar_data.dat");
    scalar_data.precision(15);
    scalar_data << m_kkk << " " << residual << " " << residual;

    if(m_bUsingNullVecs)
    {
      scalar_data <<  " " << RB_correct;
    }

    scalar_data << std::endl;

    scalar_data.close();

    // Export the convergence data
    scalar_data.open(m_scratch_folder_path + "/FETI_convergence.dat");
    scalar_data.precision(15);
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
    scalar_data.precision(15);
    scalar_data << m_kkk + 1 << " " << m_rho_0 << " " << m_current_rho;

    if(m_bUsingNullVecs)
    {
      scalar_data <<  " " << m_current_RB_mode_corr;
    }

    scalar_data << std::endl;

    scalar_data.close();

    // Export the convergence data
    scalar_data.open(m_scratch_folder_path + "/FETI_convergence.dat",std::ofstream::app);
    scalar_data.precision(15);
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
    scalar_data.precision(15);
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

//[STAT]Export the final coupled solution
void FETI_Operations::export_coupled_solution(std::string output_base)
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bCoupled_sols_set,"Coupled solutions not calculated yet!");

  write_PETSC_vector(m_coupled_sol_BIG,output_base + "/coupled_sol_A.petscvec",m_comm.rank(),m_comm.get());
  write_PETSC_vector(m_coupled_sol_micro,output_base + "/coupled_sol_B.petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_coupled_sol_BIG,output_base + "/coupled_sol_A.m",m_comm.get());
  write_PETSC_vector_MATLAB(m_coupled_sol_micro,output_base + "/coupled_sol_B.m",m_comm.get());
#endif
}

//[DYN-CG]Export the final coupled solution
void FETI_Operations::export_dynamic_coupled_solution(std::string output_base)
{
  homemade_assert_msg(m_bScratchFolderSet,"Scratch folder not set yet!");
  homemade_assert_msg(m_bCoupled_sols_set,"Coupled solutions not calculated yet!");

  write_PETSC_vector(m_coupled_sol_BIG,output_base + "/this_acc_A_link_sys_sol_vec.petscvec",m_comm.rank(),m_comm.get());
  write_PETSC_vector(m_coupled_sol_micro,output_base + "/this_acc_B_link_sys_sol_vec.petscvec",m_comm.rank(),m_comm.get());

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  write_PETSC_vector_MATLAB(m_coupled_sol_BIG,output_base + "/this_acc_A_link_sys_sol_vec.m",m_comm.get());
  write_PETSC_vector_MATLAB(m_coupled_sol_micro,output_base + "/this_acc_B_link_sys_sol_vec.m",m_comm.get());
#endif
}

void FETI_Operations::print_previous_iters_conv(int nb_of_iters)
{
  std::cout << std::endl;
  std::cout << " --> Iteration no. " << m_kkk + 1 << " : ";

  if(m_bConvResidualAbs)
  {
    std::cout << "   > Abs. residual convergence : rho(kkk+1) < " << m_abs_residual_conv << std::endl;
  }
  if(m_bConvResidualRel)
  {
    std::cout << "   > Rel. residual convergence : rho(kkk+1) < " << m_rel_residual_conv << " * rho(0) " << std::endl;
  }
  if(m_bDivResidualRel)
  {
    std::cout << "   > Rel. residual DIVERGENCE : rho(kkk+1) > " << m_rel_residual_div << " * rho(0) " << std::endl;
  }
  if(m_bDivResidualNeg)
  {
    std::cout << "   > Negative residual DIVERGENCE" << std::endl;
  }
  if(m_bConvResidualAbs)
  {
    std::cout << "   > Iter. DIVERGENCE : kkk + 1 > " << m_max_iter_div << std::endl;
  }
  if(m_bUsingNullVecs)
  { 
    if(m_bConvRBCorrRel)
    {
      std::cout << "   > Rel. RB mode correction convergence : abs( RB(kkk+1) - RB_(kkk) ) / RB(kkk+1) < " << m_rb_modes_conv << std::endl;
    }
  }
  std::cout << std::endl;

  std::cout << " --> Previous " << std::min(m_kkk+2,nb_of_iters) << " iterations convergence parameter : " << std::endl;
  if(m_bUsingNullVecs)
  { 
    std::cout << "[ kkk ] [ rho(kkk) ] [ | RB_corr(kkk) | ]" << std::endl;
  }
  else
  {
    std::cout << "[ kkk ] [ rho(kkk) ]" << std::endl;
  }

  std::string command_string = "tail -n " + std::to_string(nb_of_iters) + " " + m_scratch_folder_path + "/FETI_convergence.dat";

  if(m_comm.rank() == 0)
  {
    std::cout << exec_command(command_string) << std::endl;
  }
  if( m_bConv ) {
    std::cout << " --> Converged!" << std::endl;
  } else if ( m_bDiv ) {
    std::cout << " --> DIVERGED!" << std::endl;
  } else {
    std::cout << " --> Iterating ..." << std::endl;
  }
  std::cout << std::endl;
}

}
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */
