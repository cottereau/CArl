/*
 * FETI_dyn_operations.h
 *
 *  Created on: Nov 23,2021
 *      Author: Chensheng Luo
 */

#ifndef FETI_DYN_OPERATIONS_H_
#define FETI_DYN_OPERATIONS_H_

#include "carl_headers.h"
#include "common_enums.h"
#include "PETSC_matrix_operations.h"

namespace carl
{

class FETI_Dyn_Operations
{
protected:
    libMesh::Parallel::Communicator& m_comm;

    bool        m_bScratchFolderSet;
    std::string m_scratch_folder_path;      ///< Scratch folder path

    /// Default constructor
  FETI_Dyn_Operations();



  void copy_PETSC_vector(std::string copy_path, std::string paste_path);
public:
  /// Constructor with input parameter and libMesh communicator
  FETI_Dyn_Operations(libMesh::Parallel::Communicator& comm, const std::string& scratch_folder_path) :
    m_comm { comm },
    m_bScratchFolderSet{true},
    m_scratch_folder_path{scratch_folder_path}
  {
  };

  
  /// Destructor, deallocates the PETSc 
  ~FETI_Dyn_Operations()
  {
  };


  // -- setup
    void init_prepare_rhs_vector(std::string& rhs_vector_A_path, 
      std::string& rhs_vector_B_path);
    

  // -- operation

    void prepare_rhs_vector(double beta,
      double deltat,
      std::string force_vector_folder, 
      int index,
      std::string this_acc_path,
      std::string this_speed_path,
      std::string this_displace_path,
      std::string stiffness_matrix_path,
      std::string output_rhs_path);
    
    void Newmark_speed_free(double gamma,
      double deltat,
      std::string prev_acc_path,
      std::string this_acc_path,
      std::string input_speed_path,
      std::string output_speed_path);
    void Newmark_displacement_free(double beta,
      double deltat,
      std::string prev_acc_path,
      std::string this_acc_path,
      std::string prev_speed_path,
      std::string input_displace_path,
      std::string output_displace_path);

    void interpolate_A_acceleration(int jjj,int m);
    void rhs_interpolation(std::string coupling_matrix_A_path,std::string coupling_matrix_B_path);
    void rhs_link(std::string coupling_matrix_path,std::string interpolation_vector_path,std::string output_path);

    void Newmark_speed_link(double gamma,
      double deltat,
      std::string this_acc_path,
      std::string output_speed_path);
    void Newmark_displacement_link(double beta,
      double deltat,
      std::string this_acc_path,
      std::string output_displace_path);

    void add_free_link(std::string free_term_path,
      std::string link_term_path,
      std::string output_path);

  // -- IO
    void output_B_result(std::string result_folder_path,
      int index);

    void move_to_prev_B();

    void output_A_result(std::string result_folder_path,
      int index);

    void move_to_prev_A();

  
};
}

#endif /* FETI_DYN_OPERATIONS_H_ */

