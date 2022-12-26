/*
 * FETI_dyn_operations.h
 *
 *  Created on: Nov 23,2021
 *      Author: Chensheng Luo
 */

#ifndef FETI_DYN_OPERATIONS_H_
#define FETI_DYN_OPERATIONS_H_

#include "carl_headers.h"
#include "PETSC_matrix_operations.h"
#include "carl_loop_dyn_force_parser.h"
#include "newmark_param_parser.h"

#include <ctime>

namespace carl
{

/** \brief **DYN** Class containing the operations needed for the Dynamic solver.
 *
 *  This class is used by the several `CArl_dyn` programs to do operations of the Dynamic solver,
 *  including matrix and vector I/O and iteration operations. Due to the need to read vectors and 
 *  matrices and the usage of several PETSc operations for which the libMesh interface was not 
 *  implemented, direct PETSc `Vec`'s and `Mat`'s are used instead of their libMesh interfaces.
 *  libMesh's parallel communicators are still used, though.
 *
 *
 *  This class methods' implementations are separated into three different files: 
 *  - FETI_dyn_operations_setup.cpp, for setup calculations.
 *  - FETI_dyn_operations_operations.cpp, for FETI operations during the calculate, mainly for Newmark calculation and rhs preparation.
 *  - FETI_dyn_operations_IO.cpp, for result output, copy/paste and others.
 *
 */
class FETI_Dyn_Operations
{
protected:
  libMesh::Parallel::Communicator& m_comm;

  bool        m_bScratchFolderSet;
  std::string m_scratch_folder_path;      ///< Scratch folder path
  bool        m_bResultFolderSet;
  std::string m_result_folder_path;      ///< Result folder path


  bool        m_bAMovedToPrev;
  bool        m_bBMovedToPrev;
  bool        m_bForceAPrepared;
  bool        m_bForceBPrepared;

  
  /// Default constructor
  FETI_Dyn_Operations();

  // --IO
  void copy_PETSC_vector(std::string copy_path, std::string paste_path);                        ///< To easily copy a PETSC vector
  void scale_copy_PETSC_matrix(std::string copy_path, std::string paste_path,PetscScalar scale);///< To easily copy a PETSC vector by scaling
  void scale_copy_PETSC_vector(std::string copy_path, std::string paste_path,PetscScalar scale);///< To easily copy a PETSC matrix by scaling
  void set_NAN_vector(std::string& vector_path);                                                ///< To set a vector to NAN(thus void)


  // --operation
  void add_free_link_one(std::string& free_term_path,
      std::string& link_term_path,
      std::string& output_path);

  // --setup
  void prepare_force_vector_by_modal_and_constant(std::string force_path,
    carl::DynSystemVectorPath* vectors,
    double amplitude);      ///< Calculate  \f$ \text{modal} \times \text{Amplitude} \f$ for all moment, put the result in a folder

  void prepare_force_vector_by_modal_and_sinus(std::string force_path,
    carl::DynSystemVectorPath* vectors,
    double amplitude,
    double frequency,
    double initialPhase,
    double small_deltat,
    int index,
    int timestep);      ///< Calculate \f$ \text{modal} \times \text{Amplitude}\sin(2\pi\times\text{Frequency}\times t+\text{IniitalPhase}) \f$ for all moment, put the result in a folder
  
  void prepare_force_vector_by_modal_and_slope(std::string force_path,
    carl::DynSystemVectorPath* vectors,
    double slope,
    double saturation,
    double offset,
    double small_deltat,
    int index,
    int timestep);     ///< Calculate  \f$ \text{modal} \times \text{Slope} \times t \f$ for all moment, put the result in a folder

public:

  struct carl::DynSystemVectorPath vector_A; ///< All scratch vectors of A
  struct carl::DynSystemVectorPath vector_B; ///< All scratch vectors of A
  std::string m_rhs_interpolation; ///< Interpolation matrix path \f% H \f$(even CG solver is applied, but of no sense)
  std::string m_coupling_vector; ///< Coupling vector \f$ \Lambda \f$ path (even CG solver is applied, but of no sense)

  /// Constructor with input parameter and libMesh communicator
  FETI_Dyn_Operations(libMesh::Parallel::Communicator& comm,
   const std::string scratch_folder_path, 
   const std::string result_folder_path):
    m_comm { comm },
    m_bScratchFolderSet{true},
    m_scratch_folder_path{scratch_folder_path},
    m_bAMovedToPrev{false},
    m_bBMovedToPrev{false},
    m_bForceAPrepared{false},
    m_bForceBPrepared{false},
    m_result_folder_path{result_folder_path},
    m_bResultFolderSet{true},
    m_rhs_interpolation{m_scratch_folder_path+"/rhs_interpolation_vec.petscvec"},
    m_coupling_vector{m_scratch_folder_path+"/coupling_interpolation_vec_sys_sol_vec.petscvec"}
  {
    vector_A.prev_acc = m_scratch_folder_path+"/prev_acc_A.petscvec";
    vector_A.prev_disp = m_scratch_folder_path+"/prev_disp_A.petscvec";
    vector_A.prev_disp_free = m_scratch_folder_path+"/prev_disp_A_free.petscvec";
    vector_A.prev_speed = m_scratch_folder_path+"/prev_speed_A.petscvec";
    vector_A.prev_speed_free = m_scratch_folder_path+"/prev_speed_A_free.petscvec";

    vector_A.rhs_free = m_scratch_folder_path+"/rhs_vec_A_free.petscvec";
    vector_A.rhs_link = m_scratch_folder_path+"/rhs_vec_A_link.petscvec";

    vector_A.this_acc = m_scratch_folder_path+"/this_acc_A.petscvec";
    vector_A.this_acc_free = m_scratch_folder_path+"/this_acc_A_free_sys_sol_vec.petscvec";
    vector_A.this_acc_link = m_scratch_folder_path+"/this_acc_A_link_sys_sol_vec.petscvec";
    vector_A.this_disp = m_scratch_folder_path+"/this_disp_A.petscvec";
    vector_A.this_disp_free = m_scratch_folder_path+"/this_disp_A_free.petscvec";
    vector_A.this_disp_link = m_scratch_folder_path+"/this_disp_A_link.petscvec";
    vector_A.this_speed = m_scratch_folder_path+"/this_speed_A.petscvec";
    vector_A.this_speed_free = m_scratch_folder_path+"/this_speed_A_free.petscvec";
    vector_A.this_speed_link = m_scratch_folder_path+"/this_speed_A_link.petscvec";

    vector_A.inter_disp_free =
    m_scratch_folder_path+"/inter_disp_A_free.petscvec";
    vector_A.coupling_sign = -1; vector_A.this_force =
    m_scratch_folder_path+"/this_force_A.petscvec"; vector_A.next_force =
    m_scratch_folder_path+"/next_force_A.petscvec";

    vector_B.prev_acc = m_scratch_folder_path+"/prev_acc_B.petscvec";
    vector_B.prev_disp = m_scratch_folder_path+"/prev_disp_B.petscvec";
    vector_B.prev_disp_free = m_scratch_folder_path+"/prev_disp_B_free.petscvec";
    vector_B.prev_speed = m_scratch_folder_path+"/prev_speed_B.petscvec";
    vector_B.prev_speed_free = m_scratch_folder_path+"/prev_speed_B_free.petscvec";

    vector_B.rhs_free = m_scratch_folder_path+"/rhs_vec_B_free.petscvec";
    vector_B.rhs_link = m_scratch_folder_path+"/rhs_vec_B_link.petscvec";

    vector_B.this_acc = m_scratch_folder_path+"/this_acc_B.petscvec";
    vector_B.this_acc_free = m_scratch_folder_path+"/this_acc_B_free_sys_sol_vec.petscvec";
    vector_B.this_acc_link = m_scratch_folder_path+"/this_acc_B_link_sys_sol_vec.petscvec";
    vector_B.this_disp = m_scratch_folder_path+"/this_disp_B.petscvec";
    vector_B.this_disp_free = m_scratch_folder_path+"/this_disp_B_free.petscvec";
    vector_B.this_disp_link = m_scratch_folder_path+"/this_disp_B_link.petscvec";
    vector_B.this_speed = m_scratch_folder_path+"/this_speed_B.petscvec";
    vector_B.this_speed_free = m_scratch_folder_path+"/this_speed_B_free.petscvec";
    vector_B.this_speed_link = m_scratch_folder_path+"/this_speed_B_link.petscvec";

    vector_B.this_force = m_scratch_folder_path+"/this_force_B.petscvec";
    vector_B.next_force = m_scratch_folder_path+"/next_force_B.petscvec";

    vector_B.inter_disp_free = "";
    vector_B.coupling_sign = 1;
  };

  
  /// Destructor, deallocates the PETSc 
  ~FETI_Dyn_Operations()
  {
  };

  // -- setup
    void init_prepare_force(int prepare_mode, 
      std::string& input_file, 
      int innerTimes, 
      double small_deltat,
      carl::DynSystemVectorPath* vectorsA,
      carl::DynSystemVectorPath* vectorsB);          ///< Prepare force vectors in the folders /force_*/ of scratch folder

    void set_initial_condition(carl::DynInitialVectorPath* initial_A,
      carl::DynInitialVectorPath* initial_B,
      carl::DynSystemMatrixPath* matrix_A,
      carl::DynSystemMatrixPath* matrix_B,
      carl::NewmarkParams* newmark_A,
      carl::NewmarkParams* newmark_B,
      carl::DynSystemVectorPath* vectors_A,
      carl::DynSystemVectorPath* vectors_B,
      std::string result_folder_path);   ///< Set initial displacement/speed condition
    //TODO: ADD ACCELERATION INITIAL CONDITION!


  // -- operation

    void rhs_free(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors,
      carl::DynSystemMatrixPath* matrices);   ///< Calculate the rhs vector for free equation\f$ F-K (^{P}U) \f$
    
    void Newmark_speed_free(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors);  ///< Calculate the free speed term with Newmark Method \f$ \dot{U}^k_{j,free}=^{P}U^k_{j-1}+\gamma^k \Delta t_k \ddot{U}^k_{j,free},(k= A, B) \f$

    void Newmark_displacement_free(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors);   ///< Calculate the free displacement term with Newmark Method\f$ U^k_{j,free}=^{P}U^k_{j-1}+\beta^k \Delta t_k^{2} \ddot{U}^k_{j,free},(k= A, B) \f$

    void interpolate_A_disp(int jjj,int m);   ///< Calculate the interpolation  of A free speed
    
    void rhs_interpolation(std::string& coupling_matrix_A_path,
      std::string& coupling_matrix_B_path);   ///< Calculate the RHS of equation \f$H\Lambda_{j}=C^{A}\dot{U}^{A}_{j,free}-C^{B}\dot{U}^{B}_{j,free} \f$
    
    void rhs_link(carl::DynSystemVectorPath* vectors,
      std::string coupling_vector_path,
      carl::DynSystemMatrixPath* matrices);   ///< Calculate the RHS of equation \f$\ddot{U}^k_{j,link}=+-(C^k)^{T}\Lambda_{j},(k= A, B) \f$

    void Newmark_speed_link(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors); ///< Calculate the link speed term with Newmark Method \f$ \dot{U}^k_{j,link}= \gamma^k \Delta t_k \ddot{U}^k_{j,link},(k= A, B)  \f$

    void Newmark_displacement_link(carl::NewmarkParams* newmark_params,
      carl::DynSystemVectorPath* vectors); ///< Calculate the link displacement term with Newmark Method \f$ U^k_{j,link}= \beta^k  \Delta t_k^{2} \ddot{U}^k_{j,link},(k= A, B) \f

    void add_free_link(carl::DynSystemVectorPath* vectors); ///< Add up free and link term \f$ \cdot=\cdot_{free}+\cdot_{link} \f$

    void init_test_coupling();
    
    void test_coupling(std::string& coupling_matrix_A_path,
      std::string& coupling_matrix_B_path,
      int index);    ///< Test if output D.O.F. vector verify the equation \f$ C^A U^A-C^B U^B=0 \f$

    void prepare_A_next_force(carl::DynSystemVectorPath* vectors,
      int index,
      int step,
      int prepare_mode,
      std::string& force_prepare_file,
      double small_deltat);  ///<prepare A force to get ready for rhs_free

    void prepare_B_next_force(carl::DynSystemVectorPath* vectors,
      int index,
      int prepare_mode,
      std::string& force_prepare_file,
      double small_deltat);  ///<prepare B force to get ready for rhs_free

  // -- IO
    void output_A_result(int index);             ///< Copy the acceleration/speed/displacement of `this_A_*.petscvec` to result folder

    void move_to_prev_A();   ///< Move all `this_A_*.petscvec` file to `prev_A_* file.petscvec`

    void delete_A_this_vector();  ///< Put all "this_A_*.petscvec" vectors to NAN

    void output_B_result(int index);            ///< Copy the acceleration/speed/displacement of `this_B_*.petscvec` to result folder

    void move_to_prev_B();    ///< Move all  `this_B_*.petscvec` file to `prev_B_*.petscvec` file
    
    void delete_B_this_vector();  ///< Put all "this_B_*.petscvec" vectors to NAN

    void delete_coupling_vector(); ///< Put all coupling vectors to NAN

    // void prepare_CG_free_result(std::string& depature_path,
    //   std::string destination_path);   ///< [DYN-CG]Move all free solution of this step to prepare CG solver

    // void prepare_CG_scaled_matrix(std::string& M_path_A,std::string& M_path_B,
    //   std::string dest_path,carl::NewmarkParams* newmark_A,carl::NewmarkParams* newmark_B);      ///< [DYN-CG]Scale mass matrix to adapt static CG solver, with \f$K^k \RightArrow \frac{1}{\beta^k(\Delta t^k)^2} \tilde{M}^k\f$

    // void prepare_CG_scaled_vector(std::string& rb_path,int nb_rb,std::string forth_path,
    //   std::string destination_path,carl::NewmarkParams* newmark); ///< [DYN-CG]Scale force and rigid body vector\f$R^B \RightArrow \beta^B(\Delta t^B)^2 R^B \f$

    void export_calculation_time(int index,std::string stage);   ///< Output current moment at the file Time_data.dat
  
};
}

#endif /* FETI_DYN_OPERATIONS_H_ */

