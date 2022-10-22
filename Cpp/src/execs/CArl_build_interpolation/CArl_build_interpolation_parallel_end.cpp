/*
 *
 *  Created on: Oct 12, 2022
 *      Author: Chensheng Luo
 */

#include "CArl_build_interpolation.h"
/** \file CArl_build_interpolation_parallel_end.cpp
 * \brief **DYN-DI** Program used to calculate interpolation matrix by parallel calculation with.

 Usage: `./CArl_build_interpolation_matrix_parallel_begin -i [input file]`

The program CArl_build_interpolation_parallel_begin and CArl_build_interpolation_parallel_end will calculate the interpolation matrix \f$H=\beta^A (\Delta t^A)^2 C^A (M^A)^{-1} (C^A)^T+\beta^B (\Delta t^B)^2 C^B (M^B)^{-1} (C^B)^T \f$. 
This calculation is done in following three steps:
  1. Separating \f$ (C^i)^T \f$ into several columns: C^i_{ttt} (begin)
  2. Inverting by solving \f$ M^i X^i_{ttt}=C^i_{ttt} \f$ (solving)
  3. By putting \f$X^i_{ttt}\f$ together to form a matrix \f$G^i\f$ and then calculate \f$H=\beta^A (\Delta t^A)^2 C^A G^A+\beta^B (\Delta t^B)^2 C^B G^B \f (end)$

The input file is parsed by the get_input_params(GetPot& field_parser, carl_build_interpolation_params& input_params) function. 

 */

int main(int argc, char *argv[]) 
{
  // --- Setup part
  // - Initialize libMesh
  libMesh::LibMeshInit init(argc, argv);

  // libMesh's C++ / MPI communicator wrapper
  libMesh::Parallel::Communicator& WorldComm = init.comm();

  // Number of processors and processor rank.
  int rank = WorldComm.rank();
  int nodes = WorldComm.size();

  libMesh::PerfLog perf_log("Main program");

  // -  Set up input_params
  // Command line parser
  GetPot command_line(argc, argv);

  // File parser
  GetPot field_parser;

  // If there is an input file, parse it to get the parameters. Else, parse the command line
  std::string input_filename;
  if (command_line.search(2, "--inputfile", "-i")) 
  {
    input_filename = command_line.next(input_filename);
    field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
  } 
  else
  {
    field_parser = command_line;
  }

  // Parse file content to the input_params
  carl_build_interpolation_params input_params;
  get_input_params(field_parser, input_params);

  if( input_params.m_bParallelPossible == false){
    homemade_error_msg("[CArl Error] Parallel Not Possible, check your parameter!");
  }

  // - Load Vectors & Form Matrices

  perf_log.push("Load Vectors & Form Matrices");

  Vec cm_vector_A,cm_vector_B;
  Mat mass_invert_couplingt_A,mass_invert_couplingt_B;
  PetscInt CmA,CmB,CnA,NCmA,CnB,NCmB;

  Mat coupling_matrix_A;
  Mat coupling_matrix_B;

  Mat H_1;
  Mat H_2;

  VecCreate(PETSC_COMM_SELF,&cm_vector_A);
  VecCreate(PETSC_COMM_SELF,&cm_vector_B);
  MatCreate(WorldComm.get(),&coupling_matrix_A);
  MatCreate(WorldComm.get(),&coupling_matrix_B);

  carl::read_PETSC_vector(cm_vector_A,
    input_params.output_base + "mass_inverted_coupling_A/MCTA_"+std::to_string(0)+".petscvec_sys_sol_vec.petscvec",PETSC_COMM_SELF);
  carl::read_PETSC_vector(cm_vector_B,
    input_params.output_base + "mass_inverted_coupling_B/MCTB_"+std::to_string(0)+".petscvec_sys_sol_vec.petscvec",PETSC_COMM_SELF);
  carl::read_PETSC_matrix(coupling_matrix_A, 
                input_params.path_macro_coupling_matrix, WorldComm.get());
  carl::read_PETSC_matrix(coupling_matrix_B, 
                input_params.path_micro_coupling_matrix, WorldComm.get());

  // Examination of dimension
  

  VecGetSize(cm_vector_A, &CmA);
  VecGetSize(cm_vector_B, &CmB);
  MatGetSize(coupling_matrix_A, &CnA, &NCmA);
  MatGetSize(coupling_matrix_B, &CnB, &NCmB);
  if (NCmA != CmA || NCmB != CmB || CnA != CnB){
    homemade_error_msg("[CArl Error] Dimension Error! The coupling matrix and the mass matrix are not compatible for matrix product.");
  }

  // Prepare matrices
  MatCreateAIJ(WorldComm.get(), PETSC_DECIDE, PETSC_DECIDE, CmA, CnA, CmA, NULL, CnA, NULL, &mass_invert_couplingt_A);
  MatCreateAIJ(WorldComm.get(), PETSC_DECIDE, PETSC_DECIDE, CmB, CnB, CmB, NULL, CnB, NULL, &mass_invert_couplingt_B);

  PetscInt i;
  int j;
  PetscInt    row_indiceA[CmA];
  PetscScalar rowA[CmA];
  for(j=0 ; j<CmA ; j++){ row_indiceA[j] = j;}

  for(i=0 ; i<CnA ; i++){
    carl::read_PETSC_vector(cm_vector_A,
        input_params.output_base + "mass_inverted_coupling_A/MCTA_"+std::to_string(i)+".petscvec_sys_sol_vec.petscvec",PETSC_COMM_SELF);
    VecGetValues(cm_vector_A,CmA, row_indiceA, rowA);
    MatSetValues(mass_invert_couplingt_A,CmA,row_indiceA,1,&i,rowA,INSERT_VALUES);
  }
  MatAssemblyBegin(mass_invert_couplingt_A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mass_invert_couplingt_A, MAT_FINAL_ASSEMBLY);

  PetscInt    row_indiceB[CmB];
  PetscScalar rowB[CmB];
  for(j=0 ; j<CmB ; j++){ row_indiceB[j] = j;}

  for(i=0 ; i<CnB ; i++){
    carl::read_PETSC_vector(cm_vector_B,
        input_params.output_base + "mass_inverted_coupling_B/MCTB_"+std::to_string(i)+".petscvec_sys_sol_vec.petscvec",PETSC_COMM_SELF);
    VecGetValues(cm_vector_B,CmB, row_indiceB, rowB);
    MatSetValues(mass_invert_couplingt_B,CmB,row_indiceB,1,&i,rowB,INSERT_VALUES);
  }
  MatAssemblyBegin(mass_invert_couplingt_B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mass_invert_couplingt_B, MAT_FINAL_ASSEMBLY);

  perf_log.pop("Load Vectors & Form Matrices");

  // - Product coupling and mass matrix

  // Create result matrix
  MatCreateAIJ(WorldComm.get(),PETSC_DECIDE,PETSC_DECIDE,CnA,CnA,CnA,NULL,CnA,NULL,&H_1);
  MatCreateAIJ(WorldComm.get(),PETSC_DECIDE,PETSC_DECIDE,CnB,CnB,CnB,NULL,CnB,NULL,&H_2);

  // Matrix product
  perf_log.push("STEP1: Matrix product");
  //H_1=C_A*MCT_A
  MatMatMult(coupling_matrix_A,
        mass_invert_couplingt_A,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &H_1);

  //H_2=C_B*MCT_B
  MatMatMult(coupling_matrix_B,
        mass_invert_couplingt_B,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &H_2);
  perf_log.pop("STEP1: Matrix product");
              
  // - Linear scale and addition
  perf_log.push("STEP2: Linear scale and addition");
  //H_2=H_2*alpha_B
  MatScale(H_2,input_params.newmark_B.deltat*input_params.newmark_B.deltat*input_params.newmark_B.beta);

  //H_1=H_1*alpha_A+H_2
  MatAYPX(H_1,input_params.newmark_A.deltat*input_params.newmark_A.deltat*input_params.newmark_A.beta,H_2, SAME_NONZERO_PATTERN);
  Mat* sys_mat_PETSC_H_save;
  sys_mat_PETSC_H_save = &H_1;
  perf_log.pop("STEP2: Linear scale and addition");

  // --- Output interpolation matrix
  perf_log.push("Output matrix");
  libMesh::PetscMatrix<libMesh::Number> sys_H_mat(*sys_mat_PETSC_H_save,WorldComm);

  // Print MatLab debugging output? Variable defined at "carl_headers.h"
  #ifdef PRINT_MATLAB_DEBUG
    sys_H_mat.print_matlab(input_params.output_base + "interpolation_mat.m");
  #endif

  carl::write_PETSC_matrix(*sys_mat_PETSC_H_save, input_params.output_base + "interpolation_mat.petscmat",0,WorldComm.get(),1);
  perf_log.pop("Output matrix");

  // - Clean up

  perf_log.push("Clean up");
  MatDestroy(&mass_invert_couplingt_A);
  MatDestroy(&mass_invert_couplingt_B);
  MatDestroy(&coupling_matrix_A);
  MatDestroy(&coupling_matrix_B);
  MatDestroy(&H_1);
  MatDestroy(&H_2);
  perf_log.pop("Clean up");

  return 0;

}