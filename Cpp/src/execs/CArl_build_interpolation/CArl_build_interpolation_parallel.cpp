/*
 *
 *  Created on: Dec 1, 2022
 *      Author: Chensheng Luo
 */

#include "CArl_build_interpolation.h"
/** \file CArl_build_interpolation_parallel_begin.cpp
 * \brief **DYN-DI** Program used to calculate interpolation matrix by parallel calculation.

 Usage: `./CArl_build_interpolation_matrix_parallel_begin -i [input file]`

The program CArl_build_interpolation_parallel_begin and CArl_build_interpolation_parallel_end will calculate the interpolation matrix \f$H=\beta^A (\Delta t^A)^2 C^A (M^A)^{-1} (C^A)^T+\beta^B (\Delta t^B)^2 C^B (M^B)^{-1} (C^B)^T \f$. 
This calculation is done in following three steps:
  1. Separating \f$ (C^i)^T \f$ into several columns: C^i_{ttt} (begin)
  2. Inverting by solving \f$ M^i X^i_{ttt}=C^i_{ttt} \f$ (solving)
  3. By putting \f$X^i_{ttt}\f$ together to form a matrix \f$G^i\f$ and then calculate \f$H=\beta^A (\Delta t^A)^2 C^A G^A+\beta^B (\Delta t^B)^2 C^B G^B \f (end)$

The input file is parsed by the get_input_params(GetPot& field_parser, carl_build_interpolation_params& input_params) function. 

 */
PetscReal exam_invert_norm(Mat& matrix, Mat& inv_matrix);// Exam the inversion quality by calculating \f$ || M^k (M^k)^{-1}-I||_{\inf} \f$ and by comparing with threshold

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

  if(input_params.m_bParallelPossible==false){
    homemade_error_msg("[CArl Error] Parallel Not Possible, check your parameter!");
  }

  // - Definition of matrix

  Mat coupling_matrix_A;
  Mat coupling_matrix_B;
  Mat transpose_coupling_matrix_A;
  Mat transpose_coupling_matrix_B;
  Mat mass_A;
  Mat mass_B;
  Mat mass_invert_couplingt_A,mass_invert_couplingt_B;
  Mat H_1;
  Mat H_2;
  Vec CT_column_vec_A;
  Vec CT_column_vec_B;
  PetscInt CmA,CnA,CmB,CnB;
  PetscInt i,j;
  // - Load matrix

  perf_log.push("Load PETSc Matrices");

  MatCreate(WorldComm.get(),&coupling_matrix_A);
  MatCreate(WorldComm.get(),&coupling_matrix_B);
  MatCreate(WorldComm.get(),&transpose_coupling_matrix_A);
  MatCreate(WorldComm.get(),&transpose_coupling_matrix_B);
  MatCreate(WorldComm.get(),&mass_A);
  MatCreate(WorldComm.get(),&mass_B);
  VecCreate(WorldComm.get(),&CT_column_vec_A);
  VecCreate(WorldComm.get(),&CT_column_vec_B);
  VecSetType(CT_column_vec_A, VECMPI);
  VecSetType(CT_column_vec_B, VECMPI);

  carl::read_PETSC_matrix(coupling_matrix_A, 
                input_params.path_macro_coupling_matrix, WorldComm.get());

  carl::read_PETSC_matrix(coupling_matrix_B, 
                input_params.path_micro_coupling_matrix, WorldComm.get());

  carl::read_PETSC_matrix(mass_A,
                input_params.path_tilde_matrix_A, WorldComm.get());

  carl::read_PETSC_matrix(mass_B,
                input_params.path_tilde_matrix_B, WorldComm.get());

  perf_log.pop("Load PETSc Matrices");

  // - Transpose matrix

  perf_log.push("Transpose matrix");

  MatGetSize(coupling_matrix_A, &CnA, &CmA);
  MatGetSize(coupling_matrix_B, &CnB, &CmB);
  MatSetSizes(transpose_coupling_matrix_A, PETSC_DECIDE, PETSC_DECIDE, CmA, CnA);
  MatSetSizes(transpose_coupling_matrix_B, PETSC_DECIDE, PETSC_DECIDE, CmB, CnB);
  VecSetSizes(CT_column_vec_A, PETSC_DECIDE, CmA);
  VecSetSizes(CT_column_vec_B, PETSC_DECIDE, CmB);
  MatTranspose(coupling_matrix_A, MAT_INITIAL_MATRIX, &transpose_coupling_matrix_A);
  MatTranspose(coupling_matrix_B, MAT_INITIAL_MATRIX, &transpose_coupling_matrix_B);
  
  perf_log.pop("Transpose matrix");

  // - Calculate M-1C A
  PetscInt    local_N,local_column;
  libMesh::PetscVector<libMesh::Number> sys_sol_vec(WorldComm);

  perf_log.push("Inverse A");
  libMesh::PetscMatrix<libMesh::Number> sys_mat_A(mass_A,WorldComm);
  
  PetscScalar *rowA;
  PetscInt n,low,high;
  MatCreateAIJ(WorldComm.get(), PETSC_DECIDE, PETSC_DECIDE, CmA, CnA, CmA, NULL, CnA, NULL, &mass_invert_couplingt_A);

  for (i = 0; i < CnA; i++){
    MatGetColumnVector(transpose_coupling_matrix_A, CT_column_vec_A, i);

    libMesh::PetscVector<libMesh::Number> sys_rhs_vec_A(CT_column_vec_A,WorldComm);

    sys_sol_vec.init(sys_rhs_vec_A);
    
    MatGetLocalSize(mass_A,NULL,&local_N);

    // --- Linear solver
    libMesh::PetscLinearSolver<libMesh::Number> KSP_solver(WorldComm);
    KSP_solver.init("sys");

    // Solve!
    KSP_solver.solve(sys_mat_A,sys_sol_vec,sys_rhs_vec_A,1e-9,1000);

    Vec local_vector;

    VecGetLocalSize(sys_sol_vec.vec(), &local_column);

    VecGetOwnershipRange(sys_sol_vec.vec(), &low, &high);

    VecGetArray(sys_sol_vec.vec(), &rowA);

    std::cout << "Local column :" << local_column << std::endl;

    std::cout << "Low :"<< low<< "High:"<< high << std::endl;

    if(local_column != 0 && high == low+local_column ){
      for(j = 0; j<local_column; j++){
        MatSetValue(mass_invert_couplingt_A,j+low,i,rowA[j],INSERT_VALUES);
      }
    }

  }
  MatAssemblyBegin(mass_invert_couplingt_A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mass_invert_couplingt_A, MAT_FINAL_ASSEMBLY);

  perf_log.pop("Inverse A");

  

  perf_log.push("Inverse B");
  libMesh::PetscMatrix<libMesh::Number> sys_mat_B(mass_B,WorldComm);
  PetscScalar *rowB;
  MatCreateAIJ(WorldComm.get(), PETSC_DECIDE, PETSC_DECIDE, CmB, CnB, CmB, NULL, CnB, NULL, &mass_invert_couplingt_B);


  for (i = 0; i < CnB; i++){
    MatGetColumnVector(transpose_coupling_matrix_B, CT_column_vec_B, i);

    libMesh::PetscVector<libMesh::Number> sys_rhs_vec_B(CT_column_vec_B,WorldComm);
    
    sys_sol_vec.init(sys_rhs_vec_B);
    
    MatGetLocalSize(mass_B,NULL,&local_N);

    // --- Linear solver
    libMesh::PetscLinearSolver<libMesh::Number> KSP_solver(WorldComm);
    KSP_solver.init("sys");

    // Solve!
    KSP_solver.solve(sys_mat_B,sys_sol_vec,sys_rhs_vec_B,1e-9,1000);

    Vec local_vector;

    VecGetLocalSize(sys_sol_vec.vec(), &local_column);

    VecGetOwnershipRange(sys_sol_vec.vec(), &low, &high);

    VecGetArray(sys_sol_vec.vec(), &rowB);

    std::cout << "Local column :" << local_column << std::endl;

    std::cout << "Low :"<< low<< "High:"<< high << std::endl;

    if(local_column != 0 && high == low+local_column ){
      for(j = 0; j<local_column; j++){
        MatSetValue(mass_invert_couplingt_B,j+low,i,rowB[j],INSERT_VALUES);
      }
    }
  }
  MatAssemblyBegin(mass_invert_couplingt_B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mass_invert_couplingt_B, MAT_FINAL_ASSEMBLY);

  perf_log.pop("Inverse B");

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
  MatDestroy(&coupling_matrix_A);
  MatDestroy(&coupling_matrix_B);
  MatDestroy(&mass_A);
  MatDestroy(&mass_B);
  MatDestroy(&mass_invert_couplingt_A);
  MatDestroy(&mass_invert_couplingt_B);
  MatDestroy(&transpose_coupling_matrix_A);
  MatDestroy(&transpose_coupling_matrix_B);
  MatDestroy(&H_1);
  MatDestroy(&H_2);
  VecDestroy(&CT_column_vec_A);
  VecDestroy(&CT_column_vec_B);

  perf_log.pop("Clean up");
  
  return 0;
}

PetscReal exam_invert_norm(Mat& matrix, Mat& inv_matrix){
  PetscReal exam_norm;
  Mat MinvM;

  MatProductCreate(matrix, inv_matrix, NULL,&MinvM);
  MatProductSetType(MinvM,MATPRODUCT_AB);
  MatProductSetFromOptions(MinvM);
  MatProductSymbolic(MinvM);
  MatProductNumeric(MinvM);
  MatShift(MinvM,-1);
  MatNorm(MinvM,NORM_INFINITY,&exam_norm);

  MatDestroy(&MinvM);
  return exam_norm;
}