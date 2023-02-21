/*
 *
 *  Created on: Dec 1, 2022
 *      Author: Chensheng Luo
 */

#include "CArl_build_interpolation.h"
/** \file CArl_build_interpolation.cpp
 * \brief **DYN** Program used to calculate interpolation matrix by parallel calculation.

 Usage: `./CArl_build_interpolation -i [input file]`

The program CArl_build_interpolation will calculate the interpolation matrix \f$H=\beta^A (\Delta t^A)^2 C^A (M^A)^{-1} (C^A)^T+\beta^B (\Delta t^B)^2 C^B (M^B)^{-1} (C^B)^T \f$. 
This calculation is done in following steps:
  1. Calculate \f$(M^k)^{-1} (C^k)^T \f$ by first separating \f$ (C^i)^T \f$ into several columns \f$C^i_{ttt} \f$, then solving \f$ M^i X^i_{ttt}=C^i_{ttt} \f$ and finally putting \f$X^i_{ttt}\f$ together to form a matrix \f$G^i\f$
  1. Calculate \f$H=\beta^A (\Delta t^A)^2 C^A G^A+\beta^B (\Delta t^B)^2 C^B G^B \f$

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
  Vec CT_column_vec_A,CT_column_vec_B;
  Vec sol_vec_A,sol_vec_B;
  //Vec columnsA, columnsB;
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
  VecCreate(WorldComm.get(),&sol_vec_A);
  VecCreate(WorldComm.get(),&sol_vec_B);
  //VecCreate(WorldComm.get(),&columnsA);
  //VecCreate(WorldComm.get(),&columnsB);
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

  // --- Transpose matrix

  perf_log.push("Transpose matrix");

  MatGetSize(coupling_matrix_A, &CnA, &CmA);
  MatGetSize(coupling_matrix_B, &CnB, &CmB);
  MatSetSizes(transpose_coupling_matrix_A, PETSC_DECIDE, PETSC_DECIDE, CmA, CnA);
  MatSetSizes(transpose_coupling_matrix_B, PETSC_DECIDE, PETSC_DECIDE, CmB, CnB);
  VecSetSizes(CT_column_vec_A, PETSC_DECIDE, CmA);
  VecSetSizes(CT_column_vec_B, PETSC_DECIDE, CmB);
  VecSetSizes(sol_vec_A, PETSC_DECIDE, CmA);
  VecSetSizes(sol_vec_B, PETSC_DECIDE, CmB);
  //VecSetSizes(columnsA, PETSC_DECIDE, CnA);
  //VecSetSizes(columnsB, PETSC_DECIDE, CnB);
  MatTranspose(coupling_matrix_A, MAT_INITIAL_MATRIX, &transpose_coupling_matrix_A);
  MatTranspose(coupling_matrix_B, MAT_INITIAL_MATRIX, &transpose_coupling_matrix_B);

  VecSetFromOptions(sol_vec_A);
  VecSetFromOptions(sol_vec_B);
  
  perf_log.pop("Transpose matrix");

  // --- Calculate M-1C A
  perf_log.push("Inverse A");
  
  //libMesh::PetscVector<libMesh::Number> sys_sol_vec(WorldComm);
  //libMesh::PetscMatrix<libMesh::Number> sys_mat_A(mass_A,WorldComm);
  
  PetscScalar *rowA;
  PetscInt low,high;
  PetscInt local_column;//local_N
  KSP kspA;


  MatCreateAIJ(WorldComm.get(), PETSC_DECIDE, PETSC_DECIDE, CmA, CnA, CmA, NULL, CnA, NULL, &mass_invert_couplingt_A);

  KSPCreate(WorldComm.get(),&kspA);
  KSPSetFromOptions(kspA);
  KSPSetOperators(kspA,mass_A,mass_A);
  KSPSetUp(kspA);

  for (i = 0; i < CnA; i++){

    //VecSetValues(Vec x, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
    MatGetColumnVector(transpose_coupling_matrix_A, CT_column_vec_A, i);

    // libMesh::PetscVector<libMesh::Number> sys_rhs_vec_A(CT_column_vec_A,WorldComm);

    // sys_sol_vec.init(sys_rhs_vec_A);
    
    // MatGetLocalSize(mass_A,NULL,&local_N);

    // // --- Linear solver
    // libMesh::PetscLinearSolver<libMesh::Number> KSP_solver(WorldComm);
    // KSP_solver.init("sys");

    // // Solve!
    // KSP_solver.solve(sys_mat_A,sys_sol_vec,sys_rhs_vec_A,1e-9,1000);

    KSPSolve(kspA,CT_column_vec_A,sol_vec_A);

    VecGetLocalSize(sol_vec_A, &local_column);

    VecGetOwnershipRange(sol_vec_A, &low, &high);

    VecGetArray(sol_vec_A, &rowA);

    if(local_column != 0 && high == low+local_column ){
      for(j = 0; j<local_column; j++){
        MatSetValue(mass_invert_couplingt_A,j+low,i,rowA[j],INSERT_VALUES);
      }
    }
    // sys_rhs_vec_A.clear();
    // KSP_solver.clear();

  }
  MatAssemblyBegin(mass_invert_couplingt_A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mass_invert_couplingt_A, MAT_FINAL_ASSEMBLY);
  
  perf_log.pop("Inverse A");

  

  perf_log.push("Inverse B");
  //libMesh::PetscMatrix<libMesh::Number> sys_mat_B(mass_B,WorldComm);
  PetscScalar *rowB;
  KSP kspB;
  MatCreateAIJ(WorldComm.get(), PETSC_DECIDE, PETSC_DECIDE, CmB, CnB, CmB, NULL, CnB, NULL, &mass_invert_couplingt_B);

  KSPCreate(WorldComm.get(),&kspB);
  KSPSetFromOptions(kspB);
  KSPSetOperators(kspB,mass_B,mass_B);
  KSPSetUp(kspB);

  for (i = 0; i < CnB; i++){
    MatGetColumnVector(transpose_coupling_matrix_B, CT_column_vec_B, i);

    // libMesh::PetscVector<libMesh::Number> sys_rhs_vec_B(CT_column_vec_B,WorldComm);
    
    // sys_sol_vec.init(sys_rhs_vec_B);
    
    // MatGetLocalSize(mass_B,NULL,&local_N);

    // // --- Linear solver
    // libMesh::PetscLinearSolver<libMesh::Number> KSP_solver(WorldComm);
    // KSP_solver.init("sys");

    // // Solve!
    // KSP_solver.solve(sys_mat_B,sys_sol_vec,sys_rhs_vec_B,1e-9,1000);
    
    KSPSolve(kspB,CT_column_vec_B,sol_vec_B);

    VecGetLocalSize(sol_vec_B, &local_column);

    VecGetOwnershipRange(sol_vec_B, &low, &high);

    VecGetArray(sol_vec_B, &rowB);

    if(local_column != 0 && high == low+local_column ){
      for(j = 0; j<local_column; j++){
        MatSetValue(mass_invert_couplingt_B,j+low,i,rowB[j],INSERT_VALUES);
      }
    }
    // sys_rhs_vec_B.clear();
    // KSP_solver.clear();
  }
  MatAssemblyBegin(mass_invert_couplingt_B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mass_invert_couplingt_B, MAT_FINAL_ASSEMBLY);
  

  perf_log.pop("Inverse B");

  // --- Product coupling and mass matrix

  // Create result matrix
  MatCreateAIJ(WorldComm.get(),PETSC_DECIDE,PETSC_DECIDE,CnA,CnA,CnA,NULL,CnA,NULL,&H_1);
  MatCreateAIJ(WorldComm.get(),PETSC_DECIDE,PETSC_DECIDE,CnB,CnB,CnB,NULL,CnB,NULL,&H_2);

  // - Matrix product
  perf_log.push("Matrix product");
  //H_1=C_A*MCT_A
  MatMatMult(coupling_matrix_A,
        mass_invert_couplingt_A,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &H_1);

  // H_2=C_B*MCT_B
  MatMatMult(coupling_matrix_B,
        mass_invert_couplingt_B,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &H_2);
  perf_log.pop("Matrix product");
              
  // - Linear scale and addition
  perf_log.push("Linear scale and addition");
  //H_2=H_2*alpha_B
  MatScale(H_2,input_params.newmark_B.deltat*input_params.newmark_B.deltat*input_params.newmark_B.beta);
  //H_1=H_1*alpha_A+H_2
  MatAYPX(H_1,input_params.newmark_A.deltat*input_params.newmark_A.deltat*input_params.newmark_A.beta,H_2, SAME_NONZERO_PATTERN);
  perf_log.pop("Linear scale and addition");

  // --- Output interpolation matrix
  perf_log.push("Output matrix");
  Mat* sys_mat_PETSC_H_save;
  sys_mat_PETSC_H_save = &H_1;
  libMesh::PetscMatrix<libMesh::Number> sys_H_mat(*sys_mat_PETSC_H_save,WorldComm);

  // Print MatLab debugging output? Variable defined at "carl_headers.h"
  #ifdef PRINT_MATLAB_DEBUG
    sys_H_mat.print_matlab(input_params.output_base + ".m");
  #endif

  carl::write_PETSC_matrix(*sys_mat_PETSC_H_save, input_params.output_base+".petscmat",0,WorldComm.get(),1);
  perf_log.pop("Output matrix");


  // --- Clean up

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
  VecDestroy(&sol_vec_A);
  VecDestroy(&sol_vec_B);
  KSPDestroy(&kspA);
  KSPDestroy(&kspB);

  perf_log.pop("Clean up");
  
  return 0;
}