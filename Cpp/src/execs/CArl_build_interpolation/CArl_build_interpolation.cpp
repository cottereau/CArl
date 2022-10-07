/*
 *
 *  Created on: Sep 26，2021
 *      Author: Chensheng Luo
 */

#include "CArl_build_interpolation.h"
/** \file CArl_build_interpolation.cpp
 * \brief **DYN-DI** Program used to calculate interpolation matrix.

 Usage: `./CArl_build_interpolation_matrix -i [input file]`

This program will calculate the interpolation matrix \f$H=\beta^A (\Delta t^A)^2 C^A (M^A)^{-1} (C^A)^T+\beta^B (\Delta t^B)^2 C^B (M^B)^{-1} (C^B)^T \f$. 
This calculation is done in following three steps:
  1. Inversion of mass matrix (first in sequential then converted to parallel)
  2. Examination of inversion quality
  3. Other calculations

The input file is parsed by the get_input_params(GetPot& field_parser, carl_build_interpolation_params& input_params) function, and it contains the following parameters. 

 */

void load_PETSC_matrices(std::string path_dynamic_matrices, 
        std::string path_coupling_matix,
        libMesh::Parallel::Communicator& WorldComm,
        Mat& M_tilde,
        Mat& dense_M_tilde,
        Mat& Coupling);// Load mass matrix in parallel and in sequantial, and load coupling matrix

void new_invert_PETSC_matrix(libMesh::Parallel::Communicator& WorldComm,
        Mat& dense_matrix,
        Mat& inv_matrix,
        std::string inv_path);// Invert sequantially with LU factorisation and read it in parallel

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

  // - Definition of matrix

  // Mass matrix and its dense matrix(the inverse of mass matrix is done by the dense matrix)
  Mat m_tilde_A;
  Mat m_tilde_B;
  Mat dense_m_tilde_A;
  Mat dense_m_tilde_B;

  // Inversed matrix
  Mat inv_m_tilde_A;
  Mat inv_m_tilde_B;

  Mat coupling_matrix_A;
  Mat coupling_matrix_B;

  Mat H_1;
  Mat H_2;

  PetscReal exam_norm1;
  PetscReal exam_norm2;

  // - Load matrix

  perf_log.push("Load PETSc Matrices");

  load_PETSC_matrices(input_params.path_tilde_matrix_A,
                input_params.path_macro_coupling_matrix,
                WorldComm,
                m_tilde_A, 
                dense_m_tilde_A,
                coupling_matrix_A);

  load_PETSC_matrices(input_params.path_tilde_matrix_B,
                input_params.path_micro_coupling_matrix,
                WorldComm,
                m_tilde_B,
                dense_m_tilde_B,
                coupling_matrix_B);

  perf_log.pop("Load PETSc Matrices");

  // --- Compute part
  // - Invert mass matrix

  perf_log.push("Invert matrix");

  new_invert_PETSC_matrix(WorldComm,
                dense_m_tilde_B, 
                inv_m_tilde_B,
                input_params.output_base+input_params.path_inv_matrix_B);

  new_invert_PETSC_matrix(WorldComm,
                dense_m_tilde_A, 
                inv_m_tilde_A,
                input_params.output_base+input_params.path_inv_matrix_A);

  perf_log.pop("Invert matrix");

  // - Invert examination: we calculate \f$ || M^k (M^k)^{-1}-I||_{\inf} \f$

  perf_log.push("Invert examination");

  exam_norm1= exam_invert_norm(m_tilde_A, inv_m_tilde_A);
  exam_norm2= exam_invert_norm(m_tilde_B, inv_m_tilde_B);

  MatDestroy(&m_tilde_A);
  MatDestroy(&m_tilde_B);

  std::ofstream result_file;
  result_file.open(input_params.output_base + "result_test_inversion.txt");
  result_file.precision(15);
  result_file << "Test Inversion A:" << exam_norm1 << std::endl;
  result_file << "Test Inversion B:" << exam_norm2 << std::endl;
  result_file << "If one of these two values is too high, use CArl-Dyn-CG solver!" << std::endl;
  result_file.close();

  std::string error_message="Inversion norm A is bigger than the threshold! \n Please check result_test_inversion.txt in your output base to know the residual norm. \n This may be due to an extreme small residual or the mass matrix isn't invertible with  LU Factorisation. Try to choose a smaller residual or use directly CArl-Dyn-CG solver!\n";
  if (exam_norm1 > input_params.invert_threshold){
    homemade_error_msg(error_message);
  }

  if (exam_norm2 > input_params.invert_threshold){
    homemade_error_msg(error_message);
  }


  perf_log.pop("Invert examination");

  // - Product coupling and mass matrix

  // Examination of dimension
  PetscInt        petsc_m=0, petsc_n=0,petsc_a=0,petsc_b=0,petsc_c=0,petsc_d=0;
  MatGetSize(coupling_matrix_A, &petsc_m, &petsc_n);
  MatGetSize(inv_m_tilde_A, &petsc_a, &petsc_b);
  //MatGetSize(t_coupling_matrix_A, &petsc_c, &petsc_d);
  if (petsc_a != petsc_n || petsc_b != petsc_n){
    homemade_error_msg("Dimension Error! The coupling matrix and the mass matrix are not compatible for matrix product");
  }

  // Create result matrix
  MatCreateAIJ(WorldComm.get(),PETSC_DECIDE,PETSC_DECIDE,petsc_m,petsc_m,petsc_m,NULL,petsc_m,NULL,&H_1);
  MatCreateAIJ(WorldComm.get(),PETSC_DECIDE,PETSC_DECIDE,petsc_m,petsc_m,petsc_m,NULL,petsc_m,NULL,&H_2);

  // - Mass product
  perf_log.push("STEP1: Mass product");
  //H_1=C_A*M_A*C_A^T
  MatRARt(inv_m_tilde_A, 
        coupling_matrix_A,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &H_1);

  //H_2=C_B*M_B*C_B^T
  MatRARt(inv_m_tilde_B, 
        coupling_matrix_B,
        MAT_INITIAL_MATRIX,
        PETSC_DEFAULT,
        &H_2);

  perf_log.pop("STEP1: Mass product");
              
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

  perf_log.push("Clean up");
  MatDestroy(&inv_m_tilde_A);
  MatDestroy(&inv_m_tilde_B);
  MatDestroy(&coupling_matrix_A);
  MatDestroy(&coupling_matrix_B);
  MatDestroy(&H_1);
  MatDestroy(&H_2);

  perf_log.pop("Clean up");

  return 0;
}


void load_PETSC_matrices(std::string path_dynamic_matrices, 
        std::string path_coupling_matix,
        libMesh::Parallel::Communicator& WorldComm,
        Mat& M_tilde,
        Mat& dense_M_tilde,
        Mat& Coupling)
{
  // Create matrix
  MatCreate(PETSC_COMM_SELF,&dense_M_tilde);
  MatSetType(dense_M_tilde,MATSEQAIJ);
  MatCreate(WorldComm.get(),&M_tilde);
  MatCreate(WorldComm.get(),&Coupling);

  // Read Matrix
  carl::read_PETSC_matrix(Coupling, 
                path_coupling_matix, WorldComm.get());
  carl::read_PETSC_matrix(M_tilde, 
                path_dynamic_matrices, WorldComm.get());
  carl::read_PETSC_matrix(dense_M_tilde, 
                path_dynamic_matrices,PETSC_COMM_SELF);

  
}

void new_invert_PETSC_matrix(libMesh::Parallel::Communicator& WorldComm,
          Mat& dense_matrix, 
          Mat& inv_matrix,
          std::string inv_path)
{
  Mat             Ones,inv_dense_matrix;
  PetscErrorCode  ierr=0;
  MatFactorInfo   factor_info;
  IS        rperm, cperm;

  //Create identity and inverse of the dense matrix
  MatDuplicate(dense_matrix,MAT_DO_NOT_COPY_VALUES,&inv_dense_matrix);
  MatDuplicate(dense_matrix,MAT_DO_NOT_COPY_VALUES,&Ones);
  MatSetType(inv_dense_matrix, MATSEQDENSE);
  MatSeqDenseSetPreallocation(inv_dense_matrix,NULL);
  MatSetType(Ones, MATSEQDENSE);
  MatSeqDenseSetPreallocation(Ones,NULL);
  MatZeroEntries(Ones);
  MatShift(Ones,1);
  
  // Factor input matrix
  MatGetOrdering(dense_matrix,MATORDERINGNATURAL, &rperm,  &cperm);
  MatFactorInfoInitialize(&factor_info);
  MatLUFactor(dense_matrix,rperm,cperm,&factor_info);
  
  // Calculate inverse of the dense matrix 
  ierr = MatMatSolve(dense_matrix,Ones,inv_dense_matrix);

  // Reset input's factoring
  MatSetUnfactored(dense_matrix);

  //Output inverse of the dense matrix
  carl::write_PETSC_matrix(inv_dense_matrix,inv_path+".petscmat",0,PETSC_COMM_SELF,1);

  #ifdef PRINT_MATLAB_DEBUG
    carl::write_PETSC_matrix_MATLAB(inv_dense_matrix,inv_path+".m",PETSC_COMM_SELF);
  #endif

  // Cleanup
  MatDestroy(&Ones);
  MatDestroy(&dense_matrix);
  MatDestroy(&inv_dense_matrix);

  //Load Matrix in WorldComm
  MatCreate(WorldComm.get(),&inv_matrix);
  carl::read_PETSC_matrix(inv_matrix, inv_path+".petscmat", WorldComm.get());
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