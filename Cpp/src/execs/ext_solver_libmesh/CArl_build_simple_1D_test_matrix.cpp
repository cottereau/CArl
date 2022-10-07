/*
 * \file CArl_build_simple_1D_test_matrix.cpp
 *
 *  Created on: Aug. 13, 2022
 *      Author: Chensheng Luo
 *
 * \brief  **TEST** This is only a simple 1D case used to test the code performance
 * It corresponds two beams coupled together, with a beam bounded and a traction force applied on another beam.
 * 
 *  This program's input file description can be found at the documentation of the function carl::get_input_params(GetPot& field_parser, libmesh_apply_solution_dyn_input_params& input_params). 
 */
#include "CArl_build_simple_1D_test_matrix.h"

int main(int argc, char** argv) {

  // --- Initialize libMesh
  libMesh::LibMeshInit init(argc, argv);

  // Do performance log?
  libMesh::PerfLog perf_log("Test program");

  // libMesh's C++ / MPI communicator wrapper
  libMesh::Parallel::Communicator& WorldComm = init.comm();

  // Number of processors and processor rank.
  int rank = WorldComm.rank();
  int nodes = WorldComm.size();

  // --- Set up inputs

  // Command line parser
  GetPot command_line(argc, argv);

  // File parser
  GetPot field_parser;

  // If there is an input file, parse it to get the parameters. Else, parse the command line
  std::string input_filename;
  if (command_line.search(2, "--inputfile", "-i")) {
    input_filename = command_line.next(input_filename);
    field_parser.parse_input_file(input_filename, "#", "\n", " \t\n");
  } else {
    field_parser = command_line;
  }

  carl::simple_1D_test_matrix_input_params input_params;
  carl::get_input_params(field_parser, input_params);

  // --- Set the matrix and vectors

  Mat MA;
  Mat MB;
  Mat DA;
  Mat DB;
  Mat KA;
  Mat KB;
  Mat CA;
  Mat CB;
  Vec FA;
  Vec FB;
  PetscInt N=input_params.N;
  PetscInt Ne=input_params.N*input_params.e/input_params.L;
  int i;

  MatCreate(WorldComm.get(),&MA);
  MatCreate(WorldComm.get(),&MB);
  MatCreate(WorldComm.get(),&DA);
  MatCreate(WorldComm.get(),&DB);
  MatCreate(WorldComm.get(),&KA);
  MatCreate(WorldComm.get(),&KB);
  MatCreate(WorldComm.get(),&CA);
  MatCreate(WorldComm.get(),&CB);
  VecCreate(WorldComm.get(),&FA);
  VecCreate(WorldComm.get(),&FB);

  
  MatSetSizes(MA,PETSC_DECIDE,PETSC_DECIDE,N,N);
  MatSetSizes(KA,PETSC_DECIDE,PETSC_DECIDE,N,N);
  MatSetSizes(DA,PETSC_DECIDE,PETSC_DECIDE,N,N);

  MatSetSizes(MB,PETSC_DECIDE,PETSC_DECIDE,N+1,N+1);
  MatSetSizes(KB,PETSC_DECIDE,PETSC_DECIDE,N+1,N+1);
  MatSetSizes(DB,PETSC_DECIDE,PETSC_DECIDE,N+1,N+1);

  MatSetSizes(CA,PETSC_DECIDE,PETSC_DECIDE,Ne+1,N);
  MatSetSizes(CB,PETSC_DECIDE,PETSC_DECIDE,Ne+1,N+1);

  VecSetSizes(FA,PETSC_DECIDE,N);
  VecSetSizes(FB,PETSC_DECIDE,N+1);
  VecSetFromOptions(FA);
  VecSetFromOptions(FB);

  MatSetType(MA,MATAIJ);
  MatSetType(MB,MATAIJ);
  MatSetType(CA,MATAIJ);
  MatSetType(CB,MATAIJ);
  MatSetType(DA,MATAIJ);
  MatSetType(DB,MATAIJ);
  MatSetType(KA,MATAIJ);
  MatSetType(KB,MATAIJ);


  MatSetUp(MA);
  MatSetUp(MB);
  MatSetUp(KA);
  MatSetUp(KB);
  MatSetUp(CA);
  MatSetUp(CB);
  MatSetUp(DA);
  MatSetUp(DB);
  

  // --- Assembly the matrix and vectors

  PetscScalar RhoLOverN= input_params.rho*input_params.L/N;

  PetscScalar massBeginA[2]={2.0/3*RhoLOverN,1.0/6*RhoLOverN};
  PetscScalar massIntersectA[3]={1.0/6*RhoLOverN,1.0/2*RhoLOverN,1.0/12*RhoLOverN};
  PetscScalar massEndA[2]={1.0/12*RhoLOverN,1.0/6*RhoLOverN};

  PetscScalar massBeginB[2]={1.0/6*RhoLOverN,1.0/12*RhoLOverN};
  PetscScalar massIntersectB[3]={1.0/12*RhoLOverN,1.0/2*RhoLOverN,1.0/6*RhoLOverN};
  PetscScalar massEndB[2]={1.0/6*RhoLOverN,1.0/3*RhoLOverN};

  PetscScalar massMiddle[3]={1.0/6*RhoLOverN,2.0/3*RhoLOverN,1.0/6*RhoLOverN};
  PetscScalar massCoupling[3]={1.0/12*RhoLOverN,1.0/3*RhoLOverN,1.0/12*RhoLOverN};
  

  PetscScalar ENOverL= input_params.E*input_params.N/input_params.L;

  PetscScalar stiffnessBeginA[2]={2.0*ENOverL,-1.0*ENOverL};
  PetscScalar stiffnessIntersectA[3]={-1.0*ENOverL,3.0/2*ENOverL,-1.0/2*ENOverL};
  PetscScalar stiffnessEndA[2]={-1.0/2*ENOverL,1.0/2*ENOverL};

  PetscScalar stiffnessBeginB[2]={1.0/2*ENOverL,-1.0/2*ENOverL};
  PetscScalar stiffnessIntersectB[3]={-1.0/2*ENOverL,3.0/2*ENOverL,-1.0*ENOverL};
  PetscScalar stiffnessEndB[2]={-1.0*ENOverL,1.0*ENOverL};

  PetscScalar stiffnessMiddle[3]={-1.0*ENOverL,2.0*ENOverL,-1.0*ENOverL};
  PetscScalar stiffnessCoupling[3]={-1.0/2*ENOverL,1.0*ENOverL,-1.0/2*ENOverL};
  

  PetscInt    cols2[]= {0,1};
  PetscInt    cols1[]={0};
  MatSetValues(MA,1,cols1,2,cols2,massBeginA,INSERT_VALUES);
  MatSetValues(MB,1,cols1,2,cols2,massBeginB,INSERT_VALUES);
  MatSetValues(KA,1,cols1,2,cols2,stiffnessBeginA,INSERT_VALUES);
  MatSetValues(KB,1,cols1,2,cols2,stiffnessBeginB,INSERT_VALUES);

  PetscInt    cols3[] = {0,1,2};
  for(i=1;i<=(N-2);i++){
    cols3[0]=i-1;
    cols3[1]=i;
    cols3[2]=i+1;
    cols1[0]=i;
    if(i<=(N-Ne-2)){
        MatSetValues(MA,1,cols1,3,cols3,massMiddle,INSERT_VALUES);
        MatSetValues(KA,1,cols1,3,cols3,stiffnessMiddle,INSERT_VALUES);
    }else if(i==(N-Ne-1)){
        MatSetValues(MA,1,cols1,3,cols3,massIntersectA,INSERT_VALUES);
        MatSetValues(KA,1,cols1,3,cols3,stiffnessIntersectA,INSERT_VALUES);      
    }else{
        MatSetValues(MA,1,cols1,3,cols3,massCoupling,INSERT_VALUES);
        MatSetValues(KA,1,cols1,3,cols3,stiffnessCoupling,INSERT_VALUES);  
    }
  }

  for(i=1;i<=(N-1);i++){
    cols3[0]=i-1;
    cols3[1]=i;
    cols3[2]=i+1;
    cols1[0]=i;
    if(i<=(Ne-1)){
        MatSetValues(MB,1,cols1,3,cols3,massCoupling,INSERT_VALUES);
        MatSetValues(KB,1,cols1,3,cols3,stiffnessCoupling,INSERT_VALUES);
    }else if(i==Ne){
        MatSetValues(MB,1,cols1,3,cols3,massIntersectB,INSERT_VALUES);
        MatSetValues(KB,1,cols1,3,cols3,stiffnessIntersectB,INSERT_VALUES);      
    }else{
        MatSetValues(MB,1,cols1,3,cols3,massMiddle,INSERT_VALUES);
        MatSetValues(KB,1,cols1,3,cols3,stiffnessMiddle,INSERT_VALUES);  
    }
  }

  cols2[0]=N-2;
  cols2[1]=N-1;
  cols1[0]=N-1;
  MatSetValues(MA,1,cols1,2,cols2,massEndA,INSERT_VALUES);
  MatSetValues(KA,1,cols1,2,cols2,stiffnessEndA,INSERT_VALUES);

  cols2[0]=N-1;
  cols2[1]=N;
  cols1[0]=N;
  MatSetValues(MB,1,cols1,2,cols2,massEndB,INSERT_VALUES);
  MatSetValues(KB,1,cols1,2,cols2,stiffnessEndB,INSERT_VALUES);



  PetscScalar couplingBegin[2]={input_params.kappa*((N/input_params.L)+1.0/(input_params.es*input_params.es)*input_params.L/(3.0*N)),
    input_params.kappa*(-(N/input_params.L)+1.0/(input_params.es*input_params.es)*input_params.L/(6.0*N))};
  PetscScalar couplingMiddle[3]={input_params.kappa*(-(N/input_params.L)+1.0/(input_params.es*input_params.es)*input_params.L/(6.0*N)),
    input_params.kappa*(2.0*(N/input_params.L)+1.0/(input_params.es*input_params.es)*2.0*input_params.L/(3.0*N)),
    input_params.kappa*(-(N/input_params.L)+1.0/(input_params.es*input_params.es)*input_params.L/(6.0*N))};
  PetscScalar couplingEndA[2]={input_params.kappa*(-(N/input_params.L)+1.0/(input_params.es*input_params.es)*input_params.L/(6.0*N)),
    input_params.kappa*((N/input_params.L)+1.0/(input_params.es*input_params.es)*input_params.L/(3.0*N))}; 

  cols2[0]=N-Ne-1;
  cols2[1]=N-Ne;
  cols1[0]=0;
  MatSetValues(CA,1,cols1,2,cols2,couplingBegin,INSERT_VALUES);
  for(i=1;i<Ne;i++){
    cols3[0]=i+N-Ne-2;
    cols3[1]=i+N-Ne-1;
    cols3[2]=i+N-Ne;
    cols1[0]=i;
    MatSetValues(CA,1,cols1,3,cols3,couplingMiddle,INSERT_VALUES);
  }
  cols2[0]=N-2;
  cols2[1]=N-1;
  cols1[0]=Ne;
  MatSetValues(CA,1,cols1,2,cols2,couplingEndA,INSERT_VALUES);

  cols2[0]=0;
  cols2[1]=1;
  cols1[0]=0;
  MatSetValues(CB,1,cols1,2,cols2,couplingBegin,INSERT_VALUES);
  for(i=1;i<=Ne;i++){
    cols3[0]=i-1;
    cols3[1]=i;
    cols3[2]=i+1;
    cols1[0]=i;
    MatSetValues(CB,1,cols1,3,cols3,couplingMiddle,INSERT_VALUES);
  }

  VecZeroEntries(FA);
  VecZeroEntries(FB);
  cols1[0]=N;
  PetscScalar force[]={10.0};
  VecSetValues(FB,1,cols1,force,INSERT_VALUES);

  MatAssemblyBegin(MA,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(MB,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(KA,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(KB,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(CA,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(CB,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(DA,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(DB,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(FA);
  VecAssemblyBegin(FB);

  MatAssemblyEnd(MA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(MB,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(KA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(KB,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(CA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(CB,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(DA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(DB,MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(FA);
  VecAssemblyEnd(FB);

  MatCopy(MA,DA,DIFFERENT_NONZERO_PATTERN);
  MatCopy(MB,DB,DIFFERENT_NONZERO_PATTERN);
  MatScale(DA,input_params.CM);
  MatScale(DB,input_params.CM);
  MatAXPY(DA,input_params.CK,KA,DIFFERENT_NONZERO_PATTERN);
  MatAXPY(DB,input_params.CK,KB,DIFFERENT_NONZERO_PATTERN);

  MatAXPY(MA,input_params.newmark_A.gamma*input_params.newmark_A.deltat,DA,DIFFERENT_NONZERO_PATTERN);
  MatAXPY(MB,input_params.newmark_B.gamma*input_params.newmark_B.deltat,DB,DIFFERENT_NONZERO_PATTERN);
  MatAXPY(MA,input_params.newmark_A.beta*input_params.newmark_A.deltat*input_params.newmark_A.deltat,KA,SAME_NONZERO_PATTERN);
  MatAXPY(MB,input_params.newmark_B.beta*input_params.newmark_B.deltat*input_params.newmark_B.deltat,KB,SAME_NONZERO_PATTERN);

  //---Output

  carl::write_PETSC_matrix(MA, input_params.output_baseA + "_sys_mat.petscmat",0,WorldComm.get(),1);
  carl::write_PETSC_matrix(KA, input_params.output_baseA + "_K_mat.petscmat",0,WorldComm.get(),1);
  carl::write_PETSC_matrix(DA, input_params.output_baseA + "_damping.petscmat",0,WorldComm.get(),1);
  carl::write_PETSC_vector(FA, input_params.output_baseA + "_sys_rhs_vec.petscvec",0,WorldComm.get(),1);
  carl::write_PETSC_matrix(CA, input_params.coupling_base + "/coupling_matrix_macro.petscmat",0,WorldComm.get(),1);

  carl::write_PETSC_matrix(MB, input_params.output_baseB + "_sys_mat.petscmat",0,WorldComm.get(),1);
  carl::write_PETSC_matrix(KB, input_params.output_baseB + "_K_mat.petscmat",0,WorldComm.get(),1);
  carl::write_PETSC_matrix(DB, input_params.output_baseB + "_damping.petscmat",0,WorldComm.get(),1);
  carl::write_PETSC_vector(FB, input_params.output_baseB + "_sys_rhs_vec.petscvec",0,WorldComm.get(),1);
  carl::write_PETSC_matrix(CB, input_params.coupling_base + "/coupling_matrix_micro.petscmat",0,WorldComm.get(),1);

  libMesh::PetscMatrix<libMesh::Number> MA_mat(MA,WorldComm);
  libMesh::PetscMatrix<libMesh::Number> KA_mat(KA,WorldComm);
  libMesh::PetscMatrix<libMesh::Number> DA_mat(DA,WorldComm);
  libMesh::PetscVector<libMesh::Number> FA_vec(FA,WorldComm);
  libMesh::PetscMatrix<libMesh::Number> CA_mat(CA,WorldComm);

  libMesh::PetscMatrix<libMesh::Number> MB_mat(MB,WorldComm);
  libMesh::PetscMatrix<libMesh::Number> KB_mat(KB,WorldComm);
  libMesh::PetscMatrix<libMesh::Number> DB_mat(DB,WorldComm);
  libMesh::PetscVector<libMesh::Number> FB_vec(FB,WorldComm);
  libMesh::PetscMatrix<libMesh::Number> CB_mat(CB,WorldComm);

// Print MatLab debugging output? Variable defined at "carl_headers.h"
    #ifdef PRINT_MATLAB_DEBUG
        MA_mat.print_matlab(input_params.output_baseA + "_sys_mat.m");
        KA_mat.print_matlab(input_params.output_baseA + "_K_mat.m");
        DA_mat.print_matlab(input_params.output_baseA + "_damping.m");
        FA_vec.print_matlab(input_params.output_baseA + "_sys_rhs_vec.m");
        CA_mat.print_matlab(input_params.coupling_base + "/coupling_matrix_macro.m");

        MB_mat.print_matlab(input_params.output_baseB + "_sys_mat.m");
        KB_mat.print_matlab(input_params.output_baseB + "_K_mat.m");
        DB_mat.print_matlab(input_params.output_baseB + "_damping.m");
        FB_vec.print_matlab(input_params.output_baseB + "_sys_rhs_vec.m");
        CB_mat.print_matlab(input_params.coupling_base + "/coupling_matrix_micro.m");
    #endif

  // --- Cleanup!
  MatDestroy(&MA);
  MatDestroy(&KA);
  MatDestroy(&DA);
  VecDestroy(&FA);
  MatDestroy(&CA);

  MatDestroy(&MB);
  MatDestroy(&KB);
  MatDestroy(&DB);
  VecDestroy(&FB);
  MatDestroy(&CB);


  return 0;
}