#include "test_sparse_inverse.h"

int main(int argc, char *argv[])
{
	// Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);
	libMesh::Parallel::Communicator& WorldComm = init.comm();

    libMesh::PetscMatrix<libMesh::Number> matrix_A(WorldComm);
    matrix_A.init(4,4,4,4);

    matrix_A.set(0,0,1.);
//    matrix_A.set(0,1,2.);
//    matrix_A.set(0,2,3.);
//    matrix_A.set(0,3,4.);

//    matrix_A.set(1,0,2.);
    matrix_A.set(1,1,5.);
//    matrix_A.set(1,2,3.);
//    matrix_A.set(1,3,7.);

//    matrix_A.set(2,0,3.);
//    matrix_A.set(2,1,3.);
    matrix_A.set(2,2,9.);
//    matrix_A.set(2,3,6.);

//    matrix_A.set(3,0,4.);
//    matrix_A.set(3,1,7.);
//    matrix_A.set(3,2,6.);
    matrix_A.set(3,3,1.);

    matrix_A.close();

    Mat dummy_inv_A;
    MatCreate(PETSC_COMM_WORLD,&dummy_inv_A);
    MatSetType(dummy_inv_A,MATMPIAIJ);
    MatSetSizes(dummy_inv_A,PETSC_DECIDE,PETSC_DECIDE,4,4);
    MatMPIAIJSetPreallocation(dummy_inv_A,2,NULL,0,NULL);
    MatSetUp(dummy_inv_A);

    // Dummy matrices
//    Mat dummy_A, dummy_inv_A;
//
//
    libMesh::PetscVector<libMesh::Number> vector_unity(WorldComm,4,4);
    libMesh::PetscVector<libMesh::Number> vector_dummy_answer(WorldComm,4,4);

	VecSet(vector_unity.vec(),1);
	vector_unity.close();

	VecSet(vector_dummy_answer.vec(),0);
	vector_dummy_answer.close();

    // Solver
//	libMesh::PetscLinearSolver<libMesh::Number> KSP_dummy_solver(WorldComm);
//	KSP_dummy_solver.init(&matrix_A);

//	KSPSetOperators(KSP_dummy_solver.ksp(),matrix_A.mat(),NULL);

	KSP ksp;
	PC pc;

	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators(ksp, matrix_A.mat(), matrix_A.mat());
	KSPGetPC(ksp,&pc);
	PCSetFromOptions(pc);
	PCType dummy_type;
	PCGetType(pc,&dummy_type);
	std::cout << std::endl << dummy_type << std::endl << std::endl;
//	PCSetType(pc,PCSPAI);
//	PCHYPRESetType(pc,"parasails");
	KSPSetUp(ksp);
	KSPSolve(ksp,vector_unity.vec(),vector_dummy_answer.vec());
	PCComputeExplicitOperator(pc,&dummy_inv_A);
//	KSPGetOperators(KSP_dummy_solver.ksp(),&dummy_A,&dummy_inv_A);

	libMesh::PetscMatrix<libMesh::Number> matrix_invA(dummy_inv_A,WorldComm);
	matrix_invA.close();
//
//	//    KSP_dummy_solver.solve(matrix_A,vector_dummy_answer,vector_unity,1E-5,10000);
//
//	vector_dummy_answer.print_matlab();
//

//	libMesh::PetscMatrix<libMesh::Number> product_mat(WorldComm);
	matrix_A.print_matlab();
	matrix_invA.print_matlab();
	vector_dummy_answer.print_matlab();
	return 0;
}
