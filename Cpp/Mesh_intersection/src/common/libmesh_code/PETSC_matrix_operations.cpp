#include "PETSC_matrix_operations.h"

void carl::lump_matrix(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
							libMesh::PetscMatrix<libMesh::Number>& matrixOutput)
{
	if(matrixOutput.initialized())
	{
		matrixOutput.clear();
	}

	int M = matrixInput.m();
	int N = matrixInput.n();

	PetscInt local_M, local_N;

	MatGetLocalSize(matrixInput.mat(),&local_M,&local_N);

	// It will be a diagonal matrix, so no need of a heavy preallocation
	matrixOutput.init(M,N,local_M,local_N,1,0);

	libMesh::PetscVector<libMesh::Number> UnityVec(matrixInput.comm(),M,local_M);
	libMesh::PetscVector<libMesh::Number> DummyVector(matrixInput.comm(),M,local_M);

	VecSet(UnityVec.vec(),1);

	UnityVec.close();

	matrixInput.vector_mult(DummyVector,UnityVec);

	MatDiagonalSet(matrixOutput.mat(),DummyVector.vec(),INSERT_VALUES);

	matrixOutput.close();
}

void carl::lump_matrix_and_invert(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
										libMesh::PetscMatrix<libMesh::Number>& matrixOutput)
{
	if(matrixOutput.initialized())
	{
		matrixOutput.clear();
	}

	int M = matrixInput.m();
	int N = matrixInput.n();

	PetscInt local_M, local_N;

	MatGetLocalSize(matrixInput.mat(),&local_M,&local_N);

	// It will be a diagonal matrix, so no need of a heavy preallocation
	matrixOutput.init(M,N,local_M,local_N,1,0);

	libMesh::PetscVector<libMesh::Number> UnityVec(matrixInput.comm(),M,local_M);
	libMesh::PetscVector<libMesh::Number> DummyVector(matrixInput.comm(),M,local_M);

	VecSet(UnityVec.vec(),1);

	UnityVec.close();

	matrixInput.vector_mult(DummyVector,UnityVec);

	DummyVector.reciprocal();

	MatDiagonalSet(matrixOutput.mat(),DummyVector.vec(),INSERT_VALUES);

	matrixOutput.close();
}

void carl::lump_matrix_and_invert(		libMesh::PetscMatrix<libMesh::Number>& matrixInput,
										libMesh::PetscVector<libMesh::Number>& vecOutput)
{
	if(vecOutput.initialized())
	{
		vecOutput.clear();
	}

	int M = matrixInput.m();
	int N = matrixInput.n();

	PetscInt local_M, local_N;

	MatGetLocalSize(matrixInput.mat(),&local_M,&local_N);

	// It will be a diagonal matrix, so no need of a heavy preallocation
	vecOutput.init(M,local_M,false);

	libMesh::PetscVector<libMesh::Number> UnityVec(matrixInput.comm(),M,local_M);

	VecSet(UnityVec.vec(),1);

	UnityVec.close();

	matrixInput.vector_mult(vecOutput,UnityVec);

	vecOutput.reciprocal();
}


void carl::print_matrix(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix)
{
	libMesh::Real accumulator = 0;
	std::cout << "| M_i,j : " << CouplingTestMatrix.m() << " x " << CouplingTestMatrix.n() << std::endl;
	for(unsigned int iii = 0; iii < CouplingTestMatrix.m(); ++iii)
	{
		for(unsigned int jjj = 0; jjj < CouplingTestMatrix.n(); ++jjj)
		{
			accumulator += CouplingTestMatrix(iii,jjj);
			if(CouplingTestMatrix.m() < 101 && CouplingTestMatrix.n() < 101)
			{
				std::cout << " " << CouplingTestMatrix(iii,jjj);
			}
		}
		if(CouplingTestMatrix.m() < 101 && CouplingTestMatrix.n() < 101)
		{
			std::cout << std::endl;
		}
	}
	std::cout << "|" << std::endl;
	std::cout << "| Sum( M_i,j ) = " << accumulator << std::endl << std::endl;
}

void carl::print_matrix_dim(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix)
{
	std::cout << "| M_i,j : " << CouplingTestMatrix.m() << " x " << CouplingTestMatrix.n() << std::endl  << std::endl;
}

void carl::solve_linear_PETSC(	libMesh::PetscMatrix<libMesh::Number>& A,
								libMesh::PetscVector<libMesh::Number>& b,
								libMesh::PetscVector<libMesh::Number>& x,
								KSP& ksp, PC& pc)
{
	/*
	 * 			Solve the system A*x = b using PETSc's linear solver, using a
	 * 		Krilov method, with the linear solver context "ksp" and
	 * 		preconditioner "pc".
	 */

	// Set up inital variables
}

