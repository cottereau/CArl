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

	homemade_assert_msg(M == N, "Lumping: the matrix must be a square matrix");

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

	homemade_assert_msg(M == N, "Lumping: the matrix must be a square matrix");

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

	homemade_assert_msg(M == N, "Lumping: the matrix must be a square matrix");

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
	const libMesh::Parallel::Communicator& MatrixComm =  CouplingTestMatrix.comm();

	std::cout << "| M_i,j : " << CouplingTestMatrix.m() << " x " << CouplingTestMatrix.n() << std::endl;

	int nodes = MatrixComm.size();
	PetscInt local_M, local_N;
	MatGetLocalSize(CouplingTestMatrix.mat(),&local_M,&local_N);

	bool bPrintOnTerminal = CouplingTestMatrix.m() < 101 && CouplingTestMatrix.n() < 101 && nodes > 1;
	if(bPrintOnTerminal)
	{
		for(unsigned int iii = 0; iii < CouplingTestMatrix.m(); ++iii)
		{
			for(unsigned int jjj = 0; jjj < CouplingTestMatrix.n(); ++jjj)
			{
				std::cout << " " << CouplingTestMatrix(iii,jjj);
			}
			std::cout << std::endl;
		}
	}

	libMesh::PetscVector<libMesh::Number> dummy_vec(MatrixComm,CouplingTestMatrix.n(),local_N);
	MatGetRowSum(CouplingTestMatrix.mat(),dummy_vec.vec());
	VecSum(dummy_vec.vec(),&accumulator);
	std::cout << "|" << std::endl;
	std::cout << "| Sum( M_i,j ) = " << accumulator << std::endl << std::endl;
}

void carl::print_matrix_col_line_sum(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix, const std::string name_base)
{
	libMesh::PetscVector<libMesh::Number> col_sum(CouplingTestMatrix.comm(),CouplingTestMatrix.m());
	libMesh::PetscVector<libMesh::Number> row_sum(CouplingTestMatrix.comm(),CouplingTestMatrix.n());

	for(unsigned int iii = 0; iii < CouplingTestMatrix.m(); ++iii)
	{
		for(unsigned int jjj = 0; jjj < CouplingTestMatrix.n(); ++jjj)
		{
			col_sum.add(iii,CouplingTestMatrix(iii,jjj));
			row_sum.add(jjj,CouplingTestMatrix(iii,jjj));
		}
	}

	col_sum.print_matlab(name_base + "_col.m");
	row_sum.print_matlab(name_base + "_row.m");
}

void carl::print_matrix_matlab(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix, const std::string name_base)
{
	std::cout << "| M_i,j : " << CouplingTestMatrix.m() << " x " << CouplingTestMatrix.n() << std::endl;

	CouplingTestMatrix.print_matlab(name_base);
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

void carl::check_coupling_matrix(	libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix,
									libMesh::Mesh& IntersectionMesh,
									libMesh::Real CouplingScale,
									const std::string matrixType,
									int n_var)
{
	std::cout << "| " << matrixType << std::endl;
	libMesh::Real accumulator = 0;

	std::cout << "| M_i,j : " << CouplingTestMatrix.m() << " x " << CouplingTestMatrix.n() << std::endl;
	for(unsigned int iii = 0; iii < CouplingTestMatrix.m(); ++iii)
	{
		for(unsigned int jjj = 0; jjj < CouplingTestMatrix.n(); ++jjj)
		{
			accumulator += CouplingTestMatrix(iii,jjj);
		}
	}

	libMesh::Real vol = 0;
	libMesh::Elem* silly_elem;
	for(libMesh::MeshBase::element_iterator itBegin = IntersectionMesh.elements_begin();
											itBegin != IntersectionMesh.elements_end();
											++itBegin)
	{
		silly_elem = *itBegin;
		vol += silly_elem->volume();
	}

	libMesh::Real difference = accumulator - n_var*CouplingScale*vol;

	std::cout << "|" << std::endl;
	std::cout << "|    Sum( M_i,j )   = " << accumulator << std::endl;
	std::cout << "|    n * C * Volume = " << n_var*CouplingScale*vol << std::endl;
	std::cout << "| >  Difference     = " << difference << std::endl << std::endl;


}
