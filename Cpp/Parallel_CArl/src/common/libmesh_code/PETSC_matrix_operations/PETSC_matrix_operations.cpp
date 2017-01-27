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

void carl::print_matrix_dim(libMesh::PetscMatrix<libMesh::Number>& CouplingTestMatrix, bool bDetailed)
{
	std::cout << "| M_i,j  : " << CouplingTestMatrix.m() << " x " << CouplingTestMatrix.n() << std::endl;
//	if(bDetailed)
//	{
		MatInfo temp_info;
		MatGetInfo(CouplingTestMatrix.mat(),MAT_LOCAL,&temp_info);
		std::cout << "| LOCAL  : memory = " << temp_info.memory << std::endl;
		std::cout << "|          non-zeros used = " << (100.*temp_info.nz_used)/temp_info.nz_allocated << " % " << std::endl;

		MatGetInfo(CouplingTestMatrix.mat(),MAT_GLOBAL_SUM,&temp_info);
		std::cout << "| GLOBAL : memory = " << temp_info.memory << std::endl;
		std::cout << "|          non-zeros used = " << (100.*temp_info.nz_used)/temp_info.nz_allocated << " % " << std::endl;
		std::cout << "|          total nb. of non-zeros used = " << temp_info.nz_used << std::endl;
		std::cout << "|          total nb. of non-zeros allocated = " << temp_info.nz_allocated << std::endl;
		std::cout << "|          sparsity = " << (100.*temp_info.nz_used)/(CouplingTestMatrix.m() * CouplingTestMatrix.n()) << " % " << std::endl;

		MatGetInfo(CouplingTestMatrix.mat(),MAT_GLOBAL_MAX,&temp_info);
		std::cout << "| MAX    : memory = " << temp_info.memory << std::endl;
		std::cout << "|          non-zeros used = " << (100.*temp_info.nz_used)/temp_info.nz_allocated << " % " << std::endl << std::endl;
	//}
}

void carl::print_matrix_info(libMesh::PetscMatrix<libMesh::Number>& InputMatrix, std::ostream & os)
{
	const libMesh::Parallel::Communicator& WorldComm = InputMatrix.comm();

	unsigned int rank = WorldComm.rank();
	unsigned int nodes = WorldComm.size();

	MatInfo local_info;
	MatInfo global_info;
	MatGetInfo(InputMatrix.mat(),MAT_LOCAL,&local_info);
	MatGetInfo(InputMatrix.mat(),MAT_LOCAL,&global_info);

	// Set up local variables
	int l_nz_used = local_info.nz_used;
	int l_nz_allocated = local_info.nz_allocated;
	int l_memory = local_info.memory;

	int l_n = -1;
	int l_m = -1;

	double accumulator_n = 0;
	double accumulator_m = 0;

	MatGetLocalSize(InputMatrix.mat(),&l_n,&l_m);

	std::vector<int>	full_nz_used;
	std::vector<int>	full_nz_allocated;
	std::vector<int>	full_memory;

	std::vector<int>	full_n;
	std::vector<int>	full_m;

	if(rank == 0)
	{
		full_nz_used.resize(nodes,0);
		full_nz_allocated.resize(nodes,0);
		full_memory.resize(nodes,0);
		full_n.resize(nodes,0);
		full_m.resize(nodes,0);
	}

	WorldComm.gather(0,l_nz_used,full_nz_used);
	WorldComm.gather(0,l_nz_allocated,full_nz_allocated);
	WorldComm.gather(0,l_memory,full_memory);
	WorldComm.gather(0,l_n,full_n);
	WorldComm.gather(0,l_m,full_m);

	if(rank == 0)
	{
		os << "# rank \t local_n \t local_m \t start_n \t start_m \t nz_alloc \t nz_used \t memory  " << std::endl;
		for(unsigned int iii = 0; iii < nodes; ++iii)
		{
			os 	<< iii << " \t "
				<< full_n[iii] << " \t "
				<< full_m[iii] << " \t "
				<< accumulator_n << " \t "
				<< accumulator_m << " \t "
				<< full_nz_allocated[iii] << " \t "
				<< full_nz_used[iii] << " \t "
				<< full_memory[iii] << std::endl;
			accumulator_n += full_n[iii];
			accumulator_m += full_m[iii];
		}
	}

	WorldComm.barrier();
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

void carl::write_PETSC_vector(	libMesh::PetscVector<libMesh::Number>& input_vec,
		const std::string& filename)
{
	const libMesh::Parallel::Communicator& WorldComm = input_vec.comm();

	PetscViewer    viewer;
	PetscViewerBinaryOpen(WorldComm.get(),filename.c_str(),FILE_MODE_WRITE,&viewer);
	VecView(input_vec.vec(),viewer);

	PetscViewerDestroy(&viewer);
}

void carl::print_PETSC_vector(	libMesh::PetscVector<libMesh::Number>& input_vec,
		const std::string& filename)
{
	const libMesh::Parallel::Communicator& WorldComm = input_vec.comm();

	PetscViewer    viewer;
	PetscViewerASCIIOpen(WorldComm.get(),filename.c_str(),&viewer);
	VecView(input_vec.vec(),viewer);

	PetscViewerDestroy(&viewer);
}

void carl::read_PETSC_vector(	libMesh::PetscVector<libMesh::Number>& input_vec,
		const std::string& filename)
{
	const libMesh::Parallel::Communicator& WorldComm = input_vec.comm();

	PetscViewer    viewer;
	PetscViewerBinaryOpen(WorldComm.get(),filename.c_str(),FILE_MODE_READ,&viewer);
	VecLoad(input_vec.vec(),viewer);

	PetscViewerDestroy(&viewer);
}

void carl::create_PETSC_dense_matrix_from_vectors(const Vec *vecs_in, int nb_vecs, Mat& matrix_out)
{
	// Get the vectors' dimensions and create the matrix
	PetscInt M_local, M;
	VecGetLocalSize(vecs_in[0],&M_local);
	VecGetSize(vecs_in[0],&M);

	MatCreateDense(PETSC_COMM_WORLD,M_local,PETSC_DECIDE,M,nb_vecs,NULL,&matrix_out);

	// Get the ownership ranges and the row indexes
	PetscInt own_low, own_high;
	VecGetOwnershipRange(vecs_in[0],&own_low,&own_high);

	std::vector<PetscInt> row(M_local,0);

	for(int iii = 0; iii < M_local; ++iii)
	{
		row[iii] = own_low + iii;
	}

	// Set the values
	PetscScalar* data;
	for(int iii = 0; iii < nb_vecs; ++iii)
	{
		VecGetArray(vecs_in[iii],&data);
		MatSetValues(matrix_out,M_local,row.data(),1,&iii,data,INSERT_VALUES);
		VecRestoreArray(vecs_in[iii],&data);
	}

	// Assembly
	MatAssemblyBegin(matrix_out,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(matrix_out,MAT_FINAL_ASSEMBLY);
}

void carl::PETSC_invert_dense_matrix(Mat& matrix_in, Mat& matrix_out)
{
	// WARNING: this operation is extremely expensive for large systems !!!
	//          Only use this for small matrices (example: the 6x6 or 3x3
	//          or 3x3 matrices from the null space projectors,

	Mat				Id_mat;
	MatFactorInfo 	factor_info;
	IS				rperm, cperm;

	// Duplicate matrix data structures
	MatDuplicate(matrix_in,MAT_DO_NOT_COPY_VALUES,&matrix_out);
	MatDuplicate(matrix_in,MAT_DO_NOT_COPY_VALUES,&Id_mat);

	// Create identity
	MatZeroEntries(Id_mat);
	MatShift(Id_mat,1);

	// Factor input matrix
	MatGetOrdering(matrix_in,MATORDERINGNATURAL,  &rperm,  &cperm);
	MatFactorInfoInitialize(&factor_info);
	MatLUFactor(matrix_in,rperm,cperm,&factor_info);

	// Calculate inverse
	MatMatSolve(matrix_in,Id_mat,matrix_out);

	// Reset input's factoring
	MatSetUnfactored(matrix_in);

	// Cleanup
	MatDestroy(&Id_mat);
}
