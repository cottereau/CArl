#include "CArl_mult_coupling.h"

int main(int argc, char** argv) {

	// --- Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);

	// Do performance log?
	libMesh::PerfLog perf_log("Main program");

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

	carl::carl_mult_coupling_input_params input_params;
	get_input_params(field_parser, input_params);

	// Check libMesh installation dimension
	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// --- Set the matrix and vectors

	// Set up the PETSC versions
	Mat coupl_mat_PETSC;
	Vec input_vec_PETSC;

	MatCreate(PETSC_COMM_WORLD,&coupl_mat_PETSC);
	VecCreate(PETSC_COMM_WORLD,&input_vec_PETSC);
	
	// Read
	carl::read_PETSC_matrix(coupl_mat_PETSC, input_params.coupl_matrix_file, WorldComm.get());
	carl::read_PETSC_vector(input_vec_PETSC, input_params.input_vec_file, WorldComm.get());

	// Set up the libMesh versions
	libMesh::PetscMatrix<libMesh::Number> coupl_mat(coupl_mat_PETSC,WorldComm);
	libMesh::PetscVector<libMesh::Number> input_vec(input_vec_PETSC,WorldComm);

	// Set the output vector
	PetscInt N, local_N;
	MatGetSize(coupl_mat_PETSC,&N,NULL);
	MatGetLocalSize(coupl_mat_PETSC,&local_N,NULL);

	libMesh::PetscVector<libMesh::Number> output_vec(WorldComm,N,local_N);

	// --- Multiply ...
	coupl_mat.vector_mult(output_vec,input_vec);
	
// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
	output_vec.print_matlab(input_params.output_base + "_output_vec.m");
#endif

	// Export the solution vector
	carl::write_PETSC_vector(output_vec, input_params.output_base + "_output_vec.petscvec");

	// --- Cleanup!
	MatDestroy(&coupl_mat_PETSC);
	VecDestroy(&input_vec_PETSC);

	return 0;
}
