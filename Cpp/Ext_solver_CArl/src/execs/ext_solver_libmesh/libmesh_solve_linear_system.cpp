#include "libmesh_solve_linear_system.h"

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

	libmesh_solve_linear_system_input_params input_params;
	get_input_params(field_parser, input_params);

	// Check libMesh installation dimension
	const unsigned int dim = 3;

	libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

	// --- Set the matrix and vectors

	// Set up the PETSC versions
	Mat sys_mat_PETSC;
	Vec sys_rhs_vec_PETSC;

	MatCreate(PETSC_COMM_WORLD,&sys_mat_PETSC);
	VecCreate(PETSC_COMM_WORLD,&sys_rhs_vec_PETSC);
	
	// Read
	carl::read_PETSC_matrix(sys_mat_PETSC, input_params.sys_matrix_file, WorldComm.get());
	carl::read_PETSC_vector(sys_rhs_vec_PETSC, input_params.sys_rhs_vec_file, WorldComm.get());

	// Set up the libMesh versions
	libMesh::PetscMatrix<libMesh::Number> sys_mat(sys_mat_PETSC,WorldComm);
	libMesh::PetscVector<libMesh::Number> sys_rhs_vec(sys_rhs_vec_PETSC,WorldComm);
	libMesh::PetscVector<libMesh::Number> sys_sol_vec(WorldComm);
	sys_sol_vec.init(sys_rhs_vec);

	// --- Linear solver
	libMesh::PetscLinearSolver<libMesh::Number> KSP_solver(WorldComm);
	KSP_solver.init("sys");

	// Solve!
	KSP_solver.solve(sys_mat,sys_sol_vec,sys_rhs_vec,input_params.sys_eps,input_params.sys_iter_div);

	// Export the solution vector
	sys_sol_vec.print_matlab(input_params.output_base + "_sys_sol_vec.m");
	carl::write_PETSC_vector(sys_sol_vec.vec(), input_params.output_base + "_sys_sol_vec.petscmat",WorldComm.get());

	// --- Cleanup!
	MatDestroy(&sys_mat_PETSC);
	VecDestroy(&sys_rhs_vec_PETSC);
}
