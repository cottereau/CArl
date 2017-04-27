#include "CArl_FETI_solution.h"

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

	carl::feti_set_sol_params input_params;
	get_input_params(field_parser, input_params);

	return 0;
}
