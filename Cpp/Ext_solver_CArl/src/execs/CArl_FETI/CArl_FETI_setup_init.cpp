#include "CArl_FETI_setup_init.h"

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

	carl::feti_setup_init_params input_params;
	get_input_params(field_parser, input_params);
	carl::Solver_Files_Setup FETI_files_setup(WorldComm,input_params);
	
	// Create the scratch folder, if needed
	FETI_files_setup.set_scratch_folder();

	// [LIBMESH] Create the external solver input files
	FETI_files_setup.generate_libmesh_external_solver_inputs();

	// [LIBMESH] Create the external solver scripts
	FETI_files_setup.generate_libmesh_external_solver_scripts();

	// Create FETI input files
	FETI_files_setup.generate_FETI_inputs();

	// Create FETI script files
	FETI_files_setup.generate_FETI_scripts();

	// Create FETI lauch script files
	FETI_files_setup.generate_FETI_launch_scripts();
	
	/* --- What this program must do:
	 *
	 * 1) Generate the qsub files for the following external solvers - DONE
	 *    K_i * u_0,i  = F_i
	 *    K_i * x_0,i  = C_i^T * phi(0)
	 *    K_i * y(k)_i = C_i^T * p(k)
	 *    K_i * x_f,i  = C_i^T * phi(k+1)
	 *
	 * 2) Generate the input files for the other FETI - DONE
	 *	  CArl_FETI_setup_finish
	 *    CArl_FETI_iterate
	 *	  CArl_FETI_set_sol
	 * 
	 * 3) Generate the qsub files for the other FETI programs - DONE
	 *
	 * 4) Generate the scripts to launch the programs - DONE
	 *	  init_script
	 *    iter_script
	 *	  sol_script
	 * 
	 * 5) Read F_2, C_2 and RB to generate phi(0)
	 *
	 * >>> Create a carl::solver_files_setup class to do 1), 2), 3) and 4)
	 * >>> Create a carl::FETI_operations class to do 5) (and the other programs steps)
	 *
	 */


	return 0;
}