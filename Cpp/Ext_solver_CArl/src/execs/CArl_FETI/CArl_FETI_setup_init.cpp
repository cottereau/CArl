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
	
	// --- Crete the files / folders needed
	// Create the scratch folder
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
	
	// --- Calculate phi(0), if needed
	if(input_params.bUseRigidBodyModes)
	{
		carl::FETI_Operations feti_op(WorldComm,input_params.scratch_folder_path,input_params.coupling_folder_path);

		// --- Define if the rb modes will be used or not
		feti_op.using_rb_modes(input_params.bUseRigidBodyModes);

		// Read the matrices
		feti_op.set_coupling_matrix_R_micro();
		feti_op.set_coupling_matrix_R_BIG();

		// Set the null space
		feti_op.set_null_space(input_params.RB_vectors_base,input_params.nb_of_rb_vectors);

		// Calculate the inital solution
		feti_op.calculate_null_space_phi_0(input_params.force_micro_path);

		// Export the vectors
		feti_op.export_phi();
		feti_op.export_ext_solver_rhs_Ct_phi();
	}

	// --- Launch the "init_script.sh" script --- ONLY ON THE FIRST PROC!
	if(WorldComm.rank() == 0)
	{
		std::string init_script_command = ". " + input_params.scratch_folder_path + "/FETI_init_script.sh";
		if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
		{
			std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
			std::cout << init_script_command << std::endl << std::endl;
		} else {
			carl::exec_command(init_script_command);
		}
	}

	return 0;
}
