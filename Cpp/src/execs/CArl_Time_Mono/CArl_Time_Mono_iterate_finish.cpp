#include "CArl_Time_Mono_iterate_finish.h"

/** \file CArl_Time_Mono_iterate_finish.cpp
\brief Program responsible to initialize the Time_Mono iterate finish and launch the iterations

 */

int main(int argc, char** argv) {

    // --- Initialize libMesh
	libMesh::LibMeshInit init(argc, argv);

	// Do performance log?
	const bool MASTER_bPerfLog_carl_time_mono_iterate_finish = true;
	libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_time_mono_iterate_finish);

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

	carl::time_mono_iterate_finish_params input_params;
	get_input_params(field_parser, input_params);

	carl::Time_Mono_Iterate time_mono_it(WorldComm, input_params);
	
	// Read the state of the previous iteration
	perf_log.push("Read and update state");
	time_mono_it.read_and_update_state(); // Bien fait sur le proc 0
	perf_log.pop("Read and update state");

	// Copie des champs sortis de FETI dans le dossier résultats
	perf_log.push("Copy FETI outputs");
	time_mono_it.copy_FETI_outputs(); // Bien fait sur le proc 0
	perf_log.pop("Copy FETI outputs");

	// Check if the time loop is completed
	bool time_loop_completed = time_mono_it.check_time_loop_completed();

	if (!time_loop_completed)
	{
		// --- Launch the "Time_Mono_iterate_init.sh" script --- ONLY ON THE FIRST PROC!
		if(WorldComm.rank() == 0)
		{
			std::string iterate_init_script_command = ". " + input_params.scratch_folder_path + "/Time_Mono_iterate_init.sh";
			if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
			{
				std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
				std::cout << iterate_init_script_command << std::endl << std::endl;
			} else {
				carl::exec_command(iterate_init_script_command);
			}
		}

	}
	else
	{
		std::cout << "Time loop completed. Please check the results." << std::endl;
	}

	return 0;
}