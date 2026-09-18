#include "CArl_Time_Mono_setup.h"

/** \file CArl_Time_Mono_setup.cpp
\brief Program responsible to initialize the Time_Mono setup and launch the iterations

 */

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

	carl::time_mono_setup_params input_params;
	get_input_params(field_parser, input_params);
	carl::Time_Mono_Files_Setup Time_Mono_Files_Setup(WorldComm, input_params);
	
	// --- Crete the files / folders needed
	// Create the scratch folder
	Time_Mono_Files_Setup.set_scratch_folder();

	// Create the scratch folder
	Time_Mono_Files_Setup.generate_state_file();

	// Create the folders containing the state of each domain
	Time_Mono_Files_Setup.generate_domain_state_folders();

	// Create the results folder
	Time_Mono_Files_Setup.set_results_folder();

	// Copy the initial conditions files to the scratch folder and results folder
	Time_Mono_Files_Setup.copy_initial_conditions_files();
	
	// Generating the inputs and scripts of the external solver (libMesh)
	Time_Mono_Files_Setup.generate_libmesh_external_solver_assemble_rhs_inputs();
	Time_Mono_Files_Setup.generate_libmesh_external_solver_assemble_rhs_scripts();
	Time_Mono_Files_Setup.generate_libmesh_external_solver_update_init_cond_inputs();
	Time_Mono_Files_Setup.generate_libmesh_external_solver_update_init_cond_scripts();

	// Generating the FETI setup init script
	Time_Mono_Files_Setup.generate_FETI_setup_init_scripts();

	// Generating the Time_Mono inputs
	Time_Mono_Files_Setup.generate_Time_Mono_inputs();

	// Generating the Time_Mono scripts
	Time_Mono_Files_Setup.generate_Time_Mono_scripts();

	// Generating the Time_Mono launch scripts
	Time_Mono_Files_Setup.generate_Time_Mono_launch_scripts();

	// --- Launch the "Time_Mono_iterate_init.sh" script --- ONLY ON THE FIRST PROC!
	if(WorldComm.rank() == 0)
	{
		std::string next_script_command = ". " + input_params.scratch_folder_path + "/Time_Mono_iterate_init.sh";
		if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL) 
		{
			std::cout << " !!! LOCAL test: MPI commands cannot be launched recursivelly !!! " << std::endl;
			std::cout << "     Run the following program by hand: " << std::endl;
			std::cout << next_script_command << std::endl;
		}
		else 
		{
			carl::exec_command(next_script_command);
		}
	}

	return 0;
}