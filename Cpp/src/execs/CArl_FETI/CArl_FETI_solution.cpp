#include "CArl_FETI_solution.h"

/** \file CArl_FETI_solution.cpp
\brief **STAT** Program responsible to calculate the coupled system solution of the FETI algorithm

This program's input file description can be found at the documentation of the function carl::get_input_params(GetPot& field_parser, feti_set_sol_params& input_params).
  
It will use the following files from the `input_params.scratch_folder_path` folder:
 + solutions \f$u_{0,1}\f$ and \f$u_{0,}2\f$, from the system \f$K_l * u_{0,l} = F_l\f$. *Files*:

       ext_solver_u0_A_sys_sol_vec.petscvec
       ext_solver_u0_B_sys_sol_vec.petscvec

 + solutions \f$x_1(kkk) and \f$x_2(kkk), from the system \f$K_l * x_l(kkk) = C_l^t*phi(kkk)\f$. *Files*:

       ext_solver_A_sys_sol_vec.petscvec
       ext_solver_B_sys_sol_vec.petscvec

 + rigid body mode correction vector. *Files*:

       [RB] FETI_RB_correction.petscvec

The items marked with a [RB] are only needed if the rigid body modes projectors are used.
In the last two cases, [nb. of vectors] is the number of rigid body mode vectors (given by `input_params.nb_of_rb_vectors`)
and [iii] is an integer going from 0 to `input_params.nb_of_rb_vectors - 1` (following C++ notation).

This program outputs the solution vectors inside the `input_params.output_base` folder. *Files*:

    [input_params.output_base]_coupled_sol_A.petscvec
    [input_params.output_base]_coupled_sol_B.petscvec
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

	carl::feti_set_sol_params input_params;
	get_input_params(field_parser, input_params);

	// Object containing the FETI operations
	carl::FETI_Operations feti_op(WorldComm,input_params.scratch_folder_path);

	// --- Define if the rb modes will be used or not
	feti_op.using_rb_modes(input_params.bUseRigidBodyModes);

	// Read the decoupled solutions, K_i * u_0,i  = F_i 
	feti_op.read_decoupled_solutions();

	// Read the solutions of K_i * x_i(FINAL)  = C_i^t * phi(FINAL)
	feti_op.read_ext_solver_output();

	// Read the rb modes correction, 'RB_corr(FINAL)'
	if(input_params.bUseRigidBodyModes)
	{
		feti_op.read_rb_corr();
	}

	// Calculate the solution
	feti_op.calculate_coupled_solution();

	// Export it (finallly!)
	feti_op.export_coupled_solution(input_params.output_folder);

	return 0;
}
