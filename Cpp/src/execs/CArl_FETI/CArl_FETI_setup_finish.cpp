#include "CArl_FETI_setup_finish.h"

/** \file CArl_FETI_setup_finish.cpp
\brief **STAT** Program responsible to finish the FETI setup and launch the iterations

This program's input file description can be found at the documentation of the function 
  carl::get_input_params(GetPot& field_parser, feti_setup_finish_params& input_params) .
  
  It will use the following files ... 
 - ... from the `input_params.coupling_folder_path` folder:
   + coupling matrices \f$C_1\f$, \f$C_2\f$ and \f$C_{\mbox{Med}}\f$ (the latter used for the preconditioner). *Files*:

         coupling_matrix_macro.petscmat
         coupling_matrix_micro.petscmat
         coupling_matrix_mediator.petscmat

 -  ... from the `input_params.scratch_folder_path` folder:
   + solutions \f$u_0,1\f$ and \f$u_{0,2}\f$, from the system \f$K_l * u_{0,l} = F_l\f$. *Files*:

         ext_solver_u0_A_sys_sol_vec.petscvec
         ext_solver_u0_B_sys_sol_vec.petscvec

   + [RB] solutions \f$x_1(kkk)\f$ and \f$x_2(kkk)\f$, from the system \f$K_l * x_l(kkk) = C_l^t*\phi_0\f$. *Files*:

         ext_solver_A_sys_sol_vec.petscvec
         ext_solver_B_sys_sol_vec.petscvec

   + [RB] matrix \f$\mbox{inv}(R_I^t * R_I) = \mbox{inv}(R_2^t*C_2^t*C_2*R_2)\f$, used for the rigid body modes projections. *Files*:

         rb_inv_RITRI.petscmat

   + [RB] rigid body mode vectors multiplied by \f$C_2\f$. *Files*:

         rb_coupl_vector_[iii]_n_[nb. of vectors].petscvec

 - ... from the micro system folder (common vector path given by `input_params.RB_vectors_base`):
   + [RB] rigid body mode vectors. *Files*:

         [input_params.RB_vectors_base]_rb_vector_[iii]_n_[nb. of vectors].petscvec

The items marked with a [RB] are only needed if the rigid body modes projectors are used.
In the last two cases, [nb. of vectors] is the number of rigid body mode vectors (given by `input_params.nb_of_rb_vectors`)
and [iii] is an integer going from 0 to `input_params.nb_of_rb_vectors - 1` (following C++ notation).

This program outputs a series of files, all inside the `input_params.scratch_folder_path` folder:
  + initial iteration vectors, \f$r(0)\f$, \f$z(0)\f$, \f$p(0)\f$. *Files*:

        FETI_iter__r__0.petscvec
        FETI_iter__p__0.petscvec

  + (overwrite) vectors used as the RHS for the external solvers. *Files*:

        ext_solver_A_rhs.petscvec
        ext_solver_B_rhs.petscvec

  + (create) scalar values (iteration, residual, RB mode corrections). *Files*:

        FETI_iter_scalar_data.dat

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

	carl::feti_setup_finish_params input_params;
	get_input_params(field_parser, input_params);

	// Object containing the FETI operations
	carl::FETI_Operations feti_op(WorldComm,input_params.scratch_folder_path,input_params.coupling_folder_path);

	// --- Define if the rb modes will be used or not
	feti_op.using_rb_modes(input_params.bUseRigidBodyModes);

	// --- Read the files!

	// Read up the coupling matricesconst std::string& filename)
	feti_op.set_coupling_matrix_R_micro();
	feti_op.set_coupling_matrix_R_BIG();

	// Read the decoupled solutions, K_i * u_0,i  = F_i 
	feti_op.read_decoupled_solutions();

	// Read operations needed if we are using the rigid body modes
	if(input_params.bUseRigidBodyModes)
	{
		// Read the solutions of K_i * x_i(0)  = C_i^t * phi(0)
		feti_op.read_ext_solver_output();

		// Read the RB-related vectors and matrices
		feti_op.read_null_space_vecs(input_params.RB_vectors_base,input_params.nb_of_rb_vectors);
		feti_op.read_null_space_inv_RITRI_mat();
	}

	// --- Set up any matrices or vectors needed before calculating the outputs
	// Set up the preconditioner
	feti_op.set_preconditioner(input_params.CG_precond_type /*, initial_set = true */ );

	// --- Calculate the output vectors! All are saved internaly inside the object
	// Calculate r(0)
	feti_op.calculate_initial_r();

	// Calculate p(0)
	feti_op.calculate_initial_p();

	// Calculations needed if we are using the rigid body modes
	if(input_params.bUseRigidBodyModes)
	{
		// Calculate the rigid body modes correction RB_corr
		feti_op.calculate_rb_correction();
	}

	// --- Export output vectors!
	// Export 'r(0)' and 'p(0)'
	feti_op.export_inital_vecs();

	// Export the Ct_i * p(0) vectors
	feti_op.export_ext_solver_rhs_initial();

	// Export the scalar data, rho(0) and, if pertinent, |RB_corr|
	feti_op.export_initial_scalar_data();

	// --- Launch the "iter_script.sh" script --- ONLY ON THE FIRST PROC!
	if(WorldComm.rank() == 0)
	{
		std::string iter_script_command = ". " + input_params.scratch_folder_path + "/FETI_iter_script.sh";
		if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
		{
			std::cout << " !!! LOCAL test: MPI commands cannot be launched recursivelly !!! " << std::endl;
			std::cout << "     Run the following program by hand: " << std::endl;
			std::cout << iter_script_command << std::endl;
		} else {
			carl::exec_command(iter_script_command);
		}
	}

	return 0;
}
