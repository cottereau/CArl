#include "CArl_FETI_iterate.h"

/**	\brief Program responsible to running the FETI iterations
 *
 *	This program's input file description can be found at the documentation of the function 
 *  carl::get_input_params(GetPot& field_parser, feti_iterate_params& input_params).
 *  
 *  It will use the following files ... 
 *  * ... from the `input_params.coupling_folder_path` folder:
 *    + coupling matrices C_1 and C_2. *Files*:
 *      - `coupling_matrix_macro.petscmat`
 *      - `coupling_matrix_micro.petscmat`
 *
 *  * ... from the `input_params.scratch_folder_path` folder:
 *    + scalar values (iteration, residual, RB mode corrections). *Files*:
 *      - `FETI_iter_scalar_data.dat`
 *    + solutions x_1(kkk) and x_2(kkk), from the system K_i * x_i(kkk) = C_i^T*p(kkk). *Files*:
 *      - `ext_solver_A_sys_sol_vec.petscvec`
 *      - `ext_solver_B_sys_sol_vec.petscvec`
 *    + previous iteration vectors r(kkk) anf phi(kkk)
 *      - `FETI_iter__phi__current.petscvec`
 *	    - `FETI_iter__r__current.petscvec`
 *    + (several) previous iterations vectors p(jjj), q(jjj) (used for re-orthogonalization). *Files*:
 *      - `FETI_iter__q__[jjj].petscvec`, jjj = 0 ... kkk - 1
 *      - `FETI_iter__p__[jjj].petscvec`, jjj = 0 ... kkk
 *    + previous p(jjj).q(jjj) values (with jjj = 0 ... kkk - 1). *Files*:
 *      - `FETI_iter_p_dot_q.dat`
 *    + matrix inv (R_I^t * R_I) = inv(R_2^t*C_2^t*C_2*R_2), used for the rigid body modes projections. *Files*:
 *      - [RB] `rb_inv_RITRI.petscmat`.
 *    + rigid body mode vectors multiplied by C_2. *Files*:
 *      - [RB] `rb_coupl_vector_[iii]_n_[nb. of vectors].petscvec`
 *  * ... from the micro system folder (common vector path given by `input_params.RB_vectors_base`):
 *    + rigid body mode vectors. *Files*:
 *      - [RB] `[input_params.RB_vectors_base]_rb_vector_[iii]_n_[nb. of vectors].petscvec`
 *
 * The items marked with a [RB] are only needed if the rigid body modes projectors are used.
 * In the last two cases, [nb. of vectors] is the number of rigid body mode vectors (given by `input_params.nb_of_rb_vectors`)
 * and [iii] is an integer going from 0 to `input_params.nb_of_rb_vectors - 1` (following C++ notation).
 *
 * This program outputs a series of files, all inside the `input_params.scratch_folder_path` folder:
 *   + (append) scalar values (iteration, residual, RB mode corrections). *Files*:
 *     - `FETI_iter_scalar_data.dat`
 *   + (overwrite) vectors used as the RHS for the external solvers. *Files*:
 *     - `ext_solver_A_rhs.petscvec`
 *     - `ext_solver_B_rhs.petscvec`
 *   + next iteration vectors, phi(kkk+1), r(kkk+1), q(kkk+1), p(kkk+1). *Files*:
 *     - `FETI_iter__phi__current.petscvec`
 *     - `FETI_iter__r__current.petscvec`
 *     - `FETI_iter__q__[kkk+1].petscvec`
 *     - `FETI_iter__p__[kkk+1].petscvec`
 *   + (append) p(kkk).q(kkk) value. *Files*:
 *     - `FETI_iter_p_dot_q.dat`
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

	carl::feti_iterate_params input_params;
	get_input_params(field_parser, input_params);

	// Object containing the FETI operations
	carl::FETI_Operations feti_op(WorldComm,input_params.scratch_folder_path,input_params.coupling_folder_path);

	// --- Define if the rb modes will be used or not
	feti_op.using_rb_modes(input_params.bUseRigidBodyModes);

	// --- Read the common files: coupling matrices, null space vectors ...
	// Read up the coupling matricesconst std::string& filename)
	feti_op.set_coupling_matrix_R_micro();
	feti_op.set_coupling_matrix_R_BIG();

	// --- Set up any matrices or vectors needed before calculating the outputs
	// Set up the preconditioner
	feti_op.set_preconditioner(input_params.CG_precond_type, /* initial_set = */ false);

	// Read operations needed if we are using the rigid body modes
	if(input_params.bUseRigidBodyModes)
	{
		// Read the RB-related vectors and matrices
		feti_op.read_null_space_vecs(input_params.RB_vectors_base,input_params.nb_of_rb_vectors);
		feti_op.read_null_space_inv_RITRI_mat();
	}

	// --- Read the previous iteration files: scalar data, iteration vectors ...
	/* Read the scalar data from the previous iterations '0 ... kkk'
	 * We now have: 'kkk'
	 *				'rho(0)'
	 *				'rho(kkk)'
	 *				'| RB_corr(kkk) |'
	 *				'p(0 ... kkk - 1).q(0 ... kkk - 1)'
	 */
	feti_op.read_scalar_data();

	// Read the vector data from the previous iterations '0 ... kkk'
	/* We now have: 'r(kkk)'
	 *				'phi(kkk)'
	 *				'p(0 ... kkk)'
	 * 				'q(0 ... kkk - 1)'
	 */
	feti_op.read_vector_data();

	// Read the previous iteration 'kkk' external solver output, 'x_i(kkk)'
	feti_op.read_ext_solver_output();

	// --- Iterate!
	// Calculate 'q(kkk) = C_1 * x_1(kkk) + C_2 * x_2(kkk)' and 'p(kkk).q(kkk)'
	feti_op.calculate_q();

	// Calculate 'phi(kkk + 1) = phi(kkk) + gamma * p(kkk)'
	// gamma = rho(kkk) / ( p(kkk).q(kkk) )
	feti_op.calculate_phi();

	// Calculate 'r(kkk + 1) = r(kkk) - gamma * q(kkk)'
	feti_op.calculate_r();

	// Calculate 'z(kkk + 1)' (formula depends on preconditioner and projection settings)
	feti_op.calculate_z();

	// Calculate 'p(kkk + 1)'
	feti_op.calculate_p();

	// Calculate 'RB_corr(kkk+1)'
	if(input_params.bUseRigidBodyModes)
	{
		feti_op.calculate_rb_correction();
	}

	// Calculate the scalar data, 'rho(kkk+1)' and '| RB_corr(kkk+1) |'
	feti_op.calculate_scalar_data();

	/* Export scalar data
	 * Data to export: 'kkk+1'
	 *				   'rho(0)'
	 *				   'rho(kkk+1)'
	 *				   '| RB_corr(kkk+1) |'
	 *				   'p(kkk).q(kkk)'
	 */
	feti_op.export_scalar_data();

	/* Export the iteration vectors
	 * Vectors to export: 'r(kkk+1)'
	 *                    'phi(kkk+1)'
	 *                    'p(kkk+1)'
	 *                    'q(kkk)'
	 */
	feti_op.export_iter_vecs();

	// // --- Check the convergence
	carl::IterationStatus current_iteration_status = carl::IterationStatus::ITERATING;
	
	if(input_params.bUseRigidBodyModes)
	{
		current_iteration_status = feti_op.check_convergence(input_params.CG_coupled_conv_rel, input_params.CG_coupled_conv_abs, input_params.CG_coupled_conv_max, input_params.CG_coupled_div, input_params.CG_coupled_conv_corr);
	} else {
		current_iteration_status = feti_op.check_convergence(input_params.CG_coupled_conv_rel, input_params.CG_coupled_conv_abs, input_params.CG_coupled_conv_max, input_params.CG_coupled_div);
	}

	// Print the current values of the convergence parameters
	feti_op.print_previous_iters_conv( /* nb. of iterations = 5 */);

	switch (current_iteration_status)
	{
		case carl::IterationStatus::ITERATING :
				// --- Continue the iteration

				// Export the Ct_i * p(kkk+1) vectors
				feti_op.export_ext_solver_rhs_Ct_p();

				// --- Launch the "iter_script.sh" script --- ONLY ON THE FIRST PROC!
				if(WorldComm.rank() == 0)
				{
					std::string iter_script_command = ". " + input_params.scratch_folder_path + "/FETI_iter_script.sh";
					if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
					{
						std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
						std::cout << iter_script_command << std::endl << std::endl;
					} else {
						carl::exec_command(iter_script_command);
					}
				}
				break;
		case carl::IterationStatus::CONVERGED :
				// --- Well ... converged!

				// Export the Ct_i * phi(kkk+1) vectors
				feti_op.export_ext_solver_rhs_Ct_phi();

				// Export the rigid body modes correction vector
				feti_op.export_rb_correction_vector();
				
				// --- Launch the "sol_script.sh" script --- ONLY ON THE FIRST PROC!
				if(WorldComm.rank() == 0)
				{
					std::string sol_script_command = ". " + input_params.scratch_folder_path + "/FETI_sol_script.sh";
					if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
					{
						std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
						std::cout << sol_script_command << std::endl << std::endl;
					} else {
						carl::exec_command(sol_script_command);
					}
				}
				break;
		case carl::IterationStatus::DIVERGED :
				// --- Well, we have to stop here ...
				break;
	}

	return 0;
}
