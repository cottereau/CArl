#include "CArl_FETI_setup_init.h"

/** \file CArl_FETI_setup_init.cpp
\brief Program responsible to initialize the FETI setup and launch the iterations

This program's input file description can be found at the documentation of the function carl::get_input_params(GetPot& field_parser, feti_setup_init_params& input_params). 

If the rigid body modes are used (`input_params.bUseRigidBodyModes = true`), it will use the following files ...
 - ... from the `input_params.coupling_folder_path` folder:
   + coupling matrices \f$C_1\f$ and \f$C_2\f$. *Files*:

         coupling_matrix_macro.petscmat
         coupling_matrix_micro.petscmat

 - ... from the micro system folder 
   + external force work. *File*:

         [input_params.force_micro_path]

   + rigid body mode vectors. *Files*:

         [input_params.RB_vectors_base]_rb_vector_[iii]_n_[nb. of vectors].petscvec

It will create a scratch folder with the path `input_params.scratch_folder_path`, and create the following files:

 - Scripts to launch the other solver steps:

       FETI_init_script.sh
       FETI_iter_script.sh
       FETI_sol_script.sh

 - Scripts to launch the external solvers:

       ext_solver_A.sh
       ext_solver_B.sh
       ext_solver_u0_A.sh
       ext_solver_u0_B.sh

 - Scripts to launch the other `CArl_FETI_***` binaries:

       CArl_FETI_setup_finish.sh
       CArl_FETI_iterate.sh
       CArl_FETI_solution.sh

 - Input files for the external solvers and other `CArl_FETI_***` binaries:

       ext_solver_A.txt
       ext_solver_B.txt
       ext_solver_u0_A.txt
       ext_solver_u0_B.txt
       CArl_FETI_iterate.txt
       CArl_FETI_solution.txt
       CArl_FETI_setup_finish.txt

 - Right hand side vectors for the external solvers:

       ext_solver_A_rhs.petscvec
       ext_solver_B_rhs.petscvec

 - If the rigid body modes are used ... 
   + [RB] initial solution (the Lagrange multipliers)

         FETI_iter__phi__current.petscvec

   + [RB] products between \f$C_2\f$ and the rigid body modes vectors ([nb. of vectors] is the number of rigid body mode vectors (given by `input_params.nb_of_rb_vectors`) and [iii] is an integer going from 0 to `input_params.nb_of_rb_vectors - 1`, following C++ notation)

         rb_coupl_vector_[iii]_n_[nb. of vectors].petscvec

   + [RB] matrix \f$\mbox{inv}(R_I^t * R_I) = \mbox{inv}(R_2^t*C_2^t*C_2*R_2)\f$, used for the rigid body modes projections

         rb_inv_RITRI.petscmat


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
