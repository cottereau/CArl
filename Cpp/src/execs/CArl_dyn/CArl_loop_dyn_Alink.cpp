#include "CArl_loop_dyn.h"

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

  // --- Set up input_params

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

  carl::feti_loop_dyn_params input_params;
  get_input_params(field_parser, input_params);

  carl::FETI_Dyn_Operations feti_op(WorldComm,input_params.scratch_folder_path);
  
  // Calculate B_link_speed/displacement
  feti_op.Newmark_speed_link(input_params.gammaA,
      input_params.deltatA,
      input_params.scratch_folder_path+"/this_acc_A_link_sys_sol_vec.petscvec",
      input_params.scratch_folder_path+"/this_speed_A_link.petscvec");

  feti_op.Newmark_displacement_link(input_params.betaA,
      input_params.deltatA,
      input_params.scratch_folder_path+"/this_acc_A_link_sys_sol_vec.petscvec",
      input_params.scratch_folder_path+"/this_disp_A_link.petscvec");

  // Calculate B_speed/displacement
  feti_op.add_free_link(input_params.scratch_folder_path+"/this_acc_A_free_sys_sol_vec.petscvec",
  		input_params.scratch_folder_path+"/this_acc_A_link_sys_sol_vec.petscvec",
      	input_params.scratch_folder_path+"/this_acc_A.petscvec");
  feti_op.add_free_link(input_params.scratch_folder_path+"/this_speed_A_free.petscvec",
  		input_params.scratch_folder_path+"/this_speed_A_link.petscvec",
      	input_params.scratch_folder_path+"/this_speed_A.petscvec");
  feti_op.add_free_link(input_params.scratch_folder_path+"/this_disp_A_free.petscvec",
  		input_params.scratch_folder_path+"/this_disp_A_link.petscvec",
      	input_params.scratch_folder_path+"/this_disp_A.petscvec");

// - Judge progression
  std::string progression_filename;
  progression_filename = input_params.scratch_folder_path+"/iteration_progression.txt";
  field_parser.parse_input_file(progression_filename, "#", "\n", " \t\n");
  carl::feti_loop_dyn_iteration_progression_params progression_params;
  get_input_params(field_parser, progression_params);

  //Output result file and change all "this" to "prev"
  int index=progression_params.inner_loop_progression+progression_params.outer_loop_progression*input_params.inner_loop_times+1;
  feti_op.output_A_result(input_params.result_folder_path,progression_params.inner_loop_progression+progression_params.outer_loop_progression*input_params.inner_loop_times+1);
  feti_op.move_to_prev_A();
  feti_op.prepare_rhs_vector(input_params.betaA,
      input_params.deltatA,
      input_params.scratch_folder_path+"/force_A/force_"+std::to_string(index)+".petscvec",
      input_params.scratch_folder_path+"/this_acc_A.petscvec",
      input_params.scratch_folder_path+"/this_speed_A.petscvec",
      input_params.scratch_folder_path+"/this_disp_A.petscvec",
      input_params.stiffness_matrix_A,
      input_params.scratch_folder_path+"/rhs_vec_A_free.petscvec");

  if (progression_params.outer_loop_progression < (input_params.outer_loop_times-1)){

    //-countinue to go in outer loop!

      //Rewrite progression with outer lopp +1, inner loop go back to 0
        progression_params.outer_loop_progression = progression_params.outer_loop_progression + 1;
        progression_params.inner_loop_progression = 0;
        print_input_params(progression_filename,progression_params);

      if(WorldComm.rank() == 0){
         std::string init_script_command = ". " + input_params.scratch_folder_path + "/Bfree.sh";
         if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
         {
            std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
            std::cout << init_script_command << std::endl ;
         } else {
            carl::exec_command(init_script_command);
         }
        }
    }else{
        //Outer loop end!
       if(WorldComm.rank() == 0){
          std::cout << " The whole iteration ends!" << std::endl;
        }
    }

  
  return 0;
}