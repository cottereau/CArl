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
  
  // - Calculate B_link_speed/displacement
  feti_op.Newmark_speed_link(input_params.gammaB,
      input_params.deltatB,
      input_params.scratch_folder_path+"/this_acc_B_link_sys_sol_vec.petscvec",
      input_params.scratch_folder_path+"/this_speed_B_link.petscvec");

  feti_op.Newmark_displacement_link(input_params.betaB,
      input_params.deltatB,
      input_params.scratch_folder_path+"/this_acc_B_link_sys_sol_vec.petscvec",
      input_params.scratch_folder_path+"/this_disp_B_link.petscvec");

  // - Calculate B_speed/displacement
  feti_op.add_free_link(input_params.scratch_folder_path+"/this_acc_B_free_sys_sol_vec.petscvec",
  		input_params.scratch_folder_path+"/this_acc_B_link_sys_sol_vec.petscvec",
      	input_params.scratch_folder_path+"/this_acc_B.petscvec");
  feti_op.add_free_link(input_params.scratch_folder_path+"/this_speed_B_free.petscvec",
  		input_params.scratch_folder_path+"/this_speed_B_link.petscvec",
      	input_params.scratch_folder_path+"/this_speed_B.petscvec");
  feti_op.add_free_link(input_params.scratch_folder_path+"/this_disp_B_free.petscvec",
  		input_params.scratch_folder_path+"/this_disp_B_link.petscvec",
      	input_params.scratch_folder_path+"/this_disp_B.petscvec");

  // - Judge progression
  std::string progression_filename;
  progression_filename = input_params.scratch_folder_path+"/iteration_progression.txt";
  field_parser.parse_input_file(progression_filename, "#", "\n", " \t\n");
  carl::feti_loop_dyn_iteration_progression_params progression_params;
  get_input_params(field_parser, progression_params);

  //Output result file and change all "this" to "prev"
  int index=progression_params.inner_loop_progression+progression_params.outer_loop_progression*input_params.inner_loop_times+1;
  feti_op.output_B_result(input_params.result_folder_path,index);
  feti_op.move_to_prev_B();
  feti_op.prepare_rhs_vector(input_params.betaB,
      input_params.deltatB,
      input_params.scratch_folder_path+"/force_B/force_"+std::to_string(index)+".petscvec",
      input_params.scratch_folder_path+"/this_acc_B.petscvec",
      input_params.scratch_folder_path+"/this_speed_B.petscvec",
      input_params.scratch_folder_path+"/this_disp_B.petscvec",
      input_params.stiffness_matrix_B,
      input_params.scratch_folder_path+"/rhs_vec_B_free.petscvec");

  if (progression_params.inner_loop_progression < (input_params.inner_loop_times-1)){

    //-countinue to go in inner loop!

      //Rewrite progression with inner lopp +1
        progression_params.inner_loop_progression = progression_params.inner_loop_progression + 1;
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

        //Inner loop end!

       if(WorldComm.rank() == 0){
         std::string init_script_command = ". " + input_params.scratch_folder_path + "/Alink.sh";
         if(input_params.scheduler == carl::ClusterSchedulerType::LOCAL)
         {
            std::cout << " !!! LOCAL job 'scheduler: Run the following script manually: " << std::endl;
            std::cout << init_script_command << std::endl ;
         } else {
            carl::exec_command(init_script_command);
         }
        }

    }
  
  return 0;
}