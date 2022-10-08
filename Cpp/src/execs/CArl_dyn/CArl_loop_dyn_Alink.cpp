/*
 *
 *  Created on: Nov 17，2021
 *      Author: Chensheng Luo
 */
#include "CArl_loop_dyn.h"

/** \file CArl_loop_dyn_Alink.cpp

\brief **DYN-DI/DYN-CG** Program responsible to calculate A link speed and displacement by Newmark method, and add up to get final speed/displacement.

This program's input file description can be found at the documentation of the function carl::get_input_params(GetPot& field_parser, feti_loop_dyn_params& input_params). 

In this step, it will use the following files in scratch foler:
     - `this_acc_A_link_sys_sol_vec.petscvec` : \f$ \ddot{U}^A_{link}(t) \f$
     - `this_acc_A_free_sys_sol_vec.petscvec` : \f$ \ddot{U}^A_{free}(t) \f$

It will create:
     - `this_speed_A_link.petscvec`: \f$ \dot{U}^A_{link}(t) \f$
     - `this_disp_A_link.petscvec`: \f$ U^A_{link}(t) \f$
     - `this_acc_A.petscvec` : \f$ \ddot{U}^A(t) \f$
     - `this_speed_A.petscvec`: \f$ \dot{U}^A(t) \f$
     - `this_disp_A.petscvec`: \f$ U^A(t) \f$

At the end, it will apply a test, move all `this` to `prev` for object A, rewrite `iteration_progression.txt` and set all A/interpolation related vector to NAN.

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

  carl::FETI_Dyn_Operations Dyn_op(WorldComm,input_params.scratch_folder_path,input_params.result_folder_path);
  
  // Calculate A_link_speed/displacement
  Dyn_op.Newmark_speed_link(&input_params.newmark_A,
      &Dyn_op.vector_A);

  Dyn_op.Newmark_displacement_link(&input_params.newmark_A,
      &Dyn_op.vector_A);

  Dyn_op.add_free_link(&Dyn_op.vector_A);

// - Judge progression

  std::string progression_filename;
  progression_filename = input_params.scratch_folder_path+"/iteration_progression.txt";
  field_parser.parse_input_file(progression_filename, "#", "\n", " \t\n");
  carl::feti_loop_dyn_iteration_progression_params progression_params;
  get_input_params(field_parser, progression_params);
  

  //Output result file and change all "this" to "prev"
  int index=progression_params.inner_loop_progression+progression_params.outer_loop_progression*input_params.inner_loop_times;
  if(index%(input_params.result_times_A*input_params.inner_loop_times)){
    std::cout << "According to ResultTimeA, no output for this time step!" << std::endl;
  }else{
    Dyn_op.output_A_result(index);
  }
  Dyn_op.move_to_prev_A();
  #ifdef PRINT_CALCULATION_TIME
    Dyn_op.export_calculation_time(index,"Alink");
  #endif

  //Test coupling 
  #ifdef TEST_COUPLING_RESULT
    if(index%(input_params.result_times_A*input_params.inner_loop_times)+index%input_params.result_times_B){
        std::cout << "As one output is omitted for this index, no test could be done!" << std::endl;
    }else{
        Dyn_op.test_coupling(input_params.matrix_A.coupling,input_params.matrix_B.coupling,index);
        #ifdef PRINT_CALCULATION_TIME
          Dyn_op.export_calculation_time(index,"Test");
        #endif
    }
  #endif
  
  //delete all value of this moment
  Dyn_op.delete_A_this_vector(input_params.dyn_solver);
  Dyn_op.delete_coupling_vector(input_params.dyn_solver);
//Prepare next rhs vector
  Dyn_op.prepare_A_next_force(&Dyn_op.vector_A,
    index,
    input_params.inner_loop_times,
    input_params.force_prepare_method,
    input_params.force_prepare_params,
    input_params.newmark_B.deltat);
  Dyn_op.rhs_free(&input_params.newmark_A,
        &Dyn_op.vector_A,
        &input_params.matrix_A);

    //Rewrite progression with inner loop go back to 0
    progression_params.outer_loop_progression = progression_params.outer_loop_progression + 1;
    progression_params.inner_loop_progression = 0;
    print_input_params(progression_filename,progression_params);

  if (progression_params.outer_loop_progression < input_params.outer_loop_times){

    //-countinue to go in outer loop!


      if(WorldComm.rank() == 0){
        
         std::string init_script_command = ". " + input_params.scratch_folder_path + "/Afree.sh";
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