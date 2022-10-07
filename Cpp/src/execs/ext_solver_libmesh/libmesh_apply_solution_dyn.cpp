/*
 * \file libmesh_apply_solution_dyn.cpp
 *
 *  Created on: Feb 23, 2022
 *      Author: Chensheng Luo
 *
 * \brief  **DYN** Program responsible to transfert all petscvec result (in sense of D.O.F.) to EXODUS II mesh vector(.e)
 * 
 *  This program's input file description can be found at the documentation of the function carl::get_input_params(GetPot& field_parser, libmesh_apply_solution_dyn_input_params& input_params). 
 */

#include "libmesh_apply_solution_dyn.h"


int main(int argc, char** argv) {

  // --- Initialize libMesh
  libMesh::LibMeshInit init(argc, argv);

  // Do performance log?
  libMesh::PerfLog perf_log("Main program");

  // libMesh's C++ / MPI communicator wrapper
  libMesh::Parallel::Communicator& WorldComm = init.comm();

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

  carl::libmesh_apply_solution_dyn_input_params input_params;
  carl::get_input_params(field_parser, input_params);

  if(WorldComm.rank() == 0)
  {
    std::string command_string;

    command_string = "rm -rf " + input_params.output_mesh_folder;
    carl::exec_command(command_string.c_str());

    command_string = "mkdir -p " + input_params.output_mesh_folder;
    carl::exec_command(command_string.c_str());
    std::cout << command_string << std::endl;
  }


  // Check libMesh installation dimension
  const unsigned int dim = 3;

  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // - Parallelized meshes A
  libMesh::Mesh system_mesh(WorldComm, dim);
  system_mesh.read(input_params.input_mesh);
  system_mesh.prepare_for_use();

  // Set the equation systems object
  libMesh::EquationSystems equation_systems(system_mesh);

  // Add linear elasticity and stress
  libMesh::LinearImplicitSystem& elasticity_system
                    = add_elasticity(equation_systems,"Elasticity");
  add_stress(equation_systems);

  // Initialize the equation systems
  equation_systems.init();

  // Homogeneous physical properties
  set_homogeneous_physical_properties(equation_systems, input_params.physical_params_file);

  // Read the solution vector
  for(int i=input_params.step_loop_times;i<=input_params.total_loop_times;i=i+input_params.step_loop_times){
    libMesh::PetscVector<libMesh::Real> * sol_vec_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Real> * >(elasticity_system.solution.get());
    carl::read_PETSC_vector(sol_vec_ptr->vec(),input_params.input_vector_folder+std::to_string(i)+".petscvec", WorldComm.get());

    // Close it and update
    elasticity_system.solution->close();
    elasticity_system.update();

    // Calculate the stress
    compute_stresses(equation_systems);

    // Export solution
    // #ifdef LIBMESH_HAVE_EXODUS_API
    //   libMesh::VTKIO vtkio_interface(system_mesh);

    //   std::set<std::string> system_names;
    //   system_names.insert("Elasticity");
    //   vtkio_interface.write_equation_systems(input_params.output_mesh_folder+"."+std::to_string(i)+".vtk",equation_systems,&system_names);
    //   //vtkio_interface.write_element_data(equation_systems);
    // #endif

    #ifdef LIBMESH_HAVE_EXODUS_API
      libMesh::ExodusII_IO exo_io_interface(system_mesh, /*single_precision=*/true);

      std::set<std::string> system_names;
      system_names.insert("Elasticity");
      
      // std::string str = std::to_string(i);

      // int n=0;

      // int total=input_params.step_loop_times;
      // while(total != 0)
      // {
      //   total /= 10;
      //   ++n;
      // }
 
      // int precision = n - str.length();
      // std::string s = std::string(precision, '0').append(str);
      std::string filename=input_params.output_mesh_folder+"sol.e-s."+std::to_string(i);
      exo_io_interface.write_equation_systems(filename,equation_systems,&system_names);
      exo_io_interface.write_timestep(filename,equation_systems,i,input_params.Newmark.deltat*i,&system_names);
      exo_io_interface.write_element_data(equation_systems);
    #endif
  }

  return 0;
}
