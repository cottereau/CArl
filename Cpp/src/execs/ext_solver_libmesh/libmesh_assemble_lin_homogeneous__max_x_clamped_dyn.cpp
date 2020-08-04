// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1>FEMSystem Example 3 - Unsteady Linear Elasticity with
// FEMSystem</h1>
// \author Paul Bauman
// \date 2015
//
// This example shows how to solve the three-dimensional transient
// linear elasticity equations using the DifferentiableSystem class framework.
// This is just Systems of Equations Example 6 recast.


//Passing argument from file 
// ./libmesh_assemble_lin_homogeneous_min_x_clamped_dyn -i [input_params] -d [dynamic_params] -m [material_params]

// C++ includes
#include "libmesh_assemble_lin_homogeneous.h"

using namespace std;

#define x_scaling 1.3

// Bring in everything from the libMesh namespace
using namespace libMesh;

void set_clamped_boundary_info(libMesh::Mesh& mesh, 
                unsigned int bound_clamped_id, 
                unsigned int bound_force_id,
                bool node_enabeled);

// The main program.
int main (int argc, char ** argv)
{

  //[USER] Define boundary clamped
   const boundary_id_type bound_clamped_id = boundary_id_max_x;
  //[USER] Define boundary where force is defined
   const boundary_id_type bound_force_id   = boundary_id_max_z;
  //[USER] Define [N]
   Real force = 10000.; 

  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Do performance log?
  const bool MASTER_bPerfLog_carl_libmesh = true;
  libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);
  
  // Communicator MPI
  libMesh::Parallel::Communicator& WorldComm = init.comm();

  // Number of processors and processor rank.
  int rank = WorldComm.rank();
  int nodes = WorldComm.size();

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // This example requires 3D calculations
  libmesh_example_requires(LIBMESH_DIM > 2, "3D support");

  // We use Dirichlet boundary conditions here
  #ifndef LIBMESH_ENABLE_DIRICHLET
    libmesh_example_requires(false, "--enable-dirichlet");
  #endif

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

  // Declaration of a struct variable to catch all inputs from the input file
  libmesh_assemble_input_params input_params; // struct type 

  //Getting parsing argument from all file
  get_input_params(field_parser, input_params);

  // Initialize the cantilever mesh
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // Create a 3D mesh distributed across the default MPI communicator.
  perf_log.push("Meshes - Parallel","Read files:");
  libMesh::Mesh mesh(WorldComm, dim);
  mesh.read(input_params.mesh_file);
  mesh.prepare_for_use();

  perf_log.pop("Meshes - Parallel","Read files:");

  // --- Generate the equation systems
  perf_log.push("System setup:");;

  // Print information about the mesh on the screen
  mesh.print_info();

  // Let's add some node and edge boundary conditions.
  // Each processor should know about each boundary condition it can
  // see, so we loop over all elements, not just local elements.

  /*for (const auto & elem : mesh.element_ptr_range())
    {
      unsigned int side_min_x = 0;
      unsigned int side_max_x = 0;
      unsigned int side_max_z = 0;
      bool found_side_min_x = false;
      bool found_side_max_x = false;
      bool found_side_max_z = false;
      for (auto side : elem->side_index_range())
        {
          // if (mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_max_x))
          //   {
          //     side_max_x = side;
          //     found_side_max_x = true;
          //   }
          if(mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_min_x))
            {
              side_min_x = side;
              found_side_min_x = true;
            }
          if(mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_max_z))
            {
              side_max_z = side;
              found_side_max_z = true;
            }

         /* if (mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_min_y))
            {
              side_min_y = side;
              found_side_min_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_max_y))
            {
              side_max_y = side;
              found_side_max_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, boundary_id_max_z))
            {
              side_max_z = side;
              found_side_max_z = true;
            }
        }
        if(found_side_min_x)
          for(auto n:elem->node_index_range())
              if(elem->is_node_on_side(n, side_min_x))
                   mesh.get_boundary_info().add_node(elem->node_ptr(n), node_boundary_id);
        if(found_side_max_z)
          for(auto n:elem->node_index_range())
              if(elem->is_node_on_side(n, side_max_z))
                   mesh.get_boundary_info().add_node(elem->node_ptr(n), pressure_boundary_id);

      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X, BOUNDARY_ID_MAX_Y, BOUNDARY_ID_MAX_Z

        if(found_side_min_x)
          for(auto e : elem->edge_index_range())
            if(elem->is_edge_on_side(e, side_min_x))
              mesh.get_boundary_info().add_edge(elem, e,edge_boundary_id);
        // if(found_side_min_x)
        //   for(auto e: elem->edge_index_range())
        //     if(elem->is_edge_on_side(e,side_min_x))
        //       mesh.get_boundary_info().add_edge(elem,e,edge_boundary_id);
        if(found_side_max_z)
          for(auto e: elem->edge_index_range())
            if(elem->is_edge_on_side(e,side_max_z))
              mesh.get_boundary_info().add_edge(elem,e,pressure_boundary_id);

                  // then let's set a node boundary condition
      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X and BOUNDARY_ID_MIN_Y
      // then let's set an edge boundary condition
    }*/

  // Create an equation systems object
  EquationSystems equation_systems (mesh);

  // Declare the system "Navier-Stokes" and its variables.
  ElasticitySystem & system =
    equation_systems.add_system<ElasticitySystem> ("Linear Elasticity");

  //passing force to system
  system.setForce(force);

  set_clamped_boundary_info(
                mesh, 
                bound_clamped_id, 
                bound_force_id,
                true);


  // Reading of material parameter
  system.read_material_params(input_params.physical_params_file);

  // Solve this as a time-dependent or steady system
  std::string time_solver = std::string("newmark")/*infile("time_solver","DIE!")*/;

  ExplicitSystem * v_system;
  ExplicitSystem * a_system;

  if( time_solver == std::string("newmark") )
    {
      // Create ExplicitSystem to help output velocity
      v_system = &equation_systems.add_system<ExplicitSystem> ("Velocity");
      v_system->add_variable("u_vel", FIRST, LAGRANGE);
      v_system->add_variable("v_vel", FIRST, LAGRANGE);
      v_system->add_variable("w_vel", FIRST, LAGRANGE);

      // Create ExplicitSystem to help output acceleration
      a_system = &equation_systems.add_system<ExplicitSystem> ("Acceleration");
      a_system->add_variable("u_accel", FIRST, LAGRANGE);
      a_system->add_variable("v_accel", FIRST, LAGRANGE);
      a_system->add_variable("w_accel", FIRST, LAGRANGE);
    }

  if (time_solver == std::string("newmark"))
    {
      system.time_solver = libmesh_make_unique<NewmarkSolver>(system);
    }
  else
    libmesh_error_msg(std::string("ERROR: invalid time_solver ")+time_solver);

  // Initialize the system
  //->launch init_data() methode of ElasticitySystem
  equation_systems.init ();

  // Set the time stepping options 
  // Read step time
  const Real deltat           = static_cast<Real>(input_params.deltat);
  system.deltat = deltat;

  // And the nonlinear solver options
  DiffSolver & solver = *(system.time_solver->diff_solver().get());
  solver.quiet = input_params.solver_quiet;
  solver.verbose = !solver.quiet;
  solver.max_nonlinear_iterations = input_params.max_nonlinear_iterations;
  solver.relative_step_tolerance = static_cast<Real>(input_params.relative_step_tolerance);
  solver.relative_residual_tolerance = static_cast<Real>(input_params.relative_residual_tolerance);
  solver.absolute_residual_tolerance = static_cast<Real>(input_params.absolute_residual_tolerance);

  // And the linear solver options
  solver.max_linear_iterations = input_params.max_linear_iterations;
  solver.initial_linear_tolerance = static_cast<Real>(input_params.initial_linear_tolerance);

  // Print information about the system to the screen.
  equation_systems.print_info();

  if( time_solver == std::string("newmark") )
    {
      NewmarkSolver * newmark = cast_ptr<NewmarkSolver*>(system.time_solver.get());
      newmark->compute_initial_accel();

      // Copy over initial velocity and acceleration for output.
      // Note we can do this because of the matching variables/FE spaces
      *(v_system->solution) = system.get_vector("_old_solution_rate");
      *(a_system->solution) = system.get_vector("_old_solution_accel");
    }

#ifdef LIBMESH_HAVE_EXODUS_API
  // Output initial state
  {
    std::ostringstream file_name;

    // We write the file in the ExodusII format.
    file_name << std::string("out_mesh_max_.")+time_solver+std::string(".e-s.")
              << std::setw(3)
              << std::setfill('0')
              << std::right
              << 0;

    ExodusII_IO(mesh).write_timestep(file_name.str(),
                                     equation_systems,
                                     1, // This number indicates how many time steps
                                        // are being written to the file
                                     system.time);
  }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API


  // Now we begin the timestep loop to compute the time-accurate
  // solution of the equations.
   // Read number of step
  unsigned int n_timesteps    = input_params.n_timesteps;
  for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
    {
      // A pretty update message
      libMesh::out << "\n\nSolving time step "
                   << t_step
                   << ", time = "
                   << system.time
                   << std::endl;

      system.solve(); //-> lance la fonction ElasticitySystem::init_data
      // Advance to the next timestep in a transient problem
      system.time_solver->advance_timestep();//->launch ElasticitySystem::side_time_derivative
                                                //->launch  ElasticitySystem::element_time_derivative                

      // Copy over updated velocity and acceleration for output.
      // Note we can do this because of the matching variables/FE spaces


      if( time_solver == std::string("newmark") )
        {
          *(v_system->solution) = system.get_vector("_old_solution_rate");
          *(a_system->solution) = system.get_vector("_old_solution_accel");
        }

#ifdef LIBMESH_HAVE_EXODUS_API
    const unsigned int write_interval  = input_params.write_interval;
  #endif

#ifdef LIBMESH_HAVE_EXODUS_API
      // Write out this timestep if we're requested to
      if ((t_step+1)%write_interval == 0)
        {
          std::ostringstream file_name;

          // We write the file in the ExodusII format.
          file_name << std::string("out.")+time_solver+std::string(".e-s.")
                    << std::setw(3)
                    << std::setfill('0')
                    << std::right
                    << t_step+1;

          ExodusII_IO(mesh).write_timestep(file_name.str(),
                                           equation_systems,
                                           1, // This number indicates how many time steps
                                              // are being written to the file
                                           system.time);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }
  // Export matrix and vector
  libMesh::PetscMatrix<libMesh::Number> * temp_mat_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(system.matrix);
  libMesh::PetscVector<libMesh::Number> * temp_vec_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> * >(system.rhs);

  carl::write_PETSC_matrix(*temp_mat_ptr, input_params.output_base + "_sys_mat.petscmat");
  carl::write_PETSC_vector(*temp_vec_ptr, input_params.output_base + "_sys_rhs_vec.petscvec");

  // If needed, print rigid body vectors
  if(input_params.bCalculateRBVectors)
  {
    MatNullSpace nullsp_sys;
    build_rigid_body_vectors(system,nullsp_sys);
    write_rigid_body_vectors(nullsp_sys,input_params.output_base,WorldComm.rank());
    MatNullSpaceDestroy(&nullsp_sys);
  }
  // All done.
  //
  return 0;
}


void set_clamped_boundary_info(
                libMesh::Mesh& mesh, 
                unsigned int bound_clamped_id, 
                unsigned int bound_force_id,
                bool node_enabeled = false)

{
  unsigned int n_boundaries = 6;
    for (const auto & elem : mesh.element_ptr_range())
    {
      std::vector<unsigned int> list_side (6,0);
      std::vector<bool> list_found(6,0);

      for (auto side : elem->side_index_range())
        {          
          for(unsigned int it; it<n_boundaries; it++)
            if(mesh.get_boundary_info().has_boundary_id(elem, side, it))
              {
                list_side[it] = side;
                list_found[it] = true;
              }
        }

        if(node_enabeled)
        {
          if(list_found[bound_clamped_id])
            for(auto n:elem->node_index_range())
                if(elem->is_node_on_side(n, list_side[bound_clamped_id]))
                     mesh.get_boundary_info().add_node(elem->node_ptr(n), node_boundary_id);

          if(list_found[bound_force_id])
            for(auto n:elem->node_index_range())
                if(elem->is_node_on_side(n, list_side[bound_force_id]))
                     mesh.get_boundary_info().add_node(elem->node_ptr(n), pressure_boundary_id);
        }
        else {

          if(list_found[bound_clamped_id])
            for(auto e : elem->edge_index_range())
              if(elem->is_edge_on_side(e, list_side[bound_clamped_id]))
                mesh.get_boundary_info().add_edge(elem, e,edge_boundary_id);

          if(list_found[bound_force_id])
            for(auto e: elem->edge_index_range())
              if(elem->is_edge_on_side(e,list_side[bound_clamped_id]))
              mesh.get_boundary_info().add_edge(elem,e,pressure_boundary_id);
        }
    }
}
/* Local Variables:                                                        */
/* mode: c                                                              */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=4 et tw=80 smartindent :                               */
