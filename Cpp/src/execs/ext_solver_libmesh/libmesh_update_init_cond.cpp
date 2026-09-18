#include "libmesh_update_init_cond.h"

/** \brief Program used to update the initial conditions for a libMesh \f$x_{\mbox{Min}}\f$ face.
 * 
 *  Usage: `./libmesh_update_init_cond__min_x_clamped -i [input file]`
 *  
 */

// Functions Prototype.  
// This function load the initial conditions from files and set them in the system
void set_state(libMesh::EquationSystems & es, 
               const std::string & system_name,
               const std::string & solution_file,
               const std::string & init_cond_disp_file,
               const std::string & init_cond_vel_file,
               const std::string & init_cond_acc_file);

int main(int argc, char** argv) {
  
  // --- Initialize libMesh
  libMesh::LibMeshInit init(argc, argv);

  // Do performance log?
  const bool MASTER_bPerfLog_carl_libmesh = true;
  libMesh::PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);

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

  libmesh_update_init_cond_input_params input_params;
  get_input_params(field_parser, input_params);

  // Check libMesh installation dimension
  const unsigned int dim = 3;

  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  libMesh::Mesh system_mesh(WorldComm, dim);
  system_mesh.read(input_params.space_mesh_file);
  system_mesh.prepare_for_use();

  // -------------------------------
  // --- Generate the equation systems
  perf_log.push("System setup:");

  // Set the equation systems object
  libMesh::EquationSystems equation_systems(system_mesh);

  // Add dynamic equation and physical parameters systems
  libMesh::NewmarkSystem& dynamic_system
                  = add_dynamic(equation_systems);

  // Initialize the equation systems  
  equation_systems.init();

  // Récupération du pas de temps
  std::ifstream file_time_steps(input_params.time_mesh_file);
  double init_time;
  double delta_t;
  file_time_steps >> init_time >> delta_t;

  // Récupération des paramètres du schéma de Newmark
  std::ifstream file_newmark(input_params.newmark_params_file);
  double newmark_alpha;
  double newmark_delta;
  file_newmark >> newmark_alpha >> newmark_delta;

  // Set the time step size, and optionally the
  // Newmark parameters, so that NewmarkSystem can
  // compute integration constants.  Here we simply use
  // pass only the time step and use default values
  // for alpha=.25  and delta=.5.
  dynamic_system.set_newmark_parameters(delta_t,newmark_alpha,newmark_delta);

  // Set the state of the system (initial conditions)
  set_state(equation_systems,
            "Dynamic",
            input_params.solution_file,
            input_params.init_cond_disp_file,
            input_params.init_cond_vel_file,
            input_params.init_cond_acc_file);
  
  perf_log.pop("System setup:");  
  // -------------------------------   

  // -------------------------------
  // Compute the new initial conditions
  perf_log.push("Computing the initial conditions");

  dynamic_system.update_u_v_a();

  perf_log.pop("Computing the initial conditions");
  // -------------------------------

  //---------------------------------------

 // Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  //dynamic_system.matrix->print_matlab(input_params.output_base_matrix + "_sys_mat.m");
#endif
  
  // Export the new initial conditions
  perf_log.push("Exporting the initial conditions");
  auto & disp = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(dynamic_system.get_vector("displacement"));
  auto & vel  = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(dynamic_system.get_vector("velocity"));
  auto & acc  = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(dynamic_system.get_vector("acceleration"));

  carl::write_PETSC_vector(disp, input_params.init_cond_disp_file );
  carl::write_PETSC_vector(vel,  input_params.init_cond_vel_file  );
  carl::write_PETSC_vector(acc,  input_params.init_cond_acc_file  );
  perf_log.pop("Exporting the initial conditions");
  
  return 0;
}

// This function load the initial conditions from files and set them in the system
void set_state(libMesh::EquationSystems & es, 
               const std::string & system_name,
               const std::string & solution_file,
               const std::string & init_cond_disp_file,
               const std::string & init_cond_vel_file,
               const std::string & init_cond_acc_file)
{
  // Get a reference to our system, as before
  libMesh::NewmarkSystem & t_system = es.get_system<libMesh::NewmarkSystem> (system_name);

  auto & sol  = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(*t_system.solution);
  auto & disp = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(t_system.get_vector("displacement"));
  auto & vel  = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(t_system.get_vector("velocity"));
  auto & acc  = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(t_system.get_vector("acceleration"));

  carl::read_PETSC_vector(sol,  solution_file);
  carl::read_PETSC_vector(disp, init_cond_disp_file);
  carl::read_PETSC_vector(vel,  init_cond_vel_file);
  carl::read_PETSC_vector(acc,  init_cond_acc_file);

}