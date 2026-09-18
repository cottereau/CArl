#include "libmesh_assemble_rhs_dyn.h"

/** \brief Program used to assemble the right-hand side for a dynamic system with a clamped \f$x_{\mbox{Min}}\f$ face.
 * 
 *  Usage: `./libmesh_assemble_rhs_dyn__min_x_clamped -i [input file]`
 *  
 */

 // Functions Prototype.  
 // This function load the initial conditions from files and set them in the system
void set_state(libMesh::EquationSystems & es, 
               const std::string & system_name,
               const std::string & init_cond_disp_file,
               const std::string & init_cond_vel_file,
               const std::string & init_cond_acc_file);
// This function loads the system matrices from files.
void load_matrices(libMesh::EquationSystems & es,
                   const std::string & system_name,
                   const std::string & mass_file,
                   const std::string & damping_file,
                   const std::string & stiffness_file,
                   const std::string & system_matrix_file);
// This function applies the Dirichlet boundary conditions
void fill_rhs_dirichlet_bc(libMesh::EquationSystems & es,
                             const std::string & system_name);                   


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

  libmesh_assemble_rhs_dyn_input_params input_params;
  get_input_params(field_parser, input_params);

  // Check libMesh installation dimension
  const unsigned int dim = 3;

  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // --- Declare the three meshes to be intersected

  // - Parallelized meshes: A, B, mediator and weight
  perf_log.push("Meshes - Parallel","Read files:");
  libMesh::Mesh system_mesh(WorldComm, dim);
  system_mesh.read(input_params.space_mesh_file);
  system_mesh.prepare_for_use();

  perf_log.pop("Meshes - Parallel","Read files:");

  // --- Generate the equation systems
  perf_log.push("System setup:");

  // Set the equation systems object
  libMesh::EquationSystems equation_systems(system_mesh);
  
  //---------------------------------------
  // Ma partie adaptée de l'exemple de dynamique 

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
            input_params.scratch_folder_path + "/domain_A_state/init_cond_A_disp.petscvec",
            input_params.scratch_folder_path + "/domain_A_state/init_cond_A_vel.petscvec",
            input_params.scratch_folder_path + "/domain_A_state/init_cond_A_acc.petscvec");

  // Load the system matrices from files
  load_matrices(equation_systems,
                "Dynamic",
                input_params.domain_M_matrix_file,
                input_params.domain_C_matrix_file,
                input_params.domain_K_matrix_file,
                input_params.domain_matrix_file);

  
  // Set the system time
  // Récupération de l'indice du dernier instant pour lequel l'état du système est connu
  std::string path_string;
  path_string = input_params.scratch_folder_path + "/Time_Mono_state.dat";
  std::ifstream file_state(path_string);
  double index_time;;
  file_state >> index_time;

  dynamic_system.time = init_time + index_time * delta_t;           

  perf_log.pop("System setup:");     

  // Assembling the right-hand side vector
  perf_log.push("Assembling the right-hand side vector");

  dynamic_system.time += delta_t;
  dynamic_system.update_rhs();

  fill_rhs_dirichlet_bc(equation_systems,"Dynamic");

  perf_log.pop("Assembling the right-hand side vector");


  //---------------------------------------

 // Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  //dynamic_system.matrix->print_matlab(input_params.output_base_matrix + "_sys_mat.m");
#endif
  
  // Export RHS
  perf_log.push("Exporting the right-hand side vector");
  libMesh::PetscVector<libMesh::Number> * RHS_ptr =
  libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(dynamic_system.rhs);

  carl::write_PETSC_vector(*RHS_ptr,  input_params.scratch_folder_path + "/domain_A_state/traction_model_A_sys_rhs_vec.petscvec");
  perf_log.pop("Exporting the right-hand side vector");
  
  return 0;
}

// This function load the initial conditions from files and set them in the system
void set_state(libMesh::EquationSystems & es, 
               const std::string & system_name,
               const std::string & init_cond_disp_file,
               const std::string & init_cond_vel_file,
               const std::string & init_cond_acc_file)
{
  // Get a reference to our system, as before
  libMesh::NewmarkSystem & t_system = es.get_system<libMesh::NewmarkSystem> (system_name);


  auto & disp = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(t_system.get_vector("displacement"));
  auto & vel  = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(t_system.get_vector("velocity"));
  auto & acc  = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number>&>(t_system.get_vector("acceleration"));

  carl::read_PETSC_vector(disp, init_cond_disp_file);
  carl::read_PETSC_vector(vel,  init_cond_vel_file);
  carl::read_PETSC_vector(acc,  init_cond_acc_file);

}

// This function loads the system matrices from files.
void load_matrices(libMesh::EquationSystems & es,
                   const std::string & system_name,
                   const std::string & mass_file,
                   const std::string & damping_file,
                   const std::string & stiffness_file,
                   const std::string & system_matrix_file)
{
    // Get a reference to the Newmark system
    libMesh::NewmarkSystem & t_system = es.get_system<libMesh::NewmarkSystem>(system_name);

    auto & mass          = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>&>(t_system.get_matrix("mass"));
    auto & damping       = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>&>(t_system.get_matrix("damping"));
    auto & stiffness     = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>&>(t_system.get_matrix("stiffness"));
    auto & system_matrix = libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number>&>(*t_system.matrix);

    carl::read_PETSC_matrix(mass,          mass_file);
    carl::read_PETSC_matrix(damping,       damping_file);
    carl::read_PETSC_matrix(stiffness,     stiffness_file);
    carl::read_PETSC_matrix(system_matrix, system_matrix_file);

}

// This function applies the Dirichlet boundary conditions
void fill_rhs_dirichlet_bc(libMesh::EquationSystems & es,
                             const std::string & system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Dynamic");

  // Get a reference to our system, as before.
  libMesh::NewmarkSystem & t_system = es.get_system<libMesh::NewmarkSystem> (system_name);

  // Numéro du système dans l'equation_systems
  const unsigned int sys_num = t_system.number();

  // Get writable references to the overall matrix and vector.
  libMesh::NumericVector<libMesh::Number> & rhs    = *t_system.rhs;

  // Get a constant reference to the mesh object.
  const libMesh::MeshBase & mesh = es.get_mesh();

  // Get libMesh's pi
  const libMesh::Real pi = libMesh::pi;

  // Number of nodes in the mesh.
  unsigned int n_nodes = mesh.n_nodes();
  
  for (unsigned int n_cnt=0; n_cnt<n_nodes; n_cnt++)
    {
    
      // Get a reference to the current node.
      const libMesh::Node & curr_node = mesh.node_ref(n_cnt);
      
      // Check if Dirichlet BCs should be applied to this node.
      // Use the TOLERANCE from mesh_common.h as tolerance.
      // Here a pressure value is applied if the z-coord.
      // is equal to 4, which corresponds to one end of the
      // pipe-mesh in this directory.

      // The penalty parameter.
      const libMesh::Real penalty = 1.e10;

      if (std::abs(curr_node(2)) < 1.e-6)
        {

          // Boucle sur les ddl du noeud
          for (unsigned int d=0; d<3; d++)
          {
            // The global number of the respective degree of freedom.
            unsigned int dn = curr_node.dof_number(sys_num, d, 0);

            // Here we apply null displactement throug u, v, w
            // at one end of the barre-mesh.
            libMesh::Real disp_value = 0.0;

            // Now add the contributions to the matrix and the rhs.
            rhs.add(dn, disp_value*penalty);

          }
          
        }
    } // loop n_cnt

    rhs.close();

}