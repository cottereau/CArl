#include "libmesh_init_dyn_lin_homogeneous.h"

/** \brief Program used to initialize the dynamic system for a homogeneous model with a traction applied to \f$x_{\mbox{Max}}\f$ face.
 * 
 *  Usage: `./libmesh_init_dyn_lin_homogeneous__max_x_traction -i [input file]`
 *  
 * The input file is parsed by the get_input_params(GetPot& field_parser, libmesh_assemble_input_params& input_params) function, and it contains the following parameters. 
 *
 *  Required parameters:
 *    - `Mesh` : path to the mesh.
 *    - `PhysicalParameters` : physical parameters.
 *    - `SystemType` : parameter used to tell the assembler which weight functions must be used. *Values*: `Micro` or `Macro`.
 *    - `MeshWeight` : path to the mesh defining the domains of the Arlequin weight parameters.
 *    - `WeightIndexes` : path to the indices of the domains of the Arlequin weight parameters.
 *
 *  Optional parameter:
 *    - `OutputBase` or `--output` : base of the output files (including folders). *Default*: `test_system`.
 *
 *  Boolean flags:
 *    - `ExportRBVectors` : build and export the rigid body modes vectors.
 */

 // Function Prototype.  This function will be used to apply the
// initial conditions.
void apply_initial(libMesh::EquationSystems & es,
                   const std::string & system_name);

// Function Prototype.  This function imposes
// Dirichlet Boundary conditions via the penalty
// method after the system is assembled.
void fill_mat_dirichlet_bc(libMesh::EquationSystems & es,
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

  libmesh_init_dyn_input_params input_params;
  get_input_params(field_parser, input_params);

// Check libMesh installation dimension
  const unsigned int dim = 3;

  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // --- Declare the three meshes to be intersected

  // - Parallelized meshes: A, B, mediator and weight
  perf_log.push("Meshes - Parallel","Read files:");
  libMesh::Mesh system_mesh(WorldComm, dim);
  system_mesh.read(input_params.mesh_file);
  system_mesh.prepare_for_use();

  libMesh::Mesh mesh_weight(WorldComm, dim);
  mesh_weight.allow_renumbering(false);
  mesh_weight.read(input_params.mesh_weight_file);
  mesh_weight.prepare_for_use();

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
  std::ifstream file_time_steps(input_params.time_discretization_file);
  double dummy;
  double delta_t;
  file_time_steps >> dummy >> delta_t;

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

  // Homogeneous properties for the macro system
  set_homogeneous_physical_properties(equation_systems, input_params.physical_params_file);

  // Set the weight function object
  weight_parameter_function  system_weight(mesh_weight);
  system_weight.set_parameters(input_params.weight_domain_idx_file);

  perf_log.pop("System setup:");

  // Assemble M, C and K
  assemble_matrices_dynamic_with_weight(equation_systems,"Dynamic",system_weight,
              input_params.system_type);
  
  // Compute K~=K+a_0*M+a_1*C 
  dynamic_system.compute_matrix();

  // Ajout de la pénalisation pour les deux domaines
  fill_mat_dirichlet_bc(equation_systems, "Dynamic");

  //---------------------------------------

 // Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  dynamic_system.matrix->print_matlab(input_params.output_base_matrix + "_sys_mat.m");
#endif

  // Export initial conditions
  libMesh::PetscVector<libMesh::Number> * U_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&dynamic_system.get_vector("displacement"));
  libMesh::PetscVector<libMesh::Number> * V_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&dynamic_system.get_vector("velocity"));
  libMesh::PetscVector<libMesh::Number> * A_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&dynamic_system.get_vector("acceleration"));

  carl::write_PETSC_vector(*U_ptr, input_params.output_base_init_cond + "_disp.petscvec");
  carl::write_PETSC_vector(*V_ptr, input_params.output_base_init_cond + "_vel.petscvec");
  carl::write_PETSC_vector(*A_ptr, input_params.output_base_init_cond + "_acc.petscvec");

  // Export matrix
  libMesh::PetscMatrix<libMesh::Number> * temp_mat_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(dynamic_system.matrix);
  
  carl::write_PETSC_matrix(*temp_mat_ptr, input_params.output_base_matrix + "_sys_mat.petscmat");

  //---------------------------------
  // Ajout pour l'export de M, K et C
  libMesh::PetscMatrix<libMesh::Number> * M_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(&dynamic_system.get_matrix("mass"));
  libMesh::PetscMatrix<libMesh::Number> * K_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(&dynamic_system.get_matrix("stiffness"));
  libMesh::PetscMatrix<libMesh::Number> * C_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(&dynamic_system.get_matrix("damping"));

  //dynamic_system.matrix->print_matlab(input_params.output_base_matrix + "A_sys_mat.m");
  //dynamic_system.get_matrix("mass").print_matlab(input_params.output_base_matrix + "M_sys_mat.m");
  //dynamic_system.get_matrix("stiffness").print_matlab(input_params.output_base_matrix + "K_sys_mat.m");
  //dynamic_system.get_matrix("damping").print_matlab(input_params.output_base_matrix + "C_sys_mat.m");
 
  carl::write_PETSC_matrix(*M_ptr, input_params.output_base_matrix + "_M_sys_mat.petscmat");
  carl::write_PETSC_matrix(*K_ptr, input_params.output_base_matrix + "_K_sys_mat.petscmat");
  carl::write_PETSC_matrix(*C_ptr, input_params.output_base_matrix + "_C_sys_mat.petscmat");

  return 0;
}

// This function applies the initial conditions
void apply_initial(libMesh::EquationSystems & es,
                   const std::string & system_name)
{
  // Get a reference to our system, as before
  libMesh::NewmarkSystem & t_system = es.get_system<libMesh::NewmarkSystem> (system_name);

  // Numeric vectors for the pressure, velocity and acceleration
  // values.
  libMesh::NumericVector<libMesh::Number> & disp_vec = t_system.get_vector("displacement");
  libMesh::NumericVector<libMesh::Number> & vel_vec  = t_system.get_vector("velocity");
  libMesh::NumericVector<libMesh::Number> & acc_vec  = t_system.get_vector("acceleration");

  // Assume our fluid to be at rest, which would
  // also be the default conditions in class NewmarkSystem,
  // but let us do it explicetly here.
  disp_vec.zero();
  vel_vec.zero();
  acc_vec.zero();
}

// This function applies the Dirichlet boundary conditions
void fill_mat_dirichlet_bc(libMesh::EquationSystems & es,
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
  libMesh::SparseMatrix<libMesh::Number> & matrix = *t_system.matrix;

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
      const libMesh::Real z_coo = 5.;

      // The penalty parameter.
      const libMesh::Real penalty = 1.e10;

      if (std::abs(curr_node(2)-z_coo) < 1.e-6)
        {

          // The global number of the respective degree of freedom.
          unsigned int dn = curr_node.dof_number(sys_num, 2, 0);
          
          std::cout
          << "MAT | rank " << es.comm().rank()
          << " | node " << curr_node.id()
          << " | dof " << dn
          << " | penalty " << penalty
          << std::endl;
          
          matrix.add(dn, dn, penalty);
          
        }
    } // loop n_cnt

    matrix.close();
    
}