#include "libmesh_assemble_lin_homogeneous.h"

/** \brief Program used to assemble the rigidity matrix and the vectors of a linear, homogeneous elasticity model with a clamped \f$x_{\mbox{Min}}\f$ face.
 * 
 *  Usage: `./libmesh_assemble_lin_homogeneous__min_x_clamped -i [input file]`
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

int main(int argc, char** argv) {

  // [USER] Fixed boudary
  boundary_id_cube fixed_bound_id = boundary_id_cube::MIN_X;
  
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

  libmesh_assemble_input_params input_params;
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

  // Add linear elasticity and physical parameters systems
  libMesh::LinearImplicitSystem& elasticity_system
                    = add_elasticity(equation_systems);

  // Set clamped border
  set_clamped_border(elasticity_system, fixed_bound_id);

  // Initialize the equation systems
  equation_systems.init();

  // Homogeneous properties for the macro system
  set_homogeneous_physical_properties(equation_systems, input_params.physical_params_file);

  // Set the weight function object
  weight_parameter_function  system_weight(mesh_weight);
  system_weight.set_parameters(input_params.weight_domain_idx_file);

  perf_log.pop("System setup:");

  // Assemble!
  assemble_elasticity_with_weight(equation_systems,"Elasticity",system_weight,
              input_params.system_type);

// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  elasticity_system.matrix->print_matlab(input_params.output_base + "_sys_mat.m");
  elasticity_system.rhs->print_matlab(input_params.output_base + "_sys_rhs_vec.m");
#endif

  // Export matrix and vector
  libMesh::PetscMatrix<libMesh::Number> * temp_mat_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(elasticity_system.matrix);
  libMesh::PetscVector<libMesh::Number> * temp_vec_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> * >(elasticity_system.rhs);

  carl::write_PETSC_matrix(*temp_mat_ptr, input_params.output_base + "_sys_mat.petscmat");
  carl::write_PETSC_vector(*temp_vec_ptr, input_params.output_base + "_sys_rhs_vec.petscvec");

  // If needed, print rigid body vectors
  if(input_params.bCalculateRBVectors)
  {
    MatNullSpace nullsp_sys;
    build_rigid_body_vectors(elasticity_system,nullsp_sys);
    write_rigid_body_vectors(nullsp_sys,input_params.output_base,WorldComm.rank());
    MatNullSpaceDestroy(&nullsp_sys);
  }

  return 0;
}
