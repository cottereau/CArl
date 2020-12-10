#include "libmesh_assemble_lin_homogeneous.h"

/** \brief Program used to assemble the rigidity matrix and the vectors of a linear, homogeneous elasticity model with a traction applied to \f$x_{\mbox{Max}}\f$ face.
 * 
 *  Usage: `./libmesh_assemble_lin_homogeneous__max_x_traction -i [input file]`
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

using namespace std;
using namespace libMesh;

int main(int argc, char** argv) {

  // [USER] Traction force density
  std::vector<double> traction_density(3,0);
  traction_density[0] = 100.;
  
  // --- Initialize libMesh
  LibMeshInit init(argc, argv);

  // Do performance log?
  const bool MASTER_bPerfLog_carl_libmesh = true;
  PerfLog perf_log("Main program", MASTER_bPerfLog_carl_libmesh);

  // libMesh's C++ / MPI communicator wrapper
  Parallel::Communicator& WorldComm = init.comm();

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

  // Declaration of a struct variable to catch all inputs from the input file
  libmesh_assemble_input_params input_params;

  // Getting parsing argument from all file
  get_input_params(field_parser, input_params);

  // Check libMesh installation dimension
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // --- Declare the three meshes to be intersected

  // - Parallelized meshes: A, B, mediator and weight
  perf_log.push("Meshes - Parallel","Read files:");
  Mesh system_mesh(WorldComm, dim);
  system_mesh.read(input_params.mesh_file);
  system_mesh.prepare_for_use();

  Mesh mesh_weight(WorldComm, dim);
  mesh_weight.allow_renumbering(false);
  mesh_weight.read(input_params.mesh_weight_file);
  mesh_weight.prepare_for_use();

  perf_log.pop("Meshes - Parallel","Read files:");

  // --- Generate the equation systems
  perf_log.push("System setup:");

  // Set the equation systems object
  EquationSystems equation_systems(system_mesh);

  // Add linear elasticity and physical parameters systems
  // Create NewmarkSystem "Elasticity"
  NewmarkSystem& elasticity_system = add_dynamic_elasticity(equation_systems,
    "Dynamic Elasticity");
  //LinearImplicitSystem elasticity_system = add_dynamic_elasticity(equation_systems,
    

  // Set Newmark's parameters
  elasticity_system.set_newmark_parameters(input_params.deltat,
    input_params.beta,
    input_params.gamma);

  // Start time integration from t=0
  elasticity_system.time = 0.;

  // Set the weight function object
  weight_parameter_function  system_weight(mesh_weight);
  system_weight.set_parameters(input_params.weight_domain_idx_file);

  perf_log.pop("System setup:");

  // Initialize the equation systems
  equation_systems.init();

  // Set material parameter
  set_homogeneous_physical_properties(equation_systems, input_params.physical_params_file);
  
  // Set material parameter
  set_homogeneous_physical_properties(equation_systems, input_params.physical_params_file);
  
  // Initial conditions 
  apply_initial(equation_systems,"Dynamic Elasticity",
                true,
                input_params.dis_vec_name,
                input_params.vel_vec_name,
                input_params.acc_vec_name,
                WorldComm);
  
  // Assemble!
  assemble_dynamic_elasticity_with_weight_and_traction(equation_systems,"Dynamic Elasticity",
    system_weight, input_params.system_type, boundary_id_cube::MAX_X, traction_density, 
    input_params.deltat, input_params.beta);

  SparseMatrix < Number > * stiffness = elasticity_system.request_matrix("stiffness");
// Print MatLab debugging output? Variable defined at "carl_headers.h"
#ifdef PRINT_MATLAB_DEBUG
  elasticity_system.matrix->print_matlab(input_params.output_base + "_sys_mat.m");
  elasticity_system.rhs->print_matlab(input_params.output_base + "_sys_rhs_vec.m");
  stiffness->print_matlab(input_params.output_base + "_K_mat.m");
#endif

  // Export matrix and vector
  PetscMatrix<Number> * temp_mat_ptr = cast_ptr<PetscMatrix<Number> * >(elasticity_system.matrix);
  PetscVector<Number> * temp_vec_ptr = cast_ptr<PetscVector<Number> * >(elasticity_system.rhs);
  PetscMatrix<Number> * temp_stf_ptr = cast_ptr<PetscMatrix<Number> * >(stiffness);

  carl::write_PETSC_matrix(*temp_mat_ptr, input_params.output_base + "_sys_mat.petscmat");
  carl::write_PETSC_vector(*temp_vec_ptr, input_params.output_base + "_sys_rhs_vec.petscvec");
  carl::write_PETSC_matrix(*temp_stf_ptr, input_params.output_base + "_K_mat.petscmat");

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
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */