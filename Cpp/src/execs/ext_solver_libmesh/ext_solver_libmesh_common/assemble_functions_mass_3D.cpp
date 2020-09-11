#include "assemble_functions_elasticity_3D.h"
#include "assemble_functions_mass_3D.h"
#include "libmesh_assemble_system_input_parser.h"

using namespace libMesh;
/*
void get_input_params(GetPot& field_parser,
  libmesh_assemble_interpolation_params& input_params) 
{

  if(field_parser.search(1, "PathTildeMatrixA"))
    input_params.path_tilde_matrices_A = field_parser.next(input_params.path_tilde_matrices_A);
  else
    printf("[ERROR] Matrices system A found");

  if(field_parser.search(1, "MacroCouplingMatrix"))
    input_params.path_macro_coupling_matrix = field_parser.next(input_params.path_macro_coupling_matrix);
  else
    printf("[ERROR] Path to Macro Coupling Matrix not found\n");

  if(field_parser.search(1,"AlphaA"))
    input_params.alpha_A = field_parser.next(input_params.alpha_A);
  else
    printf("[ERROR] Alpha value not passed to the program\n");

  if(field_parser.search(1, "PathTildeMatrixB"))
    input_params.path_tilde_matrices_B = field_parser.next(input_params.path_tilde_matrices_B);
  else
    printf("[ERROR] Matrices system B found");

  if(field_parser.search(1,"MicroCouplingMatrix"))
    input_params.path_micro_coupling_matrix = field_parser.next(input_params.path_micro_coupling_matrix);
  else
    printf("[ERROR] Path to Micro Coupling Matrix not found");

  if(field_parser.search(1,"AlphaB"))
    input_params.alpha_B = field_parser.next(input_params.alpha_B);
  else
    printf("[ERROR] Alpha value not passed to the program\n");

  if(field_parser.search(1, "OutputFolder"))
    input_params.output_folder = field_parser.next(input_params.output_folder);
  else
    printf("[ERROR] Matrices system A found");
}
*/
//Matrice of rigidity
libMesh::SparseMatrix< libMesh::Number > * get_stiffness_matrix(libMesh::EquationSystems& es,
        const std::string& system_name,
        weight_parameter_function& weight_mask,
        WeightFunctionSystemType system_type)
{
  assemble_elasticity_with_weight(es,system_name,weight_mask,system_type);
  libMesh::LinearImplicitSystem& system = es.get_system<libMesh::LinearImplicitSystem>("Elasticity");
  libMesh::SparseMatrix< libMesh::Number > * stiffness = system.matrix;
  return stiffness;
}

libMesh::SparseMatrix< libMesh::Number > * get_mass_matrix(libMesh::EquationSystems& es,
        const std::string& system_name)
{
  libMesh::SparseMatrix< libMesh::Number > * mass = assemble_mass_matrix(es,system_name);
  return mass;
}

void Update_SubM( libMesh::DenseSubMatrix<libMesh::Number>& SubM,
          unsigned int qp,
          const std::vector<std::vector<libMesh::Real>> & phi,
          const unsigned int n_components,
          const unsigned int n_u_dofs,
          const std::vector<libMesh::Real>& JxW,
          double weight)
{
  for (unsigned int iii=0; iii<n_u_dofs; iii++)
  {
    for (unsigned int jjj=0; jjj<n_u_dofs; jjj++)
    {
      SubM(iii,jjj) += weight * JxW[qp]*phi[iii][qp]*phi[jjj][qp];
    }
  }
}

//
libMesh::SparseMatrix< libMesh::Number > * assemble_mass_matrix(libMesh::EquationSystems& es, const std::string& system_name)
{
  //printf("in assemble_mass_matrix functions\n");
  libmesh_assert_equal_to (system_name, "Elasticity");
  libMesh::LinearImplicitSystem& system = es.get_system<libMesh::LinearImplicitSystem>("Elasticity");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  libMesh::SparseMatrix< libMesh::Number > & mass = system.get_matrix("mass");
  
  
  // for the first (and only) variable in the system.
  const unsigned int u_var = system.variable_number("u");
  const unsigned int v_var = system.variable_number("v");
  const unsigned int w_var = system.variable_number("w");

  //printf("before dof_map\n");
  const DofMap & dof_map = system.get_dof_map();
  //printf(">> [dof_map]:%u\n", dof_map.n_dofs());
  //printf("before fe_type\n");
  libMesh::FEType fe_type = dof_map.variable_type(u_var);

  //printf("before fe\n");
  libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
  //printf("before qrule\n");
  libMesh::QGauss qrule (dim, fe_type.default_quadrature_order());
  //printf("before attach_quadrature_rule\n");
  fe->attach_quadrature_rule (&qrule);

  // -> Faces
  //printf("before fe_face\n");
  libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));
  //printf("before qface\n");
  libMesh::QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  double weight = 1000.;

  // Jacobian
  //printf("before get_JxW\n");
  const std::vector<libMesh::Real>& JxW = fe->get_JxW();

  // Shape functions derivatives
  //printf("before get_phi\n");
  const std::vector<std::vector<libMesh::Real>> & phi = fe->get_phi();
  //printf("after get_phi\n");

  // Quadrature points
  const std::vector<libMesh::Point>& qp_points = fe->get_xyz();
  //printf("after qp_points\n");
  // Weights for the Arlequin method
  double weight_alpha = 1;

  //printf("after Me\n");
  libMesh::DenseMatrix<libMesh::Number> Me;
  //printf("after Fe\n");
  libMesh::DenseVector<libMesh::Number> Fe;

  //printf("after Me\n");
  libMesh::DenseSubMatrix<libMesh::Number>
  Muu(Me), Muv(Me), Muw(Me),
  Mvu(Me), Mvv(Me), Mvw(Me),
  Mwu(Me), Mwv(Me), Mww(Me);

  //printf("after `Fe\n");
  libMesh::DenseSubVector<libMesh::Number>
  Fu(Fe),
  Fv(Fe),
  Fw(Fe);

  std::vector<libMesh::dof_id_type> dof_indices;
  std::vector<libMesh::dof_id_type> dof_indices_u;
  std::vector<libMesh::dof_id_type> dof_indices_v;
  std::vector<libMesh::dof_id_type> dof_indices_w;  

  libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  // For each element
  for ( ; el != end_el; ++el)
  {
    // Get its pointer
    const libMesh::Elem* elem = *el;

    // perf_log.push("Define DoF");
    // The total DoF indices, and those associated to each variable
    
    dof_map.dof_indices (elem, dof_indices);
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_v, v_var);
    dof_map.dof_indices (elem, dof_indices_w, w_var);
    //printf("after dof_map.dof_indices\n");

    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size();
    const unsigned int n_v_dofs = dof_indices_v.size();
    const unsigned int n_w_dofs = dof_indices_w.size();
    //printf("after dof_map.dof_indices\n");

    // perf_log.pop("Define DoF");

    // Restart the FE to the "geometry" of the element
    // -> Determines quadrature points, shape functions ...
    // !!! User can change the points to be used (can use other mesh's points
    //    instead of the quadrature points)
    fe->reinit (elem);

    // perf_log.push("Matrix manipulations");
    Me.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);

    // Set the positions of the sub-matrices
    Muu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
    Muv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
    Muw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);

    Mvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
    Mvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
    Mvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);

    Mwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
    Mwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
    Mww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);

    Fu.reposition (u_var*n_u_dofs, n_u_dofs);
    Fv.reposition (v_var*n_u_dofs, n_v_dofs);
    Fw.reposition (w_var*n_u_dofs, n_w_dofs);

    // perf_log.push("Matrix manipulations");

    // For each quadrature point determinate the sub-matrices elements
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      // perf_log.push("Rigidity","Matrix element calculations");
      
      // Internal tension
      Update_SubM(Muu, qp, phi, dim, n_u_dofs, JxW, weight);
      Update_SubM(Muv, qp, phi, dim, n_u_dofs, JxW, weight);
      Update_SubM(Muw, qp, phi, dim, n_u_dofs, JxW, weight);

      Update_SubM(Mvu, qp, phi, dim, n_u_dofs, JxW, weight);
      Update_SubM(Mvv, qp, phi, dim, n_u_dofs, JxW, weight);
      Update_SubM(Mvw, qp, phi, dim, n_u_dofs, JxW, weight);

      Update_SubM(Mwu, qp, phi, dim, n_u_dofs, JxW, weight);
      Update_SubM(Mwv, qp, phi, dim, n_u_dofs, JxW, weight);
      Update_SubM(Mww, qp, phi, dim, n_u_dofs, JxW, weight);

      // perf_log.pop("Rigidity","Matrix element calculations");
    }

    // Apply constraints
    // perf_log.push("Constraints","Matrix element calculations");
    dof_map.heterogenously_constrain_element_matrix_and_vector(Me, Fe, dof_indices);
    // perf_log.pop("Constraints","Matrix element calculations");

    // perf_log.push("Adding elements");
    mass.add_matrix(Me,dof_indices);
    //printf("after mass.add_matrix\n");
    system.matrix->add_matrix(Me, dof_indices);
    //printf("after system.matrix->add_matrix\n");
    system.rhs->add_vector(Fe, dof_indices);
    // perf_log.pop("Adding elements");
  }

  system.matrix->close();
  system.rhs->close();

  return &mass;

}


void get_mass_tilde(libMesh::Parallel::Communicator& WorldComm, 
  libmesh_assemble_input_params& input_params, 
  libMesh::EquationSystems equation_systems)
{
  libMesh::Real delta_t = input_params.deltat;
  libMesh::Real beta = 0.25;
  PetscScalar a = beta*delta_t*delta_t;

// Export matrix and vector

  
  //libMesh::PetscVector<libMesh::Number> * temp_vec_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> * >(elasticity_system.rhs);

  /*
  */

  libMesh::SparseMatrix< libMesh::Number > * mass = get_mass_matrix(equation_systems,"Linear Elasticity");
  //libMesh::SparseMatrix< libMesh::Number > * stiffness = get_stiffness_matrix(equation_systems,"Linear Elasticity");
  libMesh::PetscMatrix<libMesh::Number> * temp_mat_mass_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(mass);
  //libMesh::PetscMatrix<libMesh::Number> * temp_mat_stiffness_ptr = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number> * >(stiffness);
  //carl::write_PETSC_matrix(*temp_mat_stiffness_ptr, input_params.output_base + "_stiffness_test_sys_mat.petscmat");
  //carl::write_PETSC_matrix(*temp_mat_mass_ptr, input_params.output_base + "_mass_test_sys_mat.petscmat");
  //carl::write_PETSC_vector(*temp_vec_ptr, input_params.output_base + "_sys_rhs_vec.petscvec");
  
  Mat sys_mat_PETSC_M;
  Mat sys_mat_PETSC_K;
  MatCreate(WorldComm.get(),&sys_mat_PETSC_M);
  MatCreate(WorldComm.get(),&sys_mat_PETSC_K);

  //carl::read_PETSC_matrix(sys_mat_PETSC_M, input_params.output_base+"_mass_test_sys_mat.petscmat", WorldComm.get());
  //carl::read_PETSC_matrix(sys_mat_PETSC_K, input_params.output_base+"_stiffness_test_sys_mat.petscmat", WorldComm.get());

  //calculate de M~  = M + a*K
  MatAXPY(sys_mat_PETSC_M, a, sys_mat_PETSC_K, DIFFERENT_NONZERO_PATTERN);
  libMesh::PetscMatrix<libMesh::Number> sys_mat_PETSC_M_tilde(sys_mat_PETSC_M,WorldComm);

/*

std::string mass_tilde_path = "./GC_solver/mass_tilde/"; 
std::vector<std::string> str;
split(input_params.output_base, str,'/');

std::string name_matrix_tilde_matlab = mass_tilde_path+"mass_tilde_"+str[str.size()-1]+".m";
std::string name_matrix_tilde_petsc  = mass_tilde_path+"mass_tilde_"+str[str.size()-1]+".petscmat";

#ifdef PRINT_MATLAB_DEBUG
  sys_mat_PETSC_M_tilde.print_matlab(name_matrix_tilde_matlab);
#endif

carl::write_PETSC_matrix(sys_mat_PETSC_M_tilde,name_matrix_tilde_petsc);
*/

}
