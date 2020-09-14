#include "assemble_functions_elasticity_3D.h"
#include "assemble_functions_mass_3D.h"
#include "elasticity_system.h"
#include "PETSC_matrix_operations.h"

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
/* libMesh::SparseMatrix< libMesh::Number > * get_stiffness_matrix(libMesh::EquationSystems& es,
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
*/

void Update_SubM( libMesh::DenseSubMatrix<libMesh::Number>& SubM,
  unsigned int qp,const std::vector<std::vector<libMesh::Real>> & phi,
  const unsigned int n_components,const unsigned int n_u_dofs,
  const std::vector<libMesh::Real>& JxW,double weight)
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
void assemble_mass_tilde_with_weight(libMesh::EquationSystems& es, 
  const std::string& system_name, weight_parameter_function& weight_mask, 
  WeightFunctionSystemType system_type, libmesh_assemble_input_params& input_params)
{
  libmesh_assert_equal_to (system_name, "Elasticity");

  libMesh::PerfLog perf_log ("Mass/Stiffness Matrix Assembly ",MASTER_bPerfLog_assemble_fem);

  perf_log.push("Preamble");

  libMesh::Real delta_t = input_params.deltat;
  libMesh::Real beta = 0.25;
  PetscScalar betaDeltat2 = beta*delta_t*delta_t;

  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  // - Set up physical properties system ------------------------------------
  libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
  const unsigned int young_var = physical_param_system.variable_number ("E");
  const unsigned int mu_var = physical_param_system.variable_number ("mu");
  const unsigned int rho_var = physical_param_system.variable_number ("rho");

  const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
  std::vector<libMesh::dof_id_type> physical_dof_indices_var;

  // The DoF and values of the physical system
  const libMesh::Elem*  phys_elem   = *(mesh.active_local_elements_begin());
  physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, young_var);
  libMesh::Number localE = physical_param_system.current_solution(physical_dof_indices_var[0]);

  physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, mu_var);
  libMesh::Number localMu = physical_param_system.current_solution(physical_dof_indices_var[0]);

  physical_dof_map.dof_indices(phys_elem, physical_dof_indices_var, rho_var);
  libMesh::Number localRho = physical_param_system.current_solution(physical_dof_indices_var[0]);

  // - Set up elasticity system ---------------------------------------------
  ElasticitySystem& system = es.get_system<ElasticitySystem>("Linear Elasticity");

  libMesh::SparseMatrix< libMesh::Number > & mass = system.add_matrix("mass_tilde");
  
  const unsigned int n_components = 3;
  const unsigned int u_var = system.variable_number("u");
  const unsigned int v_var = system.variable_number("v");
  const unsigned int w_var = system.variable_number("w");

  const DofMap & dof_map = system.get_dof_map();
  libMesh::FEType fe_type = dof_map.variable_type(u_var);

  // Set up pointers to FEBase's of dimension dim and FE type fe_type
  // -> 3D elements
  libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
  libMesh::QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  // -> Faces
  libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));
  libMesh::QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  // Jacobian
  const std::vector<libMesh::Real>& JxW = fe->get_JxW();

  // Shape functions derivatives
  const std::vector<std::vector<libMesh::Real>> & phi = fe->get_phi();
  const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

  // Quadrature points
  const std::vector<libMesh::Point>& qp_points = fe->get_xyz();

  // Weights for the Arlequin method
  double weight_alpha = 1;

  libMesh::DenseMatrix<libMesh::Number> Ke;
  libMesh::DenseMatrix<libMesh::Number> Me;
  libMesh::DenseVector<libMesh::Number> Fe;

  libMesh::DenseSubMatrix<libMesh::Number>
  Kuu(Ke), Kuv(Ke), Kuw(Ke),
  Kvu(Ke), Kvv(Ke), Kvw(Ke),
  Kwu(Ke), Kwv(Ke), Kww(Ke);

  libMesh::DenseSubMatrix<libMesh::Number>
  Muu(Me), Muv(Me), Muw(Me),
  Mvu(Me), Mvv(Me), Mvw(Me),
  Mwu(Me), Mwv(Me), Mww(Me);

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

  perf_log.pop("Preamble");

  // For each element
  for ( ; el != end_el; ++el)
  {
    // Get its pointer
    const libMesh::Elem* elem = *el;

    perf_log.push("Define DoF");
    // The total DoF indices, and those associated to each variable
    
    dof_map.dof_indices (elem, dof_indices);
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_v, v_var);
    dof_map.dof_indices (elem, dof_indices_w, w_var);

    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size();
    const unsigned int n_v_dofs = dof_indices_v.size();
    const unsigned int n_w_dofs = dof_indices_w.size();

    perf_log.pop("Define DoF");

    // Restart the FE to the "geometry" of the element
    // -> Determines quadrature points, shape functions ...
    // !!! User can change the points to be used (can use other mesh's points
    //    instead of the quadrature points)
    fe->reinit (elem);

    perf_log.push("Matrix manipulations");
    Ke.resize (n_dofs, n_dofs);
    Me.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);

    // Set the positions of the sub-matrices
    Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
    Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
    Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);

    Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
    Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
    Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);

    Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
    Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
    Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
    
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

    perf_log.push("Stiffness/Mass Matrix manipulations");

    // For each quadrature point determinate the sub-matrices elements
    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {

      perf_log.pop("Mass","Matrix element calculations");
      weight_alpha = weight_mask.get_alpha(qp_points[qp],system_type);

      // Internal tension
      Update_SubK_isotropic(Kuu, qp, 0, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
      Update_SubK_isotropic(Kuv, qp, 0, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
      Update_SubK_isotropic(Kuw, qp, 0, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);

      Update_SubK_isotropic(Kvu, qp, 1, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
      Update_SubK_isotropic(Kvv, qp, 1, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
      Update_SubK_isotropic(Kvw, qp, 1, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);

      Update_SubK_isotropic(Kwu, qp, 2, 0, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
      Update_SubK_isotropic(Kwv, qp, 2, 1, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);
      Update_SubK_isotropic(Kww, qp, 2, 2, dphi, n_components, n_u_dofs, JxW, localE, localMu, weight_alpha);

      Update_SubM(Muu, qp, phi, dim, n_u_dofs, JxW, localRho);
      Update_SubM(Muv, qp, phi, dim, n_u_dofs, JxW, localRho);
      Update_SubM(Muw, qp, phi, dim, n_u_dofs, JxW, localRho);

      Update_SubM(Mvu, qp, phi, dim, n_u_dofs, JxW, localRho);
      Update_SubM(Mvv, qp, phi, dim, n_u_dofs, JxW, localRho);
      Update_SubM(Mvw, qp, phi, dim, n_u_dofs, JxW, localRho);

      Update_SubM(Mwu, qp, phi, dim, n_u_dofs, JxW, localRho);
      Update_SubM(Mwv, qp, phi, dim, n_u_dofs, JxW, localRho);
      Update_SubM(Mww, qp, phi, dim, n_u_dofs, JxW, localRho);

      perf_log.pop("Stiffness/Mass matrix","Matrix element calculations");
    }

    // Apply constraints
    perf_log.push("Constraints","Matrix element calculations");
    dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    dof_map.heterogenously_constrain_element_matrix_and_vector(Me, Fe, dof_indices);
    perf_log.pop("Constraints","Matrix element calculations");

    perf_log.push("Adding elements");
    
    Me.add(betaDeltat2,Ke);

    mass.add_matrix(Me,dof_indices);
    system.matrix->add_matrix(Me, dof_indices);

    system.rhs->add_vector(Fe, dof_indices);
    perf_log.pop("Adding elements");
  }

  system.matrix->close();
  system.rhs->close();


  libMesh::PetscMatrix<libMesh::Number> * temp_mass_tilde = libMesh::cast_ptr<libMesh::PetscMatrix<libMesh::Number>* >(system.request_matrix("mass_tilde"));
  carl::write_PETSC_matrix(*temp_mass_tilde, input_params.output_base + "_mass_tilde_with_weights.petscmat");

#ifdef PRINT_MATLAB_DEBUG
  temp_mass_tilde->print_matlab(input_params.output_base + "_mass_tilde_with_weights.m");
#endif


}

/*
void get_mass_tilde(libMesh::Parallel::Communicator& WorldComm, 
  libmesh_assemble_input_params& input_params, libMesh::EquationSystems equation_systems)
{
  

// Export matrix and vector

  
  //libMesh::PetscVector<libMesh::Number> * temp_vec_ptr = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> * >(elasticity_system.rhs);

 

  //libMesh::SparseMatrix< libMesh::Number > * mass = get_mass_matrix(equation_systems,"Linear Elasticity");
  //libMesh::SparseMatrix< libMesh::Number > * stiffness = get_stiffness_matrix(equation_systems,"Linear Elasticity");
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


std::string mass_tilde_path = "./GC_solver/mass_tilde/"; 
std::vector<std::string> str;
split(input_params.output_base, str,'/');

std::string name_matrix_tilde_matlab = mass_tilde_path+"mass_tilde_"+str[str.size()-1]+".m";
std::string name_matrix_tilde_petsc  = mass_tilde_path+"mass_tilde_"+str[str.size()-1]+".petscmat";



carl::write_PETSC_matrix(sys_mat_PETSC_M_tilde,name_matrix_tilde_petsc);

  return;

}

*/