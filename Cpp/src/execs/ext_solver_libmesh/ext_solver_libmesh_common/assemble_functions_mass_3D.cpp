#include "assemble_functions_elasticity_3D.h"
#include "assemble_functions_mass_3D.h"
#include "PETSC_matrix_operations.h"

using namespace libMesh;

void Update_SubM( libMesh::DenseSubMatrix<libMesh::Number>& SubM,
  unsigned int qp,const std::vector<std::vector<libMesh::Real>> & phi,
  const unsigned int n_components,const unsigned int n_u_dofs,
  const std::vector<libMesh::Real>& JxW, double weight)
{
  for (unsigned int iii=0; iii<n_u_dofs; iii++)
  {
    for (unsigned int jjj=0; jjj<n_u_dofs; jjj++)
    {
      SubM(iii,jjj) += weight * JxW[qp]*phi[iii][qp]*phi[jjj][qp];
    }
  }
}

// Dynamic elasticity with no external forces
void assemble_dynamic_elasticity_with_weight(libMesh::EquationSystems& es,
  const std::string& system_name, weight_parameter_function& weight_mask,
  WeightFunctionSystemType system_type, double deltat, double beta)
{

  libmesh_assert_equal_to(system_name, "Dynamic Elasticity");

  libMesh::PerfLog perf_log ("Mass/Stiffness Matrix Assembly ",MASTER_bPerfLog_assemble_fem);

  perf_log.push("Preamble");

  //libMesh::Real deltat = input_params.deltat;
  //libMesh::Real beta = input_params.beta;
  PetscScalar betaDeltat2 = beta*deltat*deltat;

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
  /* [NEWMARK] */ libMesh::NewmarkSystem& system = es.get_system<libMesh::NewmarkSystem>("Dynamic Elasticity");
  ///* [IMPLICIT]*/ libMesh::LinearImplicitSystem& system = es.get_system<libMesh::LinearImplicitSystem>("Dynamic Elasticity");

  libMesh::SparseMatrix < libMesh::Number > & stiffness = system.get_matrix("stiffness");
  libMesh::SparseMatrix < libMesh::Number > & damping   = system.get_matrix("damping");
  libMesh::SparseMatrix < libMesh::Number > & mass      = system.get_matrix("mass");
  libMesh::NumericVector< libMesh::Number > & force     = system.get_vector("force");
  libMesh::DenseMatrix< libMesh::Number > zero_matrix;
  libMesh::SparseMatrix < libMesh::Number > & mass_tilde = system.add_matrix("mass_tilde");

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
  libMesh::DenseMatrix<libMesh::Number> Ce;
  libMesh::DenseMatrix<libMesh::Number> Me;
  libMesh::DenseVector<libMesh::Number> Fe;

  libMesh::DenseSubMatrix<libMesh::Number>
  Kuu(Ke), Kuv(Ke), Kuw(Ke),
  Kvu(Ke), Kvv(Ke), Kvw(Ke),
  Kwu(Ke), Kwv(Ke), Kww(Ke);

  libMesh::DenseSubMatrix<libMesh::Number>
  Cuu(Ce), Cuv(Ce), Cuw(Ce),
  Cvu(Ce), Cvv(Ce), Cvw(Ce),
  Cwu(Ce), Cwv(Ce), Cww(Ce);

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
    Ce.resize (n_dofs, n_dofs);
    Me.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);
    //zero_matrix.resize(n_dofs, n_dofs);
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
    Cuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
    Cuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
    Cuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);

    Cvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
    Cvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
    Cvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);

    Cwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
    Cwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
    Cww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
    
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

      perf_log.push("Element Mass/Stiffness","Matrix element calculations");
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

      Update_SubM(Muu, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Muv, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Muw, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);

      Update_SubM(Mvu, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Mvv, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Mvw, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);

      Update_SubM(Mwu, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Mwv, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Mww, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);

      perf_log.pop("Element Mass/Stiffness","Matrix element calculations");
    }
    perf_log.pop("Stiffness/Mass Matrix manipulations");

    // Apply constraints
    perf_log.push("Constraints","Matrix element calculations");
    std::vector<unsigned int> dof_indicesC = dof_indices;
    std::vector<unsigned int> dof_indicesM = dof_indices;
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    dof_map.constrain_element_matrix(Ce, dof_indicesC);
    dof_map.constrain_element_matrix(Me, dof_indicesM);
    perf_log.pop("Constraints","Matrix element calculations");

    perf_log.push("Adding elements");
    stiffness.add_matrix(Ke,dof_indices);
    force.add_vector(Fe,dof_indices);
    Me.add(betaDeltat2,Ke);
    system.matrix->add_matrix(Me, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
    perf_log.pop("Adding elements");
  } // end element loop
  system.matrix->close();
  system.rhs->close();
}


// Dynamic elasticity with external traction
void assemble_dynamic_elasticity_with_weight_and_traction(libMesh::EquationSystems& es,
                const std::string& system_name, 
                weight_parameter_function& weight_mask,
                WeightFunctionSystemType system_type,
                int traction_boundary_id,
                std::vector<double> traction_density,
                double deltat, double beta)
{

  libmesh_assert_equal_to(system_name, "Dynamic Elasticity");

  libMesh::PerfLog perf_log ("Mass/Stiffness Matrix Assembly ",MASTER_bPerfLog_assemble_fem);

  perf_log.push("Preamble");

  //libMesh::Real deltat = input_params.deltat;
  //libMesh::Real beta = input_params.beta;
  //libMesh::Real gamma = input_params.gamma;
  PetscScalar betaDeltat2 = beta*deltat*deltat;

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
  /* [NEWMARK] */ libMesh::NewmarkSystem& system = es.get_system<libMesh::NewmarkSystem>("Dynamic Elasticity");
  ///* [IMPLICIT] */ libMesh::LinearImplicitSystem& system = es.get_system<libMesh::LinearImplicitSystem>("Elasticity");

  libMesh::SparseMatrix < libMesh::Number > & stiffness = system.get_matrix("stiffness");
  libMesh::SparseMatrix < libMesh::Number > & damping   = system.get_matrix("damping");
  libMesh::SparseMatrix < libMesh::Number > & mass      = system.get_matrix("mass");
  libMesh::NumericVector< libMesh::Number > & force     = system.get_vector("force");
  libMesh::DenseMatrix< libMesh::Number > zero_matrix;
  libMesh::SparseMatrix < libMesh::Number > & mass_tilde = system.add_matrix("mass_tilde");
  
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

  // -> Face traction
  libMesh::DenseVector<libMesh::Number> g_vec(3);
  g_vec(0) = traction_density[0];
  g_vec(1) = traction_density[1];
  g_vec(2) = traction_density[2];
  libMesh::Real amp_n = sqrt(g_vec(0)*g_vec(0)+g_vec(1)*g_vec(1)+g_vec(2)*g_vec(2));

  // Jacobian
  const std::vector<libMesh::Real>& JxW = fe->get_JxW();

  // Shape functions derivatives
  const std::vector<std::vector<libMesh::Real>> & phi = fe->get_phi();
  const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

  // Quadrature points
  const std::vector<libMesh::Point>& qp_points = fe->get_xyz();

  // Face Jacobian
  const std::vector<libMesh::Real> & JxW_face = fe_face->get_JxW();

  // Face shape functions
  const std::vector<std::vector<libMesh::Real> > & phi_face = fe_face->get_phi();

  // Face quadrature points
  const std::vector<libMesh::Point>& qp_face_points = fe_face->get_xyz();

  // Weights for the Arlequin method
  double weight_alpha = 1;

  libMesh::DenseMatrix<libMesh::Number> Ke;
  libMesh::DenseMatrix<libMesh::Number> Ce;
  libMesh::DenseMatrix<libMesh::Number> Me;
  libMesh::DenseVector<libMesh::Number> Fe;

  libMesh::DenseSubMatrix<libMesh::Number>
  Kuu(Ke), Kuv(Ke), Kuw(Ke),
  Kvu(Ke), Kvv(Ke), Kvw(Ke),
  Kwu(Ke), Kwv(Ke), Kww(Ke);

  libMesh::DenseSubMatrix<libMesh::Number>
  Cuu(Ce), Cuv(Ce), Cuw(Ce),
  Cvu(Ce), Cvv(Ce), Cvw(Ce),
  Cwu(Ce), Cwv(Ce), Cww(Ce);

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
    Ce.resize (n_dofs, n_dofs);
    Me.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);
    //zero_matrix.resize(n_dofs, n_dofs);

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
    Cuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
    Cuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
    Cuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);

    Cvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
    Cvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
    Cvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);

    Cwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
    Cwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
    Cww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);

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

      perf_log.push("Element Mass/Stiffness","Matrix element calculations");
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

      Update_SubM(Muu, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Muv, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Muw, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);

      Update_SubM(Mvu, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Mvv, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Mvw, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);

      Update_SubM(Mwu, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Mwv, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);
      Update_SubM(Mww, qp, phi, dim, n_u_dofs, JxW, localRho*weight_alpha);

      perf_log.pop("Element Mass/Stiffness","Matrix element calculations");
    }

    // Assemble the traction
    for (unsigned int side=0; side<elem->n_sides(); side++)
    {
      if (elem->neighbor_ptr(side) == libmesh_nullptr)
      {
        fe_face->reinit(elem, side);

        // Apply a traction
        for (unsigned int qp=0; qp<qface.n_points(); qp++)
        {
          weight_alpha = weight_mask.get_alpha(qp_face_points[qp],system_type);

          if (mesh.get_boundary_info().has_boundary_id(elem, side, traction_boundary_id))
          { 
            for (unsigned int dof_iii=0; dof_iii<n_u_dofs; dof_iii++)
            {
              Fu(dof_iii) += amp_n * weight_alpha * JxW_face[qp] * (g_vec(0) * phi_face[dof_iii][qp]);
              Fv(dof_iii) += amp_n * weight_alpha * JxW_face[qp] * (g_vec(1) * phi_face[dof_iii][qp]);
              Fw(dof_iii) += amp_n * weight_alpha * JxW_face[qp] * (g_vec(2) * phi_face[dof_iii][qp]);
            }
          }
        }
      }
    }
    perf_log.pop("Stiffness/Mass Matrix manipulations");

    // Apply constraints
    perf_log.push("Constraints","Matrix element calculations");
    std::vector<unsigned int> dof_indicesC = dof_indices;
    std::vector<unsigned int> dof_indicesM = dof_indices;
    dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    dof_map.constrain_element_matrix(Ce, dof_indicesC);
    dof_map.constrain_element_matrix(Me, dof_indicesM);
    perf_log.pop("Constraints","Matrix element calculations");

    perf_log.push("Adding elements");
    stiffness.add_matrix(Ke,dof_indices);
    force.add_vector(Fe,dof_indices);
    Me.add(betaDeltat2,Ke);
    system.matrix->add_matrix(Me, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
    perf_log.pop("Adding elements");
  }
  system.matrix->close();
  system.rhs->close();
}

// Function Prototype.  This function will be used to apply the
// initial conditions.
void apply_initial(EquationSystems & es,
                   const std::string & system_name,
                   const bool & zeroed,
                   const std::string& dis_vec_name,
                   const std::string& vel_vec_name,
                   const std::string& acc_vec_name,
                   Parallel::Communicator& WorldComm)

{
  libmesh_assert_equal_to(system_name, "Dynamic Elasticity");

  // Get a reference to our system, as before
  NewmarkSystem & t_system = es.get_system<NewmarkSystem> (system_name);
  
  NumericVector<Number> & dis_vec  = t_system.get_vector("displacement");
  NumericVector<Number> & vel_vec  = t_system.get_vector("velocity");
  NumericVector<Number> & acc_vec  = t_system.get_vector("acceleration");
  dis_vec.zero();
  vel_vec.zero();
  acc_vec.zero();
  // Numeric vectors for the pressure, velocity and acceleration
  // values.
  if (!zeroed) {
    Vec dis_vec_PETSC, vel_vec_PETSC, acc_vec_PETSC;
    VecCreate(WorldComm.get(),&dis_vec_PETSC);
    VecCreate(WorldComm.get(),&vel_vec_PETSC);
    VecCreate(WorldComm.get(),&acc_vec_PETSC);
    carl::read_PETSC_vector(dis_vec_PETSC,dis_vec_name,WorldComm.get());
    carl::read_PETSC_vector(vel_vec_PETSC,vel_vec_name,WorldComm.get());
    carl::read_PETSC_vector(acc_vec_PETSC,acc_vec_name,WorldComm.get());
    PetscVector<Number> dis_vec(dis_vec_PETSC,WorldComm);
    PetscVector<Number> vel_vec(vel_vec_PETSC,WorldComm);
    PetscVector<Number> acc_vec(acc_vec_PETSC,WorldComm);
  }
}

void update_dynamic_rhs(libMesh::EquationSystems& es,
                const std::string& system_name, 
                weight_parameter_function& weight_mask,
                WeightFunctionSystemType system_type, 
                double amp_n,
                const std::string& stf_mat_file,
                const std::string& fex_vec_name,
                const std::string& dis_vec_name,
                const std::string& vel_vec_name,
                const std::string& acc_vec_name,
                Parallel::Communicator& WorldComm,
                double deltat,
                double beta)
{
  libmesh_assert_equal_to(system_name, "Dynamic Elasticity");

  PerfLog perf_log ("Mass/Stiffness Matrix Assembly ",MASTER_bPerfLog_assemble_fem);

  perf_log.push("Preamble");
  
  //libMesh::Real deltat = input_params.deltat;
  //libMesh::Real beta = input_params.beta;
  //libMesh::Real gamma = input_params.gamma;
  PetscScalar betaDeltat2 = beta*deltat*deltat;

  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();
  
  // - Set up elasticity system ---------------------------------------------
  NewmarkSystem& system = es.get_system<NewmarkSystem>("Dynamic Elasticity");

  Vec dis_vec_PETSC;
  Vec vel_vec_PETSC;
  Vec acc_vec_PETSC;
  Vec fex_vec_PETSC;
  Mat stf_mat_PETSC;
  VecCreate(WorldComm.get(),&dis_vec_PETSC);
  VecCreate(WorldComm.get(),&vel_vec_PETSC);
  VecCreate(WorldComm.get(),&acc_vec_PETSC);
  VecCreate(WorldComm.get(),&fex_vec_PETSC);
  MatCreate(WorldComm.get(),&stf_mat_PETSC);
  carl::read_PETSC_vector(dis_vec_PETSC,dis_vec_name,WorldComm.get());
  carl::read_PETSC_vector(vel_vec_PETSC,vel_vec_name,WorldComm.get());
  carl::read_PETSC_vector(acc_vec_PETSC,acc_vec_name,WorldComm.get());
  carl::read_PETSC_vector(fex_vec_PETSC,fex_vec_name,WorldComm.get());
  carl::read_PETSC_matrix(stf_mat_PETSC,stf_mat_file,WorldComm.get());
  PetscVector<Number> dis_vec(dis_vec_PETSC,WorldComm);
  PetscVector<Number> vel_vec(vel_vec_PETSC,WorldComm);
  PetscVector<Number> acc_vec(acc_vec_PETSC,WorldComm);
  PetscVector<Number> fex_vec(fex_vec_PETSC,WorldComm);
  PetscMatrix<Number> stf_mat(stf_mat_PETSC,WorldComm);

  // Zeroes the rhs
  system.rhs->zero();

  // update external force amplitude
  fex_vec.scale(amp_n);

  // compute auxiliary vectors rhs_m
  libMesh::Real a0 = deltat*deltat*(0.5-beta);

  // Compute correction with displacement predictor
  dis_vec.scale(-1.);
  dis_vec.add(-deltat, vel_vec);
  dis_vec.add(-a0, acc_vec);

  system.rhs->add(fex_vec);
  system.rhs->add_vector(dis_vec,stf_mat);

  // --- Cleanup!
  VecDestroy(&dis_vec_PETSC);
  VecDestroy(&vel_vec_PETSC);
  VecDestroy(&acc_vec_PETSC);
  VecDestroy(&fex_vec_PETSC);
  MatDestroy(&stf_mat_PETSC);

  return;
}
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */