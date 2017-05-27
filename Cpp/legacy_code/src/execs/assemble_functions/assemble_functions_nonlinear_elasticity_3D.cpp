/*
 * assemble_functions_nonlinear_elasticity_3D.cpp
 *
 *  Created on: Aug 28, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "assemble_functions_nonlinear_elasticity_3D.h"

libMesh::NonlinearImplicitSystem& add_nonlinear_elasticity(	libMesh::EquationSystems& input_systems,
														libMesh::Order order,
														libMesh::FEFamily family)
{
	libMesh::ExplicitSystem& physical_variables =
			input_systems.add_system<libMesh::ExplicitSystem> ("PhysicalConstants");

	// Physical constants are set as constant, monomial
	physical_variables.add_variable("E", libMesh::CONSTANT, libMesh::MONOMIAL);
	physical_variables.add_variable("mu", libMesh::CONSTANT, libMesh::MONOMIAL);

	libMesh::NonlinearImplicitSystem& nonlinear_elasticity_system =
			input_systems.add_system<libMesh::NonlinearImplicitSystem> ("NonlinearElasticity");

	nonlinear_elasticity_system.add_variable("u", order, family);
	nonlinear_elasticity_system.add_variable("v", order, family);
	nonlinear_elasticity_system.add_variable("w", order, family);

	return nonlinear_elasticity_system;
};

  /**
   * Evaluate the Jacobian of the nonlinear system.
   */
  void LargeDeformationElasticity::jacobian (const libMesh::NumericVector<libMesh::Number>& soln,
									  libMesh::SparseMatrix<libMesh::Number>&  jacobian,
									  libMesh::NonlinearImplicitSystem& /*sys*/)
  {

	libMesh::Number localE = -1;
	libMesh::Number localMu = -1;

	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
	const unsigned int young_var = physical_param_system.variable_number ("E");
	const unsigned int mu_var = physical_param_system.variable_number ("mu");

	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
	std::vector<libMesh::dof_id_type> physical_dof_indices_var;

    const libMesh::MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    libMesh::NonlinearImplicitSystem& system =
      es.get_system<libMesh::NonlinearImplicitSystem>("NonlinearElasticity");

    const unsigned int u_var = system.variable_number ("u");

    const libMesh::DofMap& dof_map = system.get_dof_map();

    libMesh::FEType fe_type = dof_map.variable_type(u_var);
    libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qface (dim-1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);

    const std::vector<libMesh::Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseSubMatrix<libMesh::Number> Ke_var[3][3] =
      {
        {libMesh::DenseSubMatrix<libMesh::Number>(Ke), libMesh::DenseSubMatrix<libMesh::Number>(Ke), libMesh::DenseSubMatrix<libMesh::Number>(Ke)},
        {libMesh::DenseSubMatrix<libMesh::Number>(Ke), libMesh::DenseSubMatrix<libMesh::Number>(Ke), libMesh::DenseSubMatrix<libMesh::Number>(Ke)},
        {libMesh::DenseSubMatrix<libMesh::Number>(Ke), libMesh::DenseSubMatrix<libMesh::Number>(Ke), libMesh::DenseSubMatrix<libMesh::Number>(Ke)}
      };

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector< std::vector<libMesh::dof_id_type> > dof_indices_var(3);

    jacobian.zero();

    libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	// Quadrature points
	const std::vector<libMesh::Point>& qp_points = fe->get_xyz();

	// Weights for the Arlequin method
	double alpha_micro = 1;

    for ( ; el != end_el; ++el)
      {
        const libMesh::Elem* elem = *el;
        dof_map.dof_indices (elem, dof_indices);
        for(unsigned int var=0; var<3; var++)
          {
            dof_map.dof_indices (elem, dof_indices_var[var], var);
          }

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_var_dofs = dof_indices_var[0].size();

		physical_dof_map.dof_indices(elem, physical_dof_indices_var, young_var);
		localE = physical_param_system.current_solution(physical_dof_indices_var[0]);

		physical_dof_map.dof_indices(elem, physical_dof_indices_var, mu_var);
		localMu = physical_param_system.current_solution(physical_dof_indices_var[0]);

        fe->reinit (elem);

        Ke.resize (n_dofs,n_dofs);
        for(unsigned int var_i=0; var_i<3; var_i++)
          for(unsigned int var_j=0; var_j<3; var_j++)
            {
              Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);
            }

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
			alpha_micro = weight_mask.get_alpha_micro(qp_points[qp]);

        	libMesh::DenseVector<libMesh::Number> u_vec(3);
        	libMesh::DenseMatrix<libMesh::Number> grad_u(3,3);
            for(unsigned int var_i=0; var_i<3; var_i++)
              {
                for (unsigned int j=0; j<n_var_dofs; j++)
                  {
                    u_vec(var_i) += phi[j][qp]*soln(dof_indices_var[var_i][j]);
                  }

                for(unsigned int var_j=0; var_j<3; var_j++)
                  {
                    for (unsigned int j=0; j<n_var_dofs; j++)
                      {
                        // Row is variable u1, u2, or u3, column is x, y, or z
                        grad_u(var_i,var_j) += dphi[j][qp](var_j)*soln(dof_indices_var[var_i][j]);
                      }
                  }
              }

            libMesh::DenseMatrix<libMesh::Number> strain_tensor(3,3);
            for(unsigned int i=0; i<3; i++)
              for(unsigned int j=0; j<3; j++)
                {
                  strain_tensor(i,j) += 0.5 * ( grad_u(i,j) + grad_u(j,i) );

                  for(unsigned int k=0; k<3; k++)
                    {
                      strain_tensor(i,j) += 0.5 * grad_u(k,i)*grad_u(k,j);
                    }
                }

            // Define the deformation gradient
            libMesh::DenseMatrix<libMesh::Number> F(3,3);
            F = grad_u;
            for(unsigned int var=0; var<3; var++)
              {
                F(var,var) += 1.;
              }

            libMesh::DenseMatrix<libMesh::Number> stress_tensor(3,3);

            for(unsigned int i=0; i<3; i++)
              for(unsigned int j=0; j<3; j++)
                for(unsigned int k=0; k<3; k++)
                  for(unsigned int l=0; l<3; l++)
                    {
                      stress_tensor(i,j) +=
                    		  alpha_micro*this->eval_elasticity_tensor(i,j,k,l,localE, localMu) *
                        strain_tensor(k,l);
                    }

            for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
                {
                  for(unsigned int i=0; i<3; i++)
                    for(unsigned int j=0; j<3; j++)
                      for(unsigned int m=0; m<3; m++)
                        {
                          Ke_var[i][i](dof_i,dof_j) += JxW[qp] *
                            ( -dphi[dof_j][qp](m) * stress_tensor(m,j) * dphi[dof_i][qp](j) );
                        }

                  for(unsigned int i=0; i<3; i++)
                    for(unsigned int j=0; j<3; j++)
                      for(unsigned int k=0; k<3; k++)
                        for(unsigned int l=0; l<3; l++)
                          {
                        	libMesh::Number FxC_ijkl = 0.;
                            for(unsigned int m=0; m<3; m++)
                              {
                                FxC_ijkl +=
                                  F(i,m) *
								  alpha_micro*this->eval_elasticity_tensor(m,j,k,l,localE, localMu);
                              }

                            Ke_var[i][k](dof_i,dof_j) += JxW[qp] *
                              ( -0.5 * FxC_ijkl *
                                dphi[dof_j][qp](l) *
                                dphi[dof_i][qp](j)
                                );

                            Ke_var[i][l](dof_i,dof_j) += JxW[qp] *
                              ( -0.5 * FxC_ijkl *
                                dphi[dof_j][qp](k) *
                                dphi[dof_i][qp](j)
                                );

                            for(unsigned int n=0; n<3; n++)
                              {
                                Ke_var[i][n](dof_i,dof_j) += JxW[qp] *
                                  ( -0.5 * FxC_ijkl *
                                    ( dphi[dof_j][qp](k) * grad_u(n,l) +
                                      dphi[dof_j][qp](l) * grad_u(n,k) ) *
                                    dphi[dof_i][qp](j)
                                    );
                              }

                          }
                }

          }
        jacobian.add_matrix (Ke, dof_indices);
      }
    jacobian.close();

    std::cout << jacobian.l1_norm() << " " << jacobian.linfty_norm() << std::endl;
    // Now add the new term!
    jacobian.add(1.,*m_Coupling_Term_Matrix);
    std::cout << jacobian.l1_norm() << " " << jacobian.linfty_norm() << std::endl;
  }

  /**
   * Evaluate the residual of the nonlinear system.
   */
  void LargeDeformationElasticity::residual (const libMesh::NumericVector<libMesh::Number>& soln,
		  libMesh::NumericVector<libMesh::Number>& residual,
		  libMesh::NonlinearImplicitSystem& /*sys*/)
  {
		libMesh::Number localE = -1;
		libMesh::Number localMu = -1;

		libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
		const unsigned int young_var = physical_param_system.variable_number ("E");
		const unsigned int mu_var = physical_param_system.variable_number ("mu");

		const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
		std::vector<libMesh::dof_id_type> physical_dof_indices_var;

    const libMesh::MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    libMesh::NonlinearImplicitSystem& system =
      es.get_system<libMesh::NonlinearImplicitSystem>("NonlinearElasticity");

    const unsigned int u_var = system.variable_number ("u");

    const libMesh::DofMap& dof_map = system.get_dof_map();

    libMesh::FEType fe_type = dof_map.variable_type(u_var);
    libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qface (dim-1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);

    const std::vector<libMesh::Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

    libMesh::DenseVector<libMesh::Number> Re;
    libMesh::DenseSubVector<libMesh::Number> Re_var[3] =
      {libMesh::DenseSubVector<libMesh::Number>(Re), libMesh::DenseSubVector<libMesh::Number>(Re), libMesh::DenseSubVector<libMesh::Number>(Re)};

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector< std::vector<libMesh::dof_id_type> > dof_indices_var(3);

    residual.zero();

    libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // Quadrature points
    const std::vector<libMesh::Point>& qp_points = fe->get_xyz();

	// Weights for the Arlequin method
	double alpha_micro = 1;

    for ( ; el != end_el; ++el)
      {
        const libMesh::Elem* elem = *el;
        dof_map.dof_indices (elem, dof_indices);
        for(unsigned int var=0; var<3; var++)
          {
            dof_map.dof_indices (elem, dof_indices_var[var], var);
          }

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_var_dofs = dof_indices_var[0].size();

		physical_dof_map.dof_indices(elem, physical_dof_indices_var, young_var);
		localE = physical_param_system.current_solution(physical_dof_indices_var[0]);

		physical_dof_map.dof_indices(elem, physical_dof_indices_var, mu_var);
		localMu = physical_param_system.current_solution(physical_dof_indices_var[0]);

        fe->reinit (elem);

        Re.resize (n_dofs);
        for(unsigned int var=0; var<3; var++)
          {
            Re_var[var].reposition (var*n_var_dofs, n_var_dofs);
          }

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {

			alpha_micro = weight_mask.get_alpha_micro(qp_points[qp]);

        	libMesh::DenseVector<libMesh::Number> u_vec(3);
        	libMesh::DenseMatrix<libMesh::Number> grad_u(3,3);
            for(unsigned int var_i=0; var_i<3; var_i++)
              {
                for (unsigned int j=0; j<n_var_dofs; j++)
                  {
                    u_vec(var_i) += phi[j][qp]*soln(dof_indices_var[var_i][j]);
                  }

                for(unsigned int var_j=0; var_j<3; var_j++)
                  {
                    for (unsigned int j=0; j<n_var_dofs; j++)
                      {
                        // Row is variable u, v, or w column is x, y, or z
                        grad_u(var_i,var_j) += dphi[j][qp](var_j)*soln(dof_indices_var[var_i][j]);
                      }
                  }
              }

            libMesh::DenseMatrix<libMesh::Number> strain_tensor(3,3);
            for(unsigned int i=0; i<3; i++)
              for(unsigned int j=0; j<3; j++)
                {
                  strain_tensor(i,j) += 0.5 * ( grad_u(i,j) + grad_u(j,i) );

                  for(unsigned int k=0; k<3; k++)
                    {
                      strain_tensor(i,j) += 0.5 * grad_u(k,i)*grad_u(k,j);
                    }
                }

            // Define the deformation gradient
            libMesh::DenseMatrix<libMesh::Number> F(3,3);
            F = grad_u;
            for(unsigned int var=0; var<3; var++)
            {
              F(var,var) += 1.;
            }

            libMesh::DenseMatrix<libMesh::Number> stress_tensor(3,3);

            for(unsigned int i=0; i<3; i++)
              for(unsigned int j=0; j<3; j++)
                for(unsigned int k=0; k<3; k++)
                  for(unsigned int l=0; l<3; l++)
                    {
                      stress_tensor(i,j) +=
                    		  alpha_micro*this->eval_elasticity_tensor(i,j,k,l,localE, localMu) *
                        strain_tensor(k,l);
                    }

            for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              {
                for(unsigned int i=0; i<3; i++)
                  {
                    for(unsigned int j=0; j<3; j++)
                      {
                    	libMesh::Number FxStress_ij = 0.;
                        for(unsigned int m=0; m<3; m++)
                          {
                            FxStress_ij += F(i,m) * stress_tensor(m,j);
                          }

                        Re_var[i](dof_i) += JxW[qp] *
                          ( -FxStress_ij * dphi[dof_i][qp](j) );
                      }
                  }
              }

          }

        dof_map.constrain_element_vector (Re, dof_indices);
        residual.add_vector (Re, dof_indices);
      }
    residual.close();
    std::cout << residual.l1_norm() << " " << residual.l2_norm() << " " << residual.linfty_norm() << std::endl;
    // Now add the new term
    residual.add(*m_Coupling_Term_Vector);
    std::cout << residual.l1_norm() << " " << residual.l2_norm() << " " << residual.linfty_norm() << std::endl;
  }

  /**
   * Compute the Cauchy stress for the current solution.
   */
  void LargeDeformationElasticity::compute_stresses()
  {
		libMesh::Number localE = -1;
		libMesh::Number localMu = -1;

		libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
		const unsigned int young_var = physical_param_system.variable_number ("E");
		const unsigned int mu_var = physical_param_system.variable_number ("mu");

		const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
		std::vector<libMesh::dof_id_type> physical_dof_indices_var;

    const libMesh::MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    libMesh::NonlinearImplicitSystem& system =
      es.get_system<libMesh::NonlinearImplicitSystem>("NonlinearElasticity");

    unsigned int displacement_vars[3];
    displacement_vars[0] = system.variable_number ("u");
    displacement_vars[1] = system.variable_number ("v");
    displacement_vars[2] = system.variable_number ("w");
    const unsigned int u_var = system.variable_number ("u");

    const libMesh::DofMap& dof_map = system.get_dof_map();
    libMesh::FEType fe_type = dof_map.variable_type(u_var);
    libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    const std::vector<libMesh::Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

    // Also, get a reference to the ExplicitSystem
    libMesh::ExplicitSystem& stress_system = es.get_system<libMesh::ExplicitSystem>("StressSystem");
    const libMesh::DofMap& stress_dof_map = stress_system.get_dof_map();
    unsigned int sigma_vars[6];
    sigma_vars[0] = stress_system.variable_number ("sigma_00");
    sigma_vars[1] = stress_system.variable_number ("sigma_01");
    sigma_vars[2] = stress_system.variable_number ("sigma_02");
    sigma_vars[3] = stress_system.variable_number ("sigma_11");
    sigma_vars[4] = stress_system.variable_number ("sigma_12");
    sigma_vars[5] = stress_system.variable_number ("sigma_22");
	unsigned int vonMises_var = stress_system.variable_number ("vonMises");

    // Storage for the stress dof indices on each element
    std::vector< std::vector<libMesh::dof_id_type> > dof_indices_var(system.n_vars());
    std::vector<libMesh::dof_id_type> stress_dof_indices_var;

    // To store the stress tensor on each element
    libMesh::DenseMatrix<libMesh::Number> elem_avg_stress_tensor(3,3);

    libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
      {
        const libMesh::Elem* elem = *el;

        for(unsigned int var=0; var<3; var++)
          {
            dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);
          }

        const unsigned int n_var_dofs = dof_indices_var[0].size();

		physical_dof_map.dof_indices(elem, physical_dof_indices_var, young_var);
		localE = physical_param_system.current_solution(physical_dof_indices_var[0]);

		physical_dof_map.dof_indices(elem, physical_dof_indices_var, mu_var);
		localMu = physical_param_system.current_solution(physical_dof_indices_var[0]);

        fe->reinit (elem);

        // clear the stress tensor
        elem_avg_stress_tensor.resize(3,3);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
        	libMesh::DenseMatrix<libMesh::Number> grad_u(3,3);
            for(unsigned int var_i=0; var_i<3; var_i++)
              for(unsigned int var_j=0; var_j<3; var_j++)
                {
                  for (unsigned int j=0; j<n_var_dofs; j++)
                    {
                      // Row is variable u1, u2, or u3, column is x, y, or z
                      grad_u(var_i,var_j) += dphi[j][qp](var_j)*
                        system.current_solution(dof_indices_var[var_i][j]);
                    }
                }

            libMesh::DenseMatrix<libMesh::Number> strain_tensor(3,3);
            for(unsigned int i=0; i<3; i++)
              for(unsigned int j=0; j<3; j++)
                {
                  strain_tensor(i,j) += 0.5 * ( grad_u(i,j) + grad_u(j,i) );

                  for(unsigned int k=0; k<3; k++)
                    {
                      strain_tensor(i,j) += 0.5 * grad_u(k,i)*grad_u(k,j);
                    }
                }

            // Define the deformation gradient
            libMesh::DenseMatrix<libMesh::Number> F(3,3);
            F = grad_u;
            for(unsigned int var=0; var<3; var++)
            {
              F(var,var) += 1.;
            }

            libMesh::DenseMatrix<libMesh::Number> stress_tensor(3,3);
            for(unsigned int i=0; i<3; i++)
              for(unsigned int j=0; j<3; j++)
                for(unsigned int k=0; k<3; k++)
                  for(unsigned int l=0; l<3; l++)
                    {
                      stress_tensor(i,j) +=
                        this->eval_elasticity_tensor(i,j,k,l,localE, localMu) *
                        strain_tensor(k,l);
                    }

            // stress_tensor now holds the second Piola-Kirchoff stress (PK2) at point qp.
            // However, in this example we want to compute the Cauchy stress which is given by
            // 1/det(F) * F * PK2 * F^t, hence we now apply this transformation.
            stress_tensor.scale(1./F.det());
            stress_tensor.left_multiply(F);
            stress_tensor.right_multiply_transpose(F);

            // We want to plot the average Cauchy stress on each element, hence
            // we integrate stress_tensor
            elem_avg_stress_tensor.add(JxW[qp], stress_tensor);
          }

        // Get the average stress per element by dividing by volume
        elem_avg_stress_tensor.scale(1./elem->volume());

        // load elem_sigma data into stress_system
        unsigned int stress_var_index = 0;
        for(unsigned int i=0; i<3; i++)
          for(unsigned int j=i; j<3; j++)
            {
              stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[stress_var_index]);

              // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
              // one dof index per variable
              libMesh::dof_id_type dof_index = stress_dof_indices_var[0];

              if( (stress_system.solution->first_local_index() <= dof_index) &&
                  (dof_index < stress_system.solution->last_local_index()) )
                {
                  stress_system.solution->set(dof_index, elem_avg_stress_tensor(i,j));
                }

              stress_var_index++;
            }

    	// Calculate von Mises
    	libMesh::Number vonMises_value = std::sqrt( 0.5*( pow(elem_avg_stress_tensor(0,0) - elem_avg_stress_tensor(1,1),2.) +
    			 pow(elem_avg_stress_tensor(1,1) - elem_avg_stress_tensor(2,2),2.) +
    			 pow(elem_avg_stress_tensor(2,2) - elem_avg_stress_tensor(0,0),2.) +
    			 6.*(pow(elem_avg_stress_tensor(0,1),2.) + pow(elem_avg_stress_tensor(1,2),2.) + pow(elem_avg_stress_tensor(2,0),2.))
    			 ) );

    	// Get the DoF map
    	stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);

    	// Get the first index (CONSTANT MONOMIAL basis functions, hence
    	// only one element
    	libMesh::dof_id_type dof_index = stress_dof_indices_var[0];

    	// To be sure, test if the DoF index is inside this processor
    	if( (stress_system.solution->first_local_index() <= dof_index) &&
    		(dof_index < stress_system.solution->last_local_index()) )
    	{
    		stress_system.solution->set(dof_index, vonMises_value);
    	}
      }



    // Should call close and update when we set vector entries directly
    stress_system.solution->close();
    stress_system.update();
  }
