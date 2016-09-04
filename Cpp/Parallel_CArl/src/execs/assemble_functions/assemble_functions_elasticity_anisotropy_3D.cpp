/*
 * assemble_functions_elasticity_anisotropy.cpp
 *
 *  Created on: Sep 3, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "assemble_functions_elasticity_anisotropy_3D.h"

void Update_SubK(	libMesh::DenseSubMatrix<libMesh::Number>& SubK,
					unsigned int qp,
					unsigned int C_i,
					unsigned int C_j,
					const std::vector<std::vector<libMesh::RealGradient> >& dphi,
					const unsigned int n_components,
					const unsigned int n_u_dofs,
					const std::vector<libMesh::Real>& JxW,
					carl::anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input,
					int grain_idx,
					double cte
					)
{
	for (unsigned int iii=0; iii<n_u_dofs; iii++)
	{
		for (unsigned int jjj=0; jjj<n_u_dofs; jjj++)
		{
			for(unsigned int C_k=0; C_k<n_components; C_k++)
			{
				for(unsigned int C_l=0; C_l<n_components; C_l++)
				{
					SubK(iii,jjj) += cte * JxW[qp]*(anisotropy_obj_input.eval_elasticity_tensor(C_i,C_j,C_k,C_l,grain_idx) * dphi[iii][qp](C_k)*dphi[jjj][qp](C_l));
				}
			}
		}
	}
};

void assemble_elasticity_anisotropic_with_weight(	libMesh::EquationSystems& es,
							const std::string& system_name,
							carl::weight_parameter_function& weight_mask,
							carl::anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input)
{
	libmesh_assert_equal_to (system_name, "Elasticity");

	libMesh::PerfLog perf_log ("Matrix Assembly (Heterogeneous elasticity)",MASTER_bPerfLog_assemble_fem);

	perf_log.push("Preamble");
	// Set up mesh
	const libMesh::MeshBase& mesh = es.get_mesh();

	const unsigned int dim = mesh.mesh_dimension();

	// - Set up physical properties system ------------------------------------
	int localIdx;

	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
	const unsigned int idx_var = physical_param_system.variable_number ("Index");

	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();
	std::vector<libMesh::dof_id_type> physical_dof_indices_var;

	// - Set up elasticity system ---------------------------------------------
	libMesh::LinearImplicitSystem& system = es.get_system<libMesh::LinearImplicitSystem>("Elasticity");

	const unsigned int n_components = 3;
	const unsigned int u_var = system.variable_number ("u");
	const unsigned int v_var = system.variable_number ("v");
	const unsigned int w_var = system.variable_number ("w");

	const libMesh::DofMap& dof_map = system.get_dof_map();
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

	// Shape functions
	const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();

	// Shape functions derivatives
	const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

	// Quadrature points
	const std::vector<libMesh::Point>& qp_points = fe->get_xyz();

	// Weights for the Arlequin method
	double alpha_micro = 1;

	libMesh::DenseMatrix<libMesh::Number> Ke;
	libMesh::DenseVector<libMesh::Number> Fe;

	libMesh::DenseSubMatrix<libMesh::Number>
	Kuu(Ke), Kuv(Ke), Kuw(Ke),
	Kvu(Ke), Kvv(Ke), Kvw(Ke),
	Kwu(Ke), Kwv(Ke), Kww(Ke);

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

		perf_log.push("Define physical params");
		// The DoF and values of the physical system
		physical_dof_map.dof_indices(elem, physical_dof_indices_var, idx_var);
		localIdx = physical_param_system.current_solution(physical_dof_indices_var[0]);
		perf_log.pop("Define physical params");

		// Restart the FE to the "geometry" of the element
		// -> Determines quadrature points, shape functions ...
		// !!! User can change the points to be used (can use other mesh's points
		//		instead of the quadrature points)
		fe->reinit (elem);

		perf_log.push("Matrix manipulations");
		Ke.resize (n_dofs, n_dofs);
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

		Fu.reposition (u_var*n_u_dofs, n_u_dofs);
		Fv.reposition (v_var*n_u_dofs, n_v_dofs);
		Fw.reposition (w_var*n_u_dofs, n_w_dofs);
		perf_log.push("Matrix manipulations");

		perf_log.push("Matrix element calculations");
		// For each quadrature point determinate the sub-matrices elements
		for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		{

			perf_log.push("Rigidity","Matrix element calculations");
			// Internal tension
			alpha_micro = weight_mask.get_alpha_micro(qp_points[qp]);

			// Internal tension
			std::cout << localIdx << std::endl;
			Update_SubK(Kuu, qp, 0, 0, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);
			Update_SubK(Kuv, qp, 0, 1, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);
			Update_SubK(Kuw, qp, 0, 2, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);

			Update_SubK(Kvu, qp, 1, 0, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);
			Update_SubK(Kvv, qp, 1, 1, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);
			Update_SubK(Kvw, qp, 1, 2, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);

			Update_SubK(Kwu, qp, 2, 0, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);
			Update_SubK(Kwv, qp, 2, 1, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);
			Update_SubK(Kww, qp, 2, 2, dphi, n_components, n_u_dofs, JxW, anisotropy_obj_input, localIdx, alpha_micro);

			perf_log.pop("Rigidity","Matrix element calculations");
			// Gravity
			//		if(z_force)
			//		{
			//			for (unsigned int i=0; i<n_w_dofs; i++)
			//			  {
			//				Fw(i) -= JxW[qp] * phi[i][qp];
			//			  }
			//		}
		}

		// Apply constraints
		perf_log.push("Constraints","Matrix element calculations");
		dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
		perf_log.pop("Constraints","Matrix element calculations");

		perf_log.push("Adding elements");
		system.matrix->add_matrix (Ke, dof_indices);
		system.rhs->add_vector    (Fe, dof_indices);
		perf_log.pop("Adding elements");
	}
}
