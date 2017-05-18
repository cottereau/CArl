/*
 * non_linear_FEMSystem.cpp
 *
 *  Created on: Jul 21, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "non_linear_FEMSystem.h"

carl::NonLinearSystem::NonLinearSystem(libMesh::EquationSystems& es, const std::string& name_in,
                         const unsigned int number_in) :
						 libMesh::FEMSystem(es, name_in, number_in) {

	// Add a time solver. We are just looking at a steady state problem here.
	this->time_solver = libMesh::UniquePtr<libMesh::TimeSolver>(new libMesh::SteadySolver(*this));
}

void carl::NonLinearSystem::save_initial_mesh() {
	System & aux_sys = this->get_equation_systems().get_system("auxiliary");

	const unsigned int dim = this->get_mesh().mesh_dimension();

	// Loop over all nodes and copy the location from the current system to
	// the auxiliary system.
	const libMesh::MeshBase::const_node_iterator nd_end =
	this->get_mesh().local_nodes_end();
	for (libMesh::MeshBase::const_node_iterator nd = this->get_mesh().local_nodes_begin();
		nd != nd_end; ++nd) {
		const libMesh::Node *node = *nd;
		for (unsigned int d = 0; d < dim; ++d) {
			unsigned int source_dof = node->dof_number(this->number(), var[d], 0);
			unsigned int dest_dof = node->dof_number(aux_sys.number(), undefo_var[d],
												   0);
			libMesh::Number value = this->current_local_solution->el(source_dof);
			aux_sys.current_local_solution->set(dest_dof, value);
		}
	}
}

void carl::NonLinearSystem::init_data() {
	const unsigned int dim = this->get_mesh().mesh_dimension();

	// Get the default order of the used elements. Assumption:
	// Just one type of elements in the mesh.
	libMesh::Order order = (*(this->get_mesh().elements_begin()))->default_order();

	// Add the node positions as primary variables.
	var[0] = this->add_variable("x", order);
	var[1] = this->add_variable("y", order);
	if (dim == 3)
	var[2] = this->add_variable("z", order);
	else
	var[2] = var[1];

	// Add variables for storing the initial mesh to an auxiliary system.
	System& aux_sys = this->get_equation_systems().get_system("auxiliary");
	undefo_var[0] = aux_sys.add_variable("undefo_x", order);
	undefo_var[1] = aux_sys.add_variable("undefo_y", order);
	undefo_var[2] = aux_sys.add_variable("undefo_z", order);

	// Set the time stepping options
	this->deltat = args("schedule/dt", 0.2);

	// Do the parent's initialization after variables are defined
	FEMSystem::init_data();

	//// Tell the system to march velocity forward in time, but
	//// leave p as a constraint only
	//this->time_evolving(var[0]);
	//this->time_evolving(var[1]);
	//if (dim == 3)
	//this->time_evolving(var[2]);

	// Tell the system which variables are containing the node positions
	set_mesh_system(this);

	this->set_mesh_x_var(var[0]);
	this->set_mesh_y_var(var[1]);
	if (dim == 3)
	this->set_mesh_z_var(var[2]);

	// Fill the variables with the position of the nodes
	this->mesh_position_get();

	System::reinit();

	// Set some options for the DiffSolver
	libMesh::DiffSolver &solver = *(this->time_solver->diff_solver().get());
	solver.quiet = args("solver/quiet", false);
	solver.max_nonlinear_iterations = args(
										 "solver/nonlinear/max_nonlinear_iterations", 100);
	solver.relative_step_tolerance = args(
										"solver/nonlinear/relative_step_tolerance", 1.e-3);
	solver.relative_residual_tolerance = args(
											"solver/nonlinear/relative_residual_tolerance", 1.e-8);
	solver.absolute_residual_tolerance = args(
											"solver/nonlinear/absolute_residual_tolerance", 1.e-8);
	solver.verbose = !args("solver/quiet", false);

	((libMesh::NewtonSolver&) solver).require_residual_reduction = args(
															 "solver/nonlinear/require_reduction", false);

	// And the linear solver options
	solver.max_linear_iterations = args("max_linear_iterations", 50000);
	solver.initial_linear_tolerance = args("initial_linear_tolerance", 1.e-3);
}

void carl::NonLinearSystem::update() {
	System::update();
	this->mesh_position_set();
}

void carl::NonLinearSystem::init_context(libMesh::DiffContext &context) {
	libMesh::FEMContext &c = libMesh::cast_ref<libMesh::FEMContext&>(context);

	// Pre-request all the data needed
	libMesh::FEBase* elem_fe = NULL;
	c.get_element_fe( 0, elem_fe );

	elem_fe->get_JxW();
	elem_fe->get_phi();
	elem_fe->get_dphi();
	elem_fe->get_xyz();

	libMesh::FEBase* side_fe = NULL;
	c.get_side_fe( 0, side_fe );

	side_fe->get_JxW();
	side_fe->get_phi();
	side_fe->get_xyz();
}

/**
 * Compute contribution from internal forces in elem_residual and contribution from
 * linearized internal forces (stiffness matrix) in elem_jacobian.
 */
bool carl::NonLinearSystem::element_time_derivative(bool request_jacobian,
	libMesh::DiffContext &context) {
	libMesh::FEMContext &c = libMesh::cast_ref<libMesh::FEMContext&>(context);

	// First we get some references to cell-specific data that
	// will be used to assemble the linear system.
	libMesh::FEBase* elem_fe = NULL;
	c.get_element_fe( 0, elem_fe );

	// Element Jacobian * quadrature weights for interior integration
	const std::vector<libMesh::Real> &JxW = elem_fe->get_JxW();

	// Element basis functions
	const std::vector<std::vector<libMesh::RealGradient> > &dphi = elem_fe->get_dphi();

	// Dimension of the mesh
	const unsigned int dim = this->get_mesh().mesh_dimension();

	// The number of local degrees of freedom in each variable
	const unsigned int n_u_dofs = c.get_dof_indices(var[0]).size();
	libmesh_assert(n_u_dofs ==  c.get_dof_indices(var[1]).size());
	if (dim == 3) {
	libmesh_assert(n_u_dofs ==  c.get_dof_indices(var[2]).size());
	}

	unsigned int n_qpoints = c.get_element_qrule().n_points();

	// Some matrices and vectors for storing the results of the constitutive
	// law
	libMesh::DenseMatrix<libMesh::Real> stiff;
	libMesh::DenseVector<libMesh::Real> res;
	libMesh::VectorValue<libMesh::Gradient> grad_u;

	// Instantiate the constitutive law
	// Just calculate jacobian contribution when we need to
	NonlinearNeoHookeCurrentConfig material(dphi, args, request_jacobian);

	// Get a reference to the auxiliary system
	libMesh::TransientExplicitSystem& aux_system = this->get_equation_systems().get_system<
			libMesh::TransientExplicitSystem>("auxiliary");
	std::vector<libMesh::dof_id_type> undefo_index;

	// Assume symmetry of local stiffness matrices
	bool use_symmetry = args("assembly/use_symmetry", false);

	// Now we will build the element Jacobian and residual.
	// This must be calculated at each quadrature point by
	// summing the solution degree-of-freedom values by
	// the appropriate weight functions.
	// This class just takes care of the assembly. The matrix of
	// the jacobian and the residual vector are provided by the
	// constitutive formulation.

	for (unsigned int qp = 0; qp != n_qpoints; qp++)
	{
		// Compute the displacement gradient
		grad_u(0) = grad_u(1) = grad_u(2) = 0;
		for (unsigned int d = 0; d < dim; ++d)
		{
			std::vector<libMesh::Number> u_undefo;
			aux_system.get_dof_map().dof_indices(&c.get_elem(), undefo_index, undefo_var[d]);
			aux_system.current_local_solution->get(undefo_index, u_undefo);
			for (unsigned int l = 0; l != n_u_dofs; l++)
			{
				grad_u(d).add_scaled(dphi[l][qp], u_undefo[l]); // u_current(l)); // -
			}
		}

		// initialize the constitutive formulation with the current displacement
		// gradient
		material.init_for_qp(grad_u, qp);

		// Aquire, scale and assemble residual and stiffness
		for (unsigned int i = 0; i < n_u_dofs; i++)
		{
			res.resize(dim);
			material.get_residual(res, i);
			res.scale(JxW[qp]);
			for (unsigned int ii = 0; ii < dim; ++ii) {
				(c.get_elem_residual(ii))(i) += res(ii);
			}

			if (request_jacobian && c.elem_solution_derivative)
			{
				libmesh_assert(c.elem_solution_derivative == 1.0);
				for (unsigned int j = (use_symmetry ? i : 0); j < n_u_dofs; j++)
				{
					material.get_linearized_stiffness(stiff, i, j);
					stiff.scale(JxW[qp]);
					for (unsigned int ii = 0; ii < dim; ++ii)
					{
						for (unsigned int jj = 0; jj < dim; ++jj)
						{
							(c.get_elem_jacobian(ii,jj))(i, j) += stiff(ii, jj);
							if (use_symmetry && i != j)
							{
								(c.get_elem_jacobian(ii,jj))(j, i) += stiff(jj, ii);
							}
						}
					}
				}
			}
		}
	} // end of the quadrature point qp-loop

	return request_jacobian;
}

