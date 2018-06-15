#include "common_assemble_functions_elasticity_3D.h"
// Some boundary conditions functions
void set_displaced_border_translation(libMesh::ImplicitSystem& elasticity_system, border_displacement_values& displ, int boundary_id)
{
	// Defining the boundaries with Dirichlet conditions ...
	std::set<libMesh::boundary_id_type> boundary_id_displacement;

	boundary_id_displacement.insert(boundary_id);

	std::vector<unsigned int> variables(3);
	variables[0] = elasticity_system.variable_number("u");
	variables[1] = elasticity_system.variable_number("v");
	variables[2] = elasticity_system.variable_number("w");

	border_displacement_function move_border(	variables[0],variables[1],variables[2],
										displ.x_displ,displ.y_displ,displ.z_displ);

	// ... and set them
	libMesh::DirichletBoundary dirichlet_bc(	boundary_id_displacement,
												variables,
												&move_border);

	elasticity_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
}

void set_clamped_border(libMesh::ImplicitSystem& elasticity_system, int boundary_id)
{
	// Defining the boundaries with Dirichlet conditions ...
	std::set<libMesh::boundary_id_type> boundary_id_displacement;

	boundary_id_displacement.insert(boundary_id);

	std::vector<unsigned int> variables(3);
	variables[0] = elasticity_system.variable_number("u");
	variables[1] = elasticity_system.variable_number("v");
	variables[2] = elasticity_system.variable_number("w");

	libMesh::ZeroFunction<> zero_function;

	// ... and set them
	libMesh::DirichletBoundary dirichlet_bc(	boundary_id_displacement,
												variables,
												&zero_function);

	elasticity_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
}

libMesh::ExplicitSystem& add_stress(libMesh::EquationSystems& input_systems)
{
	libMesh::ExplicitSystem& stress_system =
			input_systems.add_system<libMesh::ExplicitSystem> ("StressSystem");

	stress_system.add_variable("sigma_00", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_01", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_02", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_10", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_11", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_12", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_20", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_21", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("sigma_22", libMesh::CONSTANT, libMesh::MONOMIAL);
	stress_system.add_variable("vonMises", libMesh::CONSTANT, libMesh::MONOMIAL);

	return stress_system;
}

libMesh::LinearImplicitSystem& add_elasticity(	libMesh::EquationSystems& input_systems,
												libMesh::Order order,
												libMesh::FEFamily family)
{
	libMesh::ExplicitSystem& physical_variables =
			input_systems.add_system<libMesh::ExplicitSystem> ("PhysicalConstants");

	// Physical constants are set as constant, monomial
	physical_variables.add_variable("E", libMesh::CONSTANT, libMesh::MONOMIAL);
	physical_variables.add_variable("mu", libMesh::CONSTANT, libMesh::MONOMIAL);
	physical_variables.add_variable("Index", libMesh::CONSTANT, libMesh::MONOMIAL);

	libMesh::LinearImplicitSystem& elasticity_system =
			input_systems.add_system<libMesh::LinearImplicitSystem> ("Elasticity");

	elasticity_system.add_variable("u", order, family);
	elasticity_system.add_variable("v", order, family);
	elasticity_system.add_variable("w", order, family);

	return elasticity_system;
}

libMesh::ExplicitSystem& add_vel_newmark(libMesh::EquationSystems& input_systems, 
    libMesh::Order order, 
    libMesh::FEFamily family)
{
    libMesh::ExplicitSystem& v_system = 
        input_systems.add_system<libMesh::ExplicitSystem>("Velocity");
    v_system.add_variable("u_vel", order, family);
    v_system.add_variable("v_vel", order, family);
    v_system.add_variable("w_vel", order, family);

    return v_system;
}

libMesh::ExplicitSystem& add_acc_newmark(libMesh::EquationSystems& input_systems, 
    libMesh::Order order,
    libMesh::FEFamily family)
{
    libMesh::ExplicitSystem& a_system = 
        input_systems.add_system<libMesh::ExplicitSystem>("Acceleration");
    a_system.add_variable("u_accel", order, family);
    a_system.add_variable("v_accel", order, family);
    a_system.add_variable("w_accel", order, family);

    return a_system;
}


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=4 et tw=80 smartindent :                               */
