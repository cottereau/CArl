#include "common_assemble_functions_elasticity_3D.h"
// Create Elasticity System for dynamic analysis
void ElasticitySystem::init_data()
{
    _u_var = this->add_variable ("Ux", libMesh::FIRST, libMesh::LAGRANGE);
    _v_var = this->add_variable ("Uy", libMesh::FIRST, libMesh::LAGRANGE);
    _w_var = this->add_variable ("Uz", libMesh::FIRST, libMesh::LAGRANGE);

    this->time_evolving(_u_var,2);
    this->time_evolving(_v_var,2);
    this->time_evolving(_w_var,2);

    std::set<libMesh::boundary_id_type> boundary_ids;
    boundary_ids.insert(BOUNDARY_ID_MIN_X);
    boundary_ids.insert(NODE_BOUNDARY_ID);
    boundary_ids.insert(EDGE_BOUNDARY_ID);

    std::vector<unsigned int> variables;
    variables.push_back(_u_var);
    variables.push_back(_v_var);
    variables.push_back(_w_var);

    libMesh::ZeroFunction<> zf;

    // Most DirichletBoundary users will want to supply a "locally
    // indexed" functor
    libMesh::DirichletBoundary dirichlet_bc(boundary_ids, variables, zf,
            libMesh::LOCAL_VARIABLE_ORDER);

    this->get_dof_map().add_dirichlet_boundary(dirichlet_bc);

    // Do the parent's initialization after variables and boundary constraints are defined
    libMesh::FEMSystem::init_data();
}

void ElasticitySystem::init_context(libMesh::DiffContext & context)
{
    libMesh::FEMContext & c = libMesh::cast_ref<libMesh::FEMContext &>(context);

    libMesh::FEBase * u_elem_fe;
    libMesh::FEBase * u_side_fe;

    c.get_element_fe(_u_var, u_elem_fe);
    c.get_side_fe(_u_var, u_side_fe);

    // We should prerequest all the data
    // we will need to build the residuals.
    u_elem_fe->get_JxW();
    u_elem_fe->get_phi();
    u_elem_fe->get_dphi();

    u_side_fe->get_JxW();
    u_side_fe->get_phi();
}

bool ElasticitySystem::element_time_derivative(bool request_jacobian,
        libMesh::DiffContext & context)
{
    libMesh::FEMContext & c = libMesh::cast_ref<libMesh::FEMContext &>(context);

    // If we have an unsteady solver, then we need to extract the corresponding
    // velocity variable. This allows us to use either a FirstOrderUnsteadySolver
    // or a SecondOrderUnsteadySolver. That is, we get back the velocity variable
    // index for FirstOrderUnsteadySolvers or, if it's a SecondOrderUnsteadySolver,
    // this is actually just giving us back the same variable index.

    // If we only wanted to use a SecondOrderUnsteadySolver, then this
    // step would be unnecessary and we would just
    // populate the _u_var, etc. blocks of the residual and Jacobian.
    unsigned int u_dot_var = this->get_second_order_dot_var(_u_var);
    unsigned int v_dot_var = this->get_second_order_dot_var(_v_var);
    unsigned int w_dot_var = this->get_second_order_dot_var(_w_var);

    libMesh::FEBase * u_elem_fe;
    c.get_element_fe(_u_var, u_elem_fe);

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> & JxW = u_elem_fe->get_JxW();

    const std::vector<std::vector<libMesh::Real>> & phi = u_elem_fe->get_phi();
    const std::vector<std::vector<libMesh::RealGradient>> & grad_phi = u_elem_fe->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> & Fu = c.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fv = c.get_elem_residual(v_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fw = c.get_elem_residual(w_dot_var);

    libMesh::DenseSubMatrix<libMesh::Number> & Kuu = c.get_elem_jacobian(u_dot_var, _u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvv = c.get_elem_jacobian(v_dot_var, _v_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kww = c.get_elem_jacobian(w_dot_var, _w_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kuv = c.get_elem_jacobian(u_dot_var, _v_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kuw = c.get_elem_jacobian(u_dot_var, _w_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvu = c.get_elem_jacobian(v_dot_var, _u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvw = c.get_elem_jacobian(v_dot_var, _w_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kwu = c.get_elem_jacobian(w_dot_var, _u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kwv = c.get_elem_jacobian(w_dot_var, _v_var);

    unsigned int n_qpoints = c.get_element_qrule().n_points();

    libMesh::Gradient body_force(0.0, 0.0, -1.0);

    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        libMesh::Gradient grad_u, grad_v, grad_w;
        c.interior_gradient(_u_var, qp, grad_u);
        c.interior_gradient(_v_var, qp, grad_v);
        c.interior_gradient(_w_var, qp, grad_w);

        // Convenience
        libMesh::Tensor grad_U (grad_u, grad_v, grad_w);

        libMesh::Tensor tau;
        for (unsigned int i = 0; i < 3; i++)
            for (unsigned int j = 0; j < 3; j++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int l = 0; l < 3; l++)
                        tau(i,j) += elasticity_tensor(i,j,k,l)*grad_U(k,l);

        for (unsigned int i=0; i != n_u_dofs; i++)
        {
            for (unsigned int alpha = 0; alpha < 3; alpha++)
            {
                Fu(i) += (tau(0,alpha)*grad_phi[i][qp](alpha) - body_force(0)*phi[i][qp])*JxW[qp];
                Fv(i) += (tau(1,alpha)*grad_phi[i][qp](alpha) - body_force(1)*phi[i][qp])*JxW[qp];
                Fw(i) += (tau(2,alpha)*grad_phi[i][qp](alpha) - body_force(2)*phi[i][qp])*JxW[qp];

                if (request_jacobian)
                {
                    for (unsigned int j=0; j != n_u_dofs; j++)
                    {
                        for (unsigned int beta = 0; beta < 3; beta++)
                        {
                            // Convenience
                            const libMesh::Real c0 = grad_phi[j][qp](beta)*c.get_elem_solution_derivative();

                            libMesh::Real dtau_uu = elasticity_tensor(0, alpha, 0, beta)*c0;
                            libMesh::Real dtau_uv = elasticity_tensor(0, alpha, 1, beta)*c0;
                            libMesh::Real dtau_uw = elasticity_tensor(0, alpha, 2, beta)*c0;
                            libMesh::Real dtau_vu = elasticity_tensor(1, alpha, 0, beta)*c0;
                            libMesh::Real dtau_vv = elasticity_tensor(1, alpha, 1, beta)*c0;
                            libMesh::Real dtau_vw = elasticity_tensor(1, alpha, 2, beta)*c0;
                            libMesh::Real dtau_wu = elasticity_tensor(2, alpha, 0, beta)*c0;
                            libMesh::Real dtau_wv = elasticity_tensor(2, alpha, 1, beta)*c0;
                            libMesh::Real dtau_ww = elasticity_tensor(2, alpha, 2, beta)*c0;

                            Kuu(i,j) += dtau_uu*grad_phi[i][qp](alpha)*JxW[qp];
                            Kuv(i,j) += dtau_uv*grad_phi[i][qp](alpha)*JxW[qp];
                            Kuw(i,j) += dtau_uw*grad_phi[i][qp](alpha)*JxW[qp];
                            Kvu(i,j) += dtau_vu*grad_phi[i][qp](alpha)*JxW[qp];
                            Kvv(i,j) += dtau_vv*grad_phi[i][qp](alpha)*JxW[qp];
                            Kvw(i,j) += dtau_vw*grad_phi[i][qp](alpha)*JxW[qp];
                            Kwu(i,j) += dtau_wu*grad_phi[i][qp](alpha)*JxW[qp];
                            Kwv(i,j) += dtau_wv*grad_phi[i][qp](alpha)*JxW[qp];
                            Kww(i,j) += dtau_ww*grad_phi[i][qp](alpha)*JxW[qp];
                        }
                    }
                }
            }
        }

    } // qp loop

    // If the Jacobian was requested, we computed it. Otherwise, we didn't.
    return request_jacobian;
}


bool ElasticitySystem::side_time_derivative (bool request_jacobian,
        libMesh::DiffContext & context)
{
    libMesh::FEMContext & c = libMesh::cast_ref<libMesh::FEMContext &>(context);

    // If we're on the correct side, apply the traction
    if (c.has_side_boundary_id(BOUNDARY_ID_MAX_X))
    {
        // If we have an unsteady solver, then we need to extract the corresponding
        // velocity variable. This allows us to use either a FirstOrderUnsteadySolver
        // or a SecondOrderUnsteadySolver. That is, we get back the velocity variable
        // index for FirstOrderUnsteadySolvers or, if it's a SecondOrderUnsteadySolver,
        // this is actually just giving us back the same variable index.

        // If we only wanted to use a SecondOrderUnsteadySolver, then this
        // step would be unnecessary and we would just
        // populate the _u_var, etc. blocks of the residual and Jacobian.
        unsigned int u_dot_var = this->get_second_order_dot_var(_u_var);
        unsigned int v_dot_var = this->get_second_order_dot_var(_v_var);
        unsigned int w_dot_var = this->get_second_order_dot_var(_w_var);

        libMesh::FEBase * u_side_fe;
        c.get_side_fe(_u_var, u_side_fe);

        // The number of local degrees of freedom in each variable
        const unsigned int n_u_dofs = c.get_dof_indices(_u_var).size();

        libMesh::DenseSubVector<libMesh::Number> & Fu = c.get_elem_residual(u_dot_var);
        libMesh::DenseSubVector<libMesh::Number> & Fv = c.get_elem_residual(v_dot_var);
        libMesh::DenseSubVector<libMesh::Number> & Fw = c.get_elem_residual(w_dot_var);

        // Element Jacobian * quadrature weights for interior integration
        const std::vector<libMesh::Real> & JxW = u_side_fe->get_JxW();

        const std::vector<std::vector<libMesh::Real>> & phi = u_side_fe->get_phi();

        unsigned int n_qpoints = c.get_side_qrule().n_points();

        libMesh::Gradient traction(0.0, 0.0, -1.0);

        for (unsigned int qp=0; qp != n_qpoints; qp++)
        {
            for (unsigned int i=0; i != n_u_dofs; i++)
            {
                Fu(i) -= traction(0)*phi[i][qp]*JxW[qp];
                Fv(i) -= traction(1)*phi[i][qp]*JxW[qp];
                Fw(i) -= traction(2)*phi[i][qp]*JxW[qp];
            }
        }
    }

    // If the Jacobian was requested, we computed it (in this case it's zero).
    // Otherwise, we didn't.
    return request_jacobian;
}

bool ElasticitySystem::mass_residual(bool request_jacobian,
        libMesh::DiffContext & context)
{
    libMesh::FEMContext & c = libMesh::cast_ref<libMesh::FEMContext &>(context);

    // We need to extract the corresponding velocity variable.
    // This allows us to use either a FirstOrderUnsteadySolver
    // or a SecondOrderUnsteadySolver. That is, we get back the velocity variable
    // index for FirstOrderUnsteadySolvers or, if it's a SecondOrderUnsteadySolver,
    // this is actually just giving us back the same variable index.

    // If we only wanted to use a SecondOrderUnsteadySolver, then this
    // step would be unnecessary and we would just
    // populate the _u_var, etc. blocks of the residual and Jacobian.
    unsigned int u_dot_var = this->get_second_order_dot_var(_u_var);
    unsigned int v_dot_var = this->get_second_order_dot_var(_v_var);
    unsigned int w_dot_var = this->get_second_order_dot_var(_w_var);

    libMesh::FEBase * u_elem_fe;
    c.get_element_fe(u_dot_var, u_elem_fe);

    // The number of local degrees of freedom in each variable
    const unsigned int n_u_dofs = c.get_dof_indices(u_dot_var).size();

    // Element Jacobian * quadrature weights for interior integration
    const std::vector<libMesh::Real> & JxW = u_elem_fe->get_JxW();

    const std::vector<std::vector<libMesh::Real>> & phi = u_elem_fe->get_phi();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> & Fu = c.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fv = c.get_elem_residual(v_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fw = c.get_elem_residual(w_dot_var);

    libMesh::DenseSubMatrix<libMesh::Number> & Kuu = c.get_elem_jacobian(u_dot_var, u_dot_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvv = c.get_elem_jacobian(v_dot_var, v_dot_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kww = c.get_elem_jacobian(w_dot_var, w_dot_var);

    unsigned int n_qpoints = c.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
        // If we only cared about using FirstOrderUnsteadySolvers for time-stepping,
        // then we could actually just use interior rate, but using interior_accel
        // allows this assembly to function for both FirstOrderUnsteadySolvers
        // and SecondOrderUnsteadySolvers
        libMesh::Number u_ddot, v_ddot, w_ddot;
        c.interior_accel(u_dot_var, qp, u_ddot);
        c.interior_accel(v_dot_var, qp, v_ddot);
        c.interior_accel(w_dot_var, qp, w_ddot);

        for (unsigned int i=0; i != n_u_dofs; i++)
        {
            Fu(i) += _rho*u_ddot*phi[i][qp]*JxW[qp];
            Fv(i) += _rho*v_ddot*phi[i][qp]*JxW[qp];
            Fw(i) += _rho*w_ddot*phi[i][qp]*JxW[qp];

            if (request_jacobian)
            {
                for (unsigned int j=0; j != n_u_dofs; j++)
                {
                    libMesh::Real jac_term = _rho*phi[i][qp]*phi[j][qp]*JxW[qp];
                    jac_term *= context.get_elem_solution_accel_derivative();

                    Kuu(i,j) += jac_term;
                    Kvv(i,j) += jac_term;
                    Kww(i,j) += jac_term;
                }
            }
        }
    }

    // If the Jacobian was requested, we computed it. Otherwise, we didn't.
    return request_jacobian;
}


libMesh::Real ElasticitySystem::elasticity_tensor(unsigned int i, unsigned int j, 
        unsigned int k, unsigned int l)
{
    // Hard code material parameters for the sake of simplicity
    const libMesh::Real poisson_ratio = 0.3;
    const libMesh::Real young_modulus = 1.0e2;

    // Define the Lame constants
    const libMesh::Real lambda_1 = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
    const libMesh::Real lambda_2 = young_modulus/(2.*(1.+poisson_ratio));

    return
        lambda_1 *  kronecker_delta(i, j) * kronecker_delta(k, l) +
        lambda_2 * (kronecker_delta(i, k) * kronecker_delta(j, l) + kronecker_delta(i, l) * kronecker_delta(j, k));
}



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
