/*
 * common_functions_elasticity_3D.h
 *
 *  Created on: Apr 17, 2017
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_FUNCTIONS_ELASTICITY_3D_H_
#define COMMON_FUNCTIONS_ELASTICITY_3D_H_

#include "common_header_ext_solver_libmesh.h"
#include "ext_solver_libmesh_enums.h"
#include "libmesh/fem_system.h"
/// 3D border displacement class, derived from libMesh::FunctionBase<libMesh::Number>.

#define BOUNDARY_ID_MIN_Z 0
#define BOUNDARY_ID_MIN_Y 1
#define BOUNDARY_ID_MAX_X 2
#define BOUNDARY_ID_MAX_Y 3
#define BOUNDARY_ID_MIN_X 4
#define BOUNDARY_ID_MAX_Z 5
#define NODE_BOUNDARY_ID 10
#define EDGE_BOUNDARY_ID 20

//using namespace libMesh;

// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class ElasticitySystem : public libMesh::FEMSystem
{
    public:
        // Constructor
        ElasticitySystem(libMesh::EquationSystems & es,
                const std::string & name_in,
                const unsigned int number_in)
            : libMesh::FEMSystem(es, name_in, number_in),
            _rho(1.0)
    {}

        // System initialization
        virtual void init_data ();

        // Context initialization
        virtual void init_context(libMesh::DiffContext & context);

        // Element residual and jacobian calculations
        // Time dependent parts
        virtual bool element_time_derivative (bool request_jacobian,
                libMesh::DiffContext & context);
        //
        virtual bool side_time_derivative (bool request_jacobian,
                libMesh::DiffContext & context);
        //
        // Mass matrix part
        virtual bool mass_residual (bool request_jacobian,
                libMesh::DiffContext & context);

    private:

        // Indices for each variable;
        unsigned int _u_var, _v_var, _w_var;

        libMesh::Real _rho;

        libMesh::Real kronecker_delta(unsigned int i, unsigned int j)
        {
            return i == j ? 1. : 0.;
        }

        libMesh::Real elasticity_tensor(unsigned int i, unsigned int j, unsigned int k, unsigned int l);
};



class border_displacement_function : public libMesh::FunctionBase<libMesh::Number>
{
    private:
        const unsigned int _u_var, _v_var, _w_var;
        const libMesh::Real _x_displ, _y_displ, _z_displ;

    public:
        border_displacement_function (	unsigned int u_var,
                unsigned int v_var,
                unsigned int w_var,
                libMesh::Real x_displ = 0,
                libMesh::Real y_displ = 0,
                libMesh::Real z_displ = 0)
            : _u_var(u_var),
            _v_var(v_var),
            _w_var(w_var),
            _x_displ(x_displ),
            _y_displ(y_displ),
            _z_displ(z_displ)
    {
        this->_initialized = true;
    }

        virtual libMesh::Number operator() (const libMesh::Point&, const libMesh::Real = 0)
        {
            libmesh_not_implemented();
        }

        virtual void operator() (	const libMesh::Point& p,
                const libMesh::Real,
                libMesh::DenseVector<libMesh::Number>& output)
        {
            output.resize(3);
            output.zero();

            output(_u_var) = _x_displ;
            output(_v_var) = _y_displ;
            output(_w_var) = _z_displ;
        }

        virtual libMesh::UniquePtr<FunctionBase<libMesh::Number> > clone() const
        {
            return libMesh::UniquePtr<FunctionBase<libMesh::Number> >
                (new border_displacement_function(	_u_var, _v_var, _w_var,
                                                    _x_displ, _y_displ, _z_displ));
        }
};

/// Small structure with the 3D border displacement values.
struct border_displacement_values
{
    double x_displ;
    double y_displ;
    double z_displ;
};

/// Set a displacement border.
void set_displaced_border_translation(libMesh::ImplicitSystem& elasticity_system, border_displacement_values& displ, int boundary_id);

/// Set a clamped border.
void set_clamped_border(libMesh::ImplicitSystem& elasticity_system, int boundary_id);

/// Add a stress libMesh::ExplicitSystem to the input libMesh::EquationSystems.
libMesh::ExplicitSystem& add_stress(libMesh::EquationSystems& input_systems);

/// Add a linear elasticity libMesh::LinearImplicitSystem to the input libMesh::EquationSystems& input_systems.
libMesh::LinearImplicitSystem& add_elasticity(libMesh::EquationSystems& input_systems,
        libMesh::Order order = libMesh::FIRST,
        libMesh::FEFamily family = libMesh::LAGRANGE);

libMesh::ExplicitSystem& add_vel_newmark(libMesh::EquationSystems& input_systems, 
        libMesh::Order order = libMesh::FIRST, 
        libMesh::FEFamily family = libMesh::LAGRANGE);

libMesh::ExplicitSystem& add_acc_newmark(libMesh::EquationSystems& input_systems, 
        libMesh::Order order = libMesh::FIRST, 
        libMesh::FEFamily family = libMesh::LAGRANGE);

#endif /* COMMON_FUNCTIONS_ELASTICITY_3D_H_ */

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=4 et tw=80 smartindent :                               */
