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

/// 3D border displacement class, derived from libMesh::FunctionBase<libMesh::Number>.
class border_displacement_function : public libMesh::FunctionBase<libMesh::Number>
{
  private:
    const unsigned int _u_var, _v_var, _w_var;
    const libMesh::Real _x_displ, _y_displ, _z_displ;

  public:
    border_displacement_function (  unsigned int u_var,
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

    virtual void operator() ( const libMesh::Point& p,
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
        (new border_displacement_function(  _u_var, _v_var, _w_var,
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
                                              const std::string& system_name,
                                              libMesh::Order order = libMesh::FIRST,
                                              libMesh::FEFamily family = libMesh::LAGRANGE);

libMesh::NewmarkSystem& add_dynamic_elasticity(libMesh::EquationSystems& input_systems,
//libMesh::LinearImplicitSystem& add_dynamic_elasticity(libMesh::EquationSystems& input_systems,
                                               const std::string& system_name,
                                               libMesh::Order order = libMesh::FIRST,
                                               libMesh::FEFamily family = libMesh::LAGRANGE);

#endif /* COMMON_FUNCTIONS_ELASTICITY_3D_H_ */
