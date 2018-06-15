// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#ifndef EXECS_ELASTICITY_SYSTEM_H_
#define EXECS_ELASTICITY_SYSTEM_H_

// DiffSystem framework files
#include "libmesh/fem_system.h"

// boundary IDs
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

//  // Element residual and jacobian calculations
//  // Time dependent parts
//  virtual bool element_time_derivative (bool request_jacobian,
//                                        libMesh::DiffContext & context);
//
//  virtual bool side_time_derivative (bool request_jacobian,
//                                     libMesh::DiffContext & context);
//
//  // Mass matrix part
//  virtual bool mass_residual (bool request_jacobian,
//                              libMesh::DiffContext & context);

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

//ElasticitySystem& add_dyn_newmark(libMesh::EquationSystems& input_systems, 
//    libMesh::Order order = libMesh::FIRST, 
//    libMesh::FEFamily family = libMesh::LAGRANGE);
#endif /* EXECS_ELASTICITY_SYSTEM_H_ */


/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=4 et tw=80 smartindent :                               */
