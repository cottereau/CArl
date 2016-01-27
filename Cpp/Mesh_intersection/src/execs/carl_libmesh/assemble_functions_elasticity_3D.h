/*
 * assemble_functions_elasticity_3D.h
 *
 *  Created on: Nov 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_
#define ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_

#include "common_header_libmesh.h"
#include "common_functions.h"

void Update_SubK(	libMesh::DenseSubMatrix<libMesh::Number>& SubK,
					unsigned int qp,
					unsigned int C_i,
					unsigned int C_k,
					const std::vector<std::vector<libMesh::RealGradient> >& dphi,

					const unsigned int n_components,
					const unsigned int n_u_dofs,
					const std::vector<libMesh::Real>& JxW,
					libMesh::Number E = 1.0,
					libMesh::Number mu = 0.4
					);

void set_physical_properties(libMesh::EquationSystems& es, std::string& physicalParamsFile, double& meanE, double& meanMu);

void set_constant_physical_properties(libMesh::EquationSystems& es, double meanE, double meanMu);

void assemble_elasticity(libMesh::EquationSystems& es,
					   const std::string& system_name);

void assemble_elasticity_heterogeneous(libMesh::EquationSystems& es,
					   const std::string& system_name);

void compute_stresses(libMesh::EquationSystems& es);

libMesh::Real eval_elasticity_tensor(unsigned int i,
						  unsigned int j,
						  unsigned int k,
						  unsigned int l,
						  libMesh::Number E = 1.,
						  libMesh::Number mu = 0.4);

#endif /* ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_ */
