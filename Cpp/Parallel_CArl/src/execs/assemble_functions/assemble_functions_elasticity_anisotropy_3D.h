/*
 * assemble_functions_elasticity_anisotropy_3D.h
 *
 *  Created on: Sep 3, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef EXECS_ASSEMBLE_FUNCTIONS_ASSEMBLE_FUNCTIONS_ELASTICITY_ANISOTROPY_3D_H_
#define EXECS_ASSEMBLE_FUNCTIONS_ASSEMBLE_FUNCTIONS_ELASTICITY_ANISOTROPY_3D_H_

#include "common_header_libmesh.h"
#include "common_functions.h"
#include "anisotropic_elasticity_cubic_sym.h"
#include "weight_parameter_function.h"

#include "assemble_functions_elasticity_3D.h"

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
					);

void assemble_elasticity_anisotropic_with_weight(	libMesh::EquationSystems& es,
							const std::string& system_name,
							carl::weight_parameter_function& weight_mask,
							carl::anisotropic_elasticity_tensor_cubic_sym& anisotropy_obj_input);

#endif /* EXECS_ASSEMBLE_FUNCTIONS_ASSEMBLE_FUNCTIONS_ELASTICITY_ANISOTROPY_3D_H_ */
