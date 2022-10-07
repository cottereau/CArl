/*
 * assemble_functions_elasticity_3D.h
 *
 *  Created on: Nov 30, 2015
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_
#define ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_

#include "common_header_ext_solver_libmesh.h"
#include "common_assemble_functions_elasticity_3D.h"
#include "weight_parameter_function.h"
#include "ext_solver_libmesh_enums.h"



const bool MASTER_bPerfLog_assemble_fem = false;

/// Set the homogeneous physical properties from a file
void set_homogeneous_physical_properties(libMesh::EquationSystems& es, std::string& physicalParamsFile);

/// Set the heterogeneous, isotropic physical properties from a file
void set_heterogeneous_physical_properties(libMesh::EquationSystems& es, std::string& physicalParamsFile);

/// Calculate lambda_1 from E and Mu
libMesh::Real eval_lambda_1(libMesh::Real E, libMesh::Real mu);

/// Calculate the elasticity tensor
libMesh::Real eval_elasticity_tensor(unsigned int i,
						  unsigned int j,
						  unsigned int k,
						  unsigned int l,
						  libMesh::Number E,
						  libMesh::Number mu);

/// Calculate the rigidity sub-matrix contribution
void Update_SubK_isotropic(	libMesh::DenseSubMatrix<libMesh::Number>& SubK,
					unsigned int qp,
					unsigned int C_i,
					unsigned int C_k,
					const std::vector<std::vector<libMesh::RealGradient> >& dphi,
					const unsigned int n_components,
					const unsigned int n_u_dofs,
					const std::vector<libMesh::Real>& JxW,
					libMesh::Number E = 1.0,
					libMesh::Number mu = 0.4,
					double cte = 1
					);

/// Assemble homogeneous elasticity with domain weights
void assemble_elasticity_with_weight(	libMesh::EquationSystems& es,
							const std::string& system_name,
							weight_parameter_function& weight_mask,
							WeightFunctionSystemType system_type);

/// Assemble homogeneous elasticity with domain weights and traction
void assemble_elasticity_with_weight_and_traction(libMesh::EquationSystems& es,
							const std::string& system_name, 
							weight_parameter_function& weight_mask,
							WeightFunctionSystemType system_type,
							int traction_boundary_id,
							std::vector<double> traction_density);

/// Assemble heterogeneous elasticity with domain weights
void assemble_elasticity_heterogeneous_with_weight(	libMesh::EquationSystems& es,
							const std::string& system_name,
							weight_parameter_function& weight_mask,
							WeightFunctionSystemType system_type);

/// Compute the stress (based on one of libMesh's examples)
void compute_stresses(libMesh::EquationSystems& es);

#endif /* ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_H_ */
