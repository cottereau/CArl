/*
 * assemble_functions_elasticity_3D_dyn.h
 *
 *  Created on: May 28, 2018
 *      Author: Filippo Gatti
 */

#ifndef ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_DYN_H_
#define ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_DYN_H_

#include "common_header_ext_solver_libmesh.h"

#include "common_assemble_functions_elasticity_3D.h"
#include "weight_parameter_function.h"
#include "ext_solver_libmesh_enums.h"

/// Set the homogeneous physical properties from a file ---> add rho
void set_homogeneous_physical_properties_dyn(libMesh::EquationSystems& es, std::string& physicalParamsFile);

/// Set the heterogeneous, isotropic physical properties from a file ---> add rho
void set_heterogeneous_physical_properties_dyn(libMesh::EquationSystems& es, std::string& physicalParamsFile);

void assemble_elasticity_with_weight_and_traction_dyn(libMesh::EquationSystems& es,\
    const std::string& system_name,weight_parameter_function& weight_mask,\
    WeightFunctionSystemType system_type,int traction_boundary_id,\
    std::vector<double> traction_density);

#endif /* ELASTICITY_3D_ASSEMBLE_FUNCTIONS_ELASTICITY_3D_DYN_H_ */
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=4 et tw=80 smartindent :                               */
