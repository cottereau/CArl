#ifndef ASSEMBLE_FUNCTIONS_MASS_3D_H_
#define ASSEMBLE_FUNCTIONS_MASS_3D_H_

#include "common_header_ext_solver_libmesh.h"
#include "common_assemble_functions_elasticity_3D.h"
#include "weight_parameter_function.h"
#include "ext_solver_libmesh_enums.h"
#include "libmesh_assemble_system_input_parser.h"

void Update_SubM(libMesh::DenseSubMatrix<libMesh::Number>& SubM,
          unsigned int qp,
          const std::vector<std::vector<libMesh::Real>> & phi,
          const unsigned int n_components,
          const unsigned int n_u_dofs,
          const std::vector<libMesh::Real>& JxW,
          double weight);

void assemble_dynamic_elasticity_with_weight(libMesh::EquationSystems& es, 
          const std::string& system_name,
          weight_parameter_function& weight_mask,
          WeightFunctionSystemType system_type,
          double deltat, double beta);

void assemble_dynamic_elasticity_with_weight_and_traction(libMesh::EquationSystems& es,
          const std::string& system_name, 
          weight_parameter_function& weight_mask,
          WeightFunctionSystemType system_type,
          int traction_boundary_id,
          std::vector<double> traction_density,
          double deltat, double beta);

void update_dynamic_rhs(libMesh::EquationSystems& es,
          const std::string& system_name, 
          weight_parameter_function& weight_mask,
          WeightFunctionSystemType system_type,
          int traction_boundary_id, double amp_n,
          std::vector<double> traction_density,
          const std::string& stiffness_file,
          const std::string& force_file,
          const std::string& displacement_file,
          const std::string& velocity_file,
          const std::string& acceleration_file,
          double deltat, double beta);

#endif /*ASSEMBLE_FUNCTIONS_MASS_3D_H_*/
