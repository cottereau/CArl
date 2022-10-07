#ifndef ASSEMBLE_FUNCTIONS_MASS_3D_H_
#define ASSEMBLE_FUNCTIONS_MASS_3D_H_

#include "common_header_ext_solver_libmesh.h"
#include "common_assemble_functions_elasticity_3D.h"
#include "weight_parameter_function.h"
#include "ext_solver_libmesh_enums.h"
#include "libmesh_assemble_system_input_parser.h"
#include "newmark_param_parser.h"

void Update_SubM(libMesh::DenseSubMatrix<libMesh::Number>& SubM,
          unsigned int qp,
          const std::vector<std::vector<libMesh::Real>> & phi,
          const unsigned int n_components,
          const unsigned int n_u_dofs,
          const unsigned int n_v_dofs,
          const std::vector<libMesh::Real>& JxW,
          double weight);


void Update_SubC(libMesh::DenseSubMatrix<libMesh::Number>& SubC,
  libMesh::DenseSubMatrix<libMesh::Number>& SubK,
  libMesh::DenseSubMatrix<libMesh::Number>& SubM,
  const unsigned int n_u_dofs,
  const unsigned int n_v_dofs,
  double CK,double CM);

void assemble_dynamic_elasticity_with_weight(libMesh::EquationSystems& es, 
          const std::string& system_name,
          weight_parameter_function& weight_mask,
          WeightFunctionSystemType system_type,
          carl::NewmarkParams* newmark,
          double CM,
          double CK);

void assemble_dynamic_elasticity_with_weight_and_traction(libMesh::EquationSystems& es,
          const std::string& system_name, 
          weight_parameter_function& weight_mask,
          WeightFunctionSystemType system_type,
          int traction_boundary_id,
          std::vector<double> traction_density,
          carl::NewmarkParams* newmark,
          double CM,
          double CK);

void assemble_dynamic_elasticity_with_weight_and_bending(libMesh::EquationSystems& es,
          const std::string& system_name, 
          weight_parameter_function& weight_mask,
          WeightFunctionSystemType system_type,
          int traction_boundary_id,
          std::vector<double> traction_density,
          carl::NewmarkParams* newmark,
          double CM,
          double CK);


void apply_initial(libMesh::EquationSystems & es,
                   const std::string & system_name,
                   const bool & zeroed,
                   const std::string& dis_vec_name,
                   const std::string& vel_vec_name,
                   const std::string& acc_vec_name,
                   libMesh::Parallel::Communicator& WorldComm);


void update_dynamic_rhs(libMesh::EquationSystems& es,
          const std::string& system_name, 
          weight_parameter_function& weight_mask,
          WeightFunctionSystemType system_type,
          double amp_n,
          const std::string& stf_mat_file,
          const std::string& fex_vec_name,
          const std::string& dis_vec_name,
          const std::string& vel_vec_name,
          const std::string& acc_vec_name,
          libMesh::Parallel::Communicator& WorldComm,
          double deltat, 
          double beta);

#endif /*ASSEMBLE_FUNCTIONS_MASS_3D_H_*/
/* Local Variables:                                                        */
/* mode: c                                                                 */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=2 ts=2 et tw=80 smartindent :                               */