#ifndef ASSEMBLE_FUNCTIONS_MASS_3D_H_
#define ASSEMBLE_FUNCTIONS_MASS_3D_H_

#include "common_header_ext_solver_libmesh.h"
#include "common_assemble_functions_elasticity_3D.h"
#include "weight_parameter_function.h"
#include "ext_solver_libmesh_enums.h"


struct libmesh_assemble_interpolation_params
{
  std::string path_tilde_matrices_A;
  std::string path_macro_coupling_matrix;
  libMesh::Real alpha_A;

  std::string path_tilde_matrices_B;
  std::string path_micro_coupling_matrix;
  libMesh::Real alpha_B;

  std::string output_folder;
};

void get_input_params(GetPot& field_parser,
     libmesh_assemble_interpolation_params& input_params);


libMesh::SparseMatrix< libMesh::Number > * get_stiffness_matrix(libMesh::EquationSystems& es,
				const std::string& system_name,
				weight_parameter_function& weight_mask,
				WeightFunctionSystemType system_type);

libMesh::SparseMatrix< libMesh::Number> * get_mass_matrix(libMesh::EquationSystems& es,
				const std::string& system_name);

void Update_SubM( libMesh::DenseSubMatrix<libMesh::Number>& SubM,
          unsigned int qp,
          const std::vector<std::vector<libMesh::Real>> & phi,
          const unsigned int n_components,
          const unsigned int n_u_dofs,
          const std::vector<libMesh::Real>& JxW,
          double weight);

libMesh::SparseMatrix< libMesh::Number > *  assemble_mass_matrix(libMesh::EquationSystems& es,
				const std::string& system_name);



#endif /*ASSEMBLE_FUNCTIONS_MASS_3D_H_*/