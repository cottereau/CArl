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

#include "weight_parameter_function.h"

#include "assemble_functions_elasticity_3D.h"

class anisotropic_elasticity_tensor_cubic_sym
{
private:
	std::vector<libMesh::DenseMatrix<libMesh::Real> > m_Rotation_matrices;
	std::vector<double>				m_angles_x;
	std::vector<double>				m_angles_y;
	std::vector<double>				m_angles_z;
	std::unordered_map<int,int>		m_domain_to_vec_map;

	double					m_c11;
	double 					m_c12;
	double					m_c44;

	double m_meanE;
	double m_meanMu;

	int 					m_nb_grains;

	void generate_rotation_matrices();

public:
	anisotropic_elasticity_tensor_cubic_sym() :
		m_c11 { -1 },
		m_c12 { -1 },
		m_c44 { -1 },
		m_meanE { -1 },
		m_meanMu { -1 },
		m_nb_grains { -1}
	{

	};

	libMesh::Real eval_internal_elasticity_tensor(unsigned int i,
							  unsigned int j,
							  unsigned int k,
							  unsigned int l);

	anisotropic_elasticity_tensor_cubic_sym(libMesh::EquationSystems& es, std::string& physicalParamsFile,
			double& meanE, double& meanMu)
	{
		this->set_parameters(es, physicalParamsFile,meanE,meanMu);
	};

	void set_parameters(libMesh::EquationSystems& es, std::string& physicalParamsFile,double& meanE, double& meanMu);

	libMesh::DenseMatrix<libMesh::Real>& get_rotation(int idx_grain);

	libMesh::Real eval_elasticity_tensor(unsigned int i,
							  unsigned int j,
							  unsigned int k,
							  unsigned int l,
							  int idx_grain);
};


#endif /* EXECS_ASSEMBLE_FUNCTIONS_ASSEMBLE_FUNCTIONS_ELASTICITY_ANISOTROPY_3D_H_ */
