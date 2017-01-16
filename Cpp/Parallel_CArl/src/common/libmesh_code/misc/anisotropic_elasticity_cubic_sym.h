/*
 * anisotropic_elasticity_cubic_sym.h
 *
 *  Created on: Sep 4, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_ANISOTROPIC_ELASTICITY_CUBIC_SYM_H_
#define COMMON_LIBMESH_CODE_ANISOTROPIC_ELASTICITY_CUBIC_SYM_H_

#include "carl_headers.h"

#include "weight_parameter_function.h"

namespace carl
{
class Order4Tensor
{
private:
	int m_dim;
	std::vector<int> m_idx_mult;
	unsigned int m_idx;

public:
	Order4Tensor() :
		m_dim { 0 },
		m_idx { 0 }
	{

	};

	Order4Tensor(int dim) :
		m_idx { 0 }
	{
		m_dim = dim;
		set_dimension(dim);
	};

	std::vector<libMesh::Real> vector;

	void set_dimension(int dim)
	{
		m_dim = dim;
		m_idx_mult.resize(4);
		vector.resize(std::pow(m_dim,4),0);

		for(int iii = 0; iii < 4; ++iii)
		{
			m_idx_mult[iii] = std::pow(m_dim,iii);
		}
	}
	libMesh::Real& operator () (unsigned int iii, unsigned int jjj, unsigned int kkk, unsigned int lll);
};

inline libMesh::Real& carl::Order4Tensor::operator () (unsigned int iii, unsigned int jjj, unsigned int kkk, unsigned int lll)
{
	m_idx =   m_idx_mult[0] * iii
			+ m_idx_mult[1] * jjj
			+ m_idx_mult[2] * kkk
			+ m_idx_mult[3] * lll;
	return vector[m_idx];
}

class anisotropic_elasticity_tensor_cubic_sym
{
private:
	std::vector<libMesh::DenseMatrix<libMesh::Real> > m_Rotation_matrices;
	std::vector<Order4Tensor>          			  m_Elasticity_tensors;
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
	void generate_elasticity_tensors();

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

	anisotropic_elasticity_tensor_cubic_sym(libMesh::EquationSystems& es, std::string& physicalParamsFile,
			double& meanE, double& meanMu)
	{
		this->set_parameters(es, physicalParamsFile,meanE,meanMu);
	};

	libMesh::Real eval_internal_elasticity_tensor(unsigned int i,
							  unsigned int j,
							  unsigned int k,
							  unsigned int l);

	void set_parameters(libMesh::EquationSystems& es, std::string& physicalParamsFile,double& meanE, double& meanMu);

	libMesh::DenseMatrix<libMesh::Real>& get_rotation(int idx_grain);
	Order4Tensor&  get_elasticity(int idx_grain);

	libMesh::Real eval_elasticity_tensor(unsigned int i,
							  unsigned int j,
							  unsigned int k,
							  unsigned int l,
							  int idx_grain);
};
}


#endif /* COMMON_LIBMESH_CODE_ANISOTROPIC_ELASTICITY_CUBIC_SYM_H_ */
