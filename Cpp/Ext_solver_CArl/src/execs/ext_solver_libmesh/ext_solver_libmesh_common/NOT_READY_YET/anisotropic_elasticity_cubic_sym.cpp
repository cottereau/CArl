/*
 * anisotropic_elasticity_cubic_sym.cpp
 *
 *  Created on: Sep 4, 2016
 *      Author: Thiago Milanetto Schlittler
 */

/*
 * assemble_functions_elasticity_anisotropy.cpp
 *
 *  Created on: Sep 3, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#include "anisotropic_elasticity_cubic_sym.h"

void carl::anisotropic_elasticity_tensor_cubic_sym::set_parameters(libMesh::EquationSystems& es, std::string& physicalParamsFile)
{
	const libMesh::Parallel::Communicator& SysComm = es.comm();
	int rank = SysComm.rank();
	int nodes = SysComm.size();
	std::vector<int> temp_idx;

	// Read the random data info
	if(rank == 0)
	{
		std::ifstream physicalParamsIFS(physicalParamsFile);
		physicalParamsIFS >> m_nb_grains >> m_c11 >> m_c12 >> m_c44;
		m_angles_x.resize(m_nb_grains);
		m_angles_y.resize(m_nb_grains);
		m_angles_z.resize(m_nb_grains);
		temp_idx.resize(m_nb_grains);
		m_domain_to_vec_map.reserve(m_nb_grains);

		for(int iii = 0; iii < m_nb_grains; ++iii)
		{
			physicalParamsIFS >> m_angles_x[iii];
			physicalParamsIFS >> m_angles_y[iii];
			physicalParamsIFS >> m_angles_z[iii];
			physicalParamsIFS >> temp_idx[iii];
		}

		physicalParamsIFS.close();
	}

	if(nodes > 1)
	{
		SysComm.broadcast(m_c11);
		SysComm.broadcast(m_c12);
		SysComm.broadcast(m_c44);
		SysComm.broadcast(m_nb_grains);

		if(rank != 0)
		{
			m_angles_x.resize(m_nb_grains);
			m_angles_y.resize(m_nb_grains);
			m_angles_z.resize(m_nb_grains);
			temp_idx.resize(m_nb_grains);
			m_domain_to_vec_map.reserve(m_nb_grains);
		}
		SysComm.broadcast(m_angles_x);
		SysComm.broadcast(m_angles_y);
		SysComm.broadcast(m_angles_z);
		SysComm.broadcast(temp_idx);
	}

	for(int iii = 0; iii < m_nb_grains; ++iii)
	{
		m_domain_to_vec_map[temp_idx[iii]] = iii;
	}

	m_Rotation_matrices.resize(m_nb_grains,libMesh::DenseMatrix<libMesh::Real>(3,3));
	m_Elasticity_tensors.resize(m_nb_grains);
	for(int iii = 0; iii < m_nb_grains; ++iii)
	{
		m_Elasticity_tensors[iii].set_dimension(3);
	}

	// Mesh pointer
	const libMesh::MeshBase& mesh = es.get_mesh();

	// Indexes
	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();

	unsigned int physical_consts[7];
	physical_consts[0] = physical_param_system.variable_number ("Index");
	physical_consts[1] = physical_param_system.variable_number ("Angle_x");
	physical_consts[2] = physical_param_system.variable_number ("Angle_y");
	physical_consts[3] = physical_param_system.variable_number ("Angle_z");
	physical_consts[4] = physical_param_system.variable_number ("color_r");
	physical_consts[5] = physical_param_system.variable_number ("color_g");
	physical_consts[6] = physical_param_system.variable_number ("color_b");

	// Calculate the rotation matrices
	this->generate_rotation_matrices();

	// And generate the elasticity tensors
	this->generate_elasticity_tensors();

	std::vector<libMesh::dof_id_type> physical_dof_indices_var;
	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	for ( ; el != end_el; ++el)
	{
		const libMesh::Elem* elem = *el;

		// Grain index
		physical_dof_map.dof_indices(elem, physical_dof_indices_var, physical_consts[0]);
		libMesh::dof_id_type dof_index = physical_dof_indices_var[0];

		if( (physical_param_system.solution->first_local_index() <= dof_index) &&
		(dof_index < physical_param_system.solution->last_local_index()) )
		{
			physical_param_system.solution->set(dof_index, elem->subdomain_id());
		}
	}

	physical_param_system.solution->close();
	physical_param_system.update();
}

void  carl::anisotropic_elasticity_tensor_cubic_sym::generate_rotation_matrices()
{
	double phi, theta, psi;
	double c_phi, c_theta, c_psi;
	double s_phi, s_theta, s_psi;
	for(int iii = 0; iii < m_nb_grains; ++iii)
	{
		phi = m_angles_x[iii];
		theta = m_angles_y[iii];
		psi = m_angles_z[iii];

		s_phi 	= sin(phi); 	c_phi 	= cos(phi);
		s_theta = sin(theta); 	c_theta = cos(theta);
		s_psi 	= sin(psi); 	c_psi 	= cos(psi);

		m_Rotation_matrices[iii](0,0) =           c_theta;
		m_Rotation_matrices[iii](0,1) =         - s_theta * c_psi;
		m_Rotation_matrices[iii](0,2) =           s_theta * s_psi;

		m_Rotation_matrices[iii](1,0) =   c_phi * s_theta;
		m_Rotation_matrices[iii](1,1) =   c_phi * c_theta * c_psi
				                        - s_phi           * s_psi;
		m_Rotation_matrices[iii](1,2) = - s_phi           * c_psi
				                        - c_phi * c_theta * s_psi;

		m_Rotation_matrices[iii](2,0) =   s_phi * s_theta;
		m_Rotation_matrices[iii](2,1) =   c_phi           * s_psi
				                        + s_phi * c_theta * c_psi;
		m_Rotation_matrices[iii](2,2) =   c_phi           * c_psi
				                        - s_phi * c_theta * s_psi;
	}
}

void  carl::anisotropic_elasticity_tensor_cubic_sym::generate_elasticity_tensors()
{
	for(int nnn = 0; nnn < m_nb_grains; ++nnn)
	{
		libMesh::DenseMatrix<libMesh::Real>& R = m_Rotation_matrices[nnn];
		Order4Tensor& C = 						 m_Elasticity_tensors[nnn];

		for(unsigned int iii = 0; iii < 3; ++iii)
		for(unsigned int jjj = 0; jjj < 3; ++jjj)
		for(unsigned int kkk = 0; kkk < 3; ++kkk)
		for(unsigned int lll = 0; lll < 3; ++lll)
		{
			// Contraction indexes
			for(unsigned int ppp = 0; ppp < 3; ++ppp)
			for(unsigned int qqq = 0; qqq < 3; ++qqq)
			for(unsigned int rrr = 0; rrr < 3; ++rrr)
			for(unsigned int sss = 0; sss < 3; ++sss)
			{
				 C(iii,jjj,kkk,lll)+= R(iii,ppp) * R(jjj,qqq) * R(kkk,rrr) * R(lll,sss) *
						this->eval_internal_elasticity_tensor(ppp,qqq,rrr,sss);
			}
		}
	}
}

libMesh::Real  carl::anisotropic_elasticity_tensor_cubic_sym::eval_internal_elasticity_tensor(unsigned int i,
		  unsigned int j,
		  unsigned int k,
		  unsigned int l)
{
	double output = 0;

	if(i == j && k == l)
	{
		if(i == k && j == l)
		{
			output = m_c11;
		}
		else
		{
			output = m_c12;
		}
	}
	else if( (i == k && j == l) || (i == l && j == k) )
	{
		output = m_c44;
	}

	return output;
}

libMesh::DenseMatrix<libMesh::Real>&  carl::anisotropic_elasticity_tensor_cubic_sym::get_rotation(int idx_grain)
{
	return m_Rotation_matrices[m_domain_to_vec_map[idx_grain]];
}

carl::Order4Tensor&  carl::anisotropic_elasticity_tensor_cubic_sym::get_elasticity(int idx_grain)
{
	return m_Elasticity_tensors[m_domain_to_vec_map[idx_grain]];
}

libMesh::Real  carl::anisotropic_elasticity_tensor_cubic_sym::eval_elasticity_tensor(unsigned int i,
						  unsigned int j,
						  unsigned int k,
						  unsigned int l,
						  int idx_grain)
{
	Order4Tensor& C = get_elasticity(idx_grain);

	return C(i,j,k,l);
}



