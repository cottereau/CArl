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

void carl::anisotropic_elasticity_tensor_cubic_sym::set_parameters(libMesh::EquationSystems& es, std::string& physicalParamsFile,
		double& meanE, double& meanMu)
{
	const libMesh::Parallel::Communicator& SysComm = es.comm();
	int rank = SysComm.rank();
	int nodes = SysComm.size();
	std::vector<int> temp_idx;

	// Read the random data info
	if(rank == 0)
	{
		std::ifstream physicalParamsIFS(physicalParamsFile);
		physicalParamsIFS >> m_nb_grains >> m_c11 >> m_c12 >> m_c44 >> m_meanE >> m_meanMu;
		std::cout << m_nb_grains << m_c11 << m_c12 << m_c44 <<  m_meanE << m_meanMu << std::endl;
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
		SysComm.broadcast(m_meanE);
		SysComm.broadcast(m_meanMu);

		meanE = m_meanE;
		meanMu = m_meanMu;

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

		for(int iii = 0; iii < m_nb_grains; ++iii)
		{
			m_domain_to_vec_map[temp_idx[iii]] = iii;
		}

		m_Rotation_matrices.resize(m_nb_grains,libMesh::DenseMatrix<libMesh::Real>(3,3));
	}

	// Mesh pointer
	const libMesh::MeshBase& mesh = es.get_mesh();

	// Indexes
	libMesh::ExplicitSystem& physical_param_system = es.get_system<libMesh::ExplicitSystem>("PhysicalConstants");
	const libMesh::DofMap& physical_dof_map = physical_param_system.get_dof_map();

	unsigned int physical_consts[1];
	physical_consts[0] = physical_param_system.variable_number ("Index");

	std::vector<libMesh::dof_id_type> physical_dof_indices_var;
	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	int currentSubdomain = -1;

	for ( ; el != end_el; ++el)
	{
		const libMesh::Elem* elem = *el;

		currentSubdomain = elem->subdomain_id()-1;

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

	// Calculate the rotation matrices
	this->generate_rotation_matrices();
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

		m_Rotation_matrices[iii](0,0) =
				c_theta * c_psi;
		m_Rotation_matrices[iii](0,1) =
				- c_phi * s_psi - s_phi * s_theta * c_psi;
		m_Rotation_matrices[iii](0,2) =
				s_phi * s_psi 	- c_phi * s_theta * c_psi;

		m_Rotation_matrices[iii](1,0) =
				c_theta * s_psi;
		m_Rotation_matrices[iii](1,1) =
				c_phi * c_psi 	- s_phi * s_theta * s_psi;
		m_Rotation_matrices[iii](1,2) =
				- s_phi * c_psi - c_phi * s_theta * s_psi;

		m_Rotation_matrices[iii](2,0) =
				- s_theta;
		m_Rotation_matrices[iii](2,1) =
				s_phi * c_theta;
		m_Rotation_matrices[iii](2,2) =
				c_phi * c_theta;
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
	else if(i == k && j == l)
	{
		output = m_c44;
	}

	return output;
}

libMesh::DenseMatrix<libMesh::Real>&  carl::anisotropic_elasticity_tensor_cubic_sym::get_rotation(int idx_grain)
{
	return m_Rotation_matrices[m_domain_to_vec_map[idx_grain]];
}

libMesh::Real  carl::anisotropic_elasticity_tensor_cubic_sym::eval_elasticity_tensor(unsigned int i,
						  unsigned int j,
						  unsigned int k,
						  unsigned int l,
						  int idx_grain)
{
	double output = 0;
	libMesh::DenseMatrix<libMesh::Real>& R = get_rotation(idx_grain);

	for(unsigned int ppp = 0; ppp < 3; ++ppp)
	{
		for(unsigned int qqq = 0; qqq < 3; ++qqq)
		{
			for(unsigned int rrr = 0; rrr < 3; ++rrr)
			{
				for(unsigned int sss = 0; sss < 3; ++sss)
				{
					output += R(i,ppp) * R(j,qqq) * R(k,rrr) * R(l,sss) *
							this->eval_internal_elasticity_tensor(ppp,qqq,rrr,sss);
				}
			}
		}
	}

	return output;
}



