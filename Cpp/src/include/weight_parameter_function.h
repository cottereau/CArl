/*
 * weight_parameter_function.h
 *
 *  Created on: Feb 2, 2016
 *      Author: Thiago Milanetto Schlittler
 */

#ifndef COMMON_LIBMESH_CODE_WEIGHT_PARAMETER_FUNCTION_H_
#define COMMON_LIBMESH_CODE_WEIGHT_PARAMETER_FUNCTION_H_

#include "carl_headers.h"

namespace carl
{

class weight_parameter_function
{
protected:
	// Members
//	const libMesh::PointLocatorBase& m_locator;
	libMesh::Mesh& m_alpha_mesh;
	std::unique_ptr<libMesh::PointLocatorBase> m_locator_unique_ptr;

	double m_alpha_eps;
	double m_alpha_coupling_BIG;
	double m_alpha_coupling_micro;

	int m_subdomain_idx_BIG;
	int m_subdomain_idx_micro;
	int m_subdomain_idx_coupling;

public:

	// Constructor
	weight_parameter_function(	libMesh::Mesh& alpha_mesh, double alpha_eps,
						double alpha_coupling_BIG,
						int subdomain_idx_BIG, int subdomain_idx_micro, int subdomain_idx_coupling) :
		m_alpha_mesh{ alpha_mesh },
		m_locator_unique_ptr { alpha_mesh.sub_point_locator() },
		m_alpha_eps { alpha_eps },
		m_alpha_coupling_BIG { alpha_coupling_BIG },
		m_alpha_coupling_micro { 1. - alpha_coupling_BIG },
		m_subdomain_idx_BIG { subdomain_idx_BIG },
		m_subdomain_idx_micro { subdomain_idx_micro },
		m_subdomain_idx_coupling { subdomain_idx_coupling }
	{
	};

	weight_parameter_function(	libMesh::Mesh& alpha_mesh ) :
		m_alpha_mesh{ alpha_mesh },
		m_locator_unique_ptr { alpha_mesh.sub_point_locator() },
		m_alpha_eps { 10E-2 },
		m_alpha_coupling_BIG { 0.5 },
		m_alpha_coupling_micro { 0.5 },
		m_subdomain_idx_BIG { -1 },
		m_subdomain_idx_micro { -1 },
		m_subdomain_idx_coupling { -1 }
	{
	};

	// Destructor
	~weight_parameter_function()
	{
		clear();
	}

	void set_parameters(double alpha_eps, double alpha_coupling_BIG,
						int subdomain_idx_BIG, int subdomain_idx_micro, int subdomain_idx_coupling)
	{
		m_alpha_eps = alpha_eps;
		m_alpha_coupling_BIG = alpha_coupling_BIG;
		m_alpha_coupling_micro = 1. - alpha_coupling_BIG;
		m_subdomain_idx_BIG  = subdomain_idx_BIG;
		m_subdomain_idx_micro = subdomain_idx_micro;
		m_subdomain_idx_coupling = subdomain_idx_coupling;
	};

	double get_alpha_BIG(const libMesh::Point& qpoint)
	{
		libMesh::PointLocatorBase& locator = *m_locator_unique_ptr.get();
		const libMesh::Elem* elem = locator(qpoint);
		double output = 0;

		if(elem->subdomain_id() == m_subdomain_idx_BIG)
		{
			output = 1;
		}
		else if(elem->subdomain_id() == m_subdomain_idx_micro)
		{
			output = m_alpha_eps;
		}
		else if(elem->subdomain_id() == m_subdomain_idx_coupling)
		{
			output = m_alpha_coupling_BIG;
		}

		return output;
	}

	double get_alpha_micro(const libMesh::Point& qpoint)
	{
		libMesh::PointLocatorBase& locator = *m_locator_unique_ptr.get();
		const libMesh::Elem* elem = locator(qpoint);
		double output = 0;

		if(elem->subdomain_id() == m_subdomain_idx_BIG)
		{
			output = 0;
		}
		else if(elem->subdomain_id() == m_subdomain_idx_micro)
		{
			output = 1 - m_alpha_eps;
		}
		else if(elem->subdomain_id() == m_subdomain_idx_coupling)
		{
			output = m_alpha_coupling_micro;
		}

		return output;
	}

	// Methods
	void clear()
	{
		m_locator_unique_ptr.reset(NULL);
		// Nothing to do ... for now
	};
};

}

#endif /* COMMON_LIBMESH_CODE_WEIGHT_PARAMETER_FUNCTION_H_ */
